from bisect import bisect_left
from os import makedirs
from os.path import exists
from subprocess import check_call
import pysam as ps
from sklearn.externals.joblib import Parallel, delayed

from core.annotations import path_to_annotations
from core.constants import data_path, dr_genes, upstream_length
from core.data_reading import path_to_ref, read_h37rv
from data_processing.parse_gatk_variants import parse_complex_event

path_to_ids = data_path + 'all_with_pheno.txt'
bam_path = data_path + 'skesa_minimap_mapq/'
# out_path_vcf = data_path + 'snps/skesa_bwa_mapq/'
out_path_vcf = data_path + 'snps/skesa_minimap_mapq/'
# out_path_variants = data_path + 'snps/skesa_bwa_mapq_raw/'
out_path_variants = data_path + 'snps/skesa_minimap_mapq_raw/'

path_to_freebayes = '/export/home/fedonin/freebayes-v1.3.1'

filter_annotation_repeats = False
filter_short_repeats = False
path_to_short_tandem_repeats = data_path + 'h37rv.fasta.2.7.7.80.10.20.10.dat'
filter_out_DR_genes = False


def bcftools(sample_id):
    if exists(bam_path + sample_id + '.bam'):
        check_call('bcftools mpileup -Ou -f ' + path_to_ref + ' ' + bam_path + sample_id +
                   '.bam | bcftools call -c --ploidy 1 -v -o ' + out_path_vcf + sample_id + '.vcf', shell=True)
        return 1
    else:
        return 0


def freebayes(sample_id):
    if exists(bam_path + sample_id + '.bam'):
        check_call(path_to_freebayes + ' --bam ' + bam_path + sample_id + '.bam --vcf ' + out_path_vcf +
                  sample_id + '.vcf --fasta-reference ' + path_to_ref + ' --pvar 0.0001 --ploidy 1 ' +
                  '--min-mapping-quality 0 --min-base-quality 0', shell=True)
        return 1
    else:
        return 0


caller = bcftools


def run_variant_calling():
    if not exists(out_path_vcf):
        makedirs(out_path_vcf)
    tasks = Parallel(n_jobs=-1)(delayed(caller)(sample_id.strip()) for sample_id in open(path_to_ids).readlines())
    c = 0
    for task in tasks:
        c += task
    print('%d bam files filtered\n' % c)


def get_intervals_to_filter_out():
    """
    Reads annotations and picks coordinates of repeats and mobile elements.

    :param filter_out_DR_genes: if True adds coordinates of DR genes from list of DR genes to the results
    :param filter_short_repeats: if True adds coordinates of short tandem repeats to the result
    :return: list of tuples (begin, end), representing intervals
    """
    coords = []
    with open(path_to_annotations, 'r') as f:
        for line in f.readlines():
            s = line.strip().split('\t')
            if filter_out_DR_genes and s[2] == 'Gene':
                strand = s[6]
                for gene_name in dr_genes:
                    if gene_name in s[-1]:
                        if strand == '+':
                            coords.append((int(s[3]) - upstream_length, int(s[4])))
                        else:
                            coords.append((int(s[3]), int(s[4]) + upstream_length))
            elif filter_annotation_repeats and s[2] == 'Repeat_region':
                coords.append((int(s[3]), int(s[4])))
            elif s[2] == 'mobile_element':
                coords.append((int(s[3]), int(s[4])))
    if filter_short_repeats:
        with open(path_to_short_tandem_repeats) as f:
            for line in f.readlines()[15:]:
                s = line.split()
                coords.append((int(s[0]), int(s[1]) + 1))
    coords.sort(key=lambda tup: tup[0])
    return coords


def in_filter_interval(pos, filter_intervals_starts, filter_intervals):
    i = bisect_left(filter_intervals_starts, pos)
    if i == 0:
        if pos == filter_intervals_starts[0]:
            return True
    else:
        if i < len(filter_intervals) and filter_intervals_starts[i] == pos:
            return True
        elif pos <= filter_intervals[i - 1][1]:
            return True
    return False


def parse(sample_name, filter_intervals):
    filter_intervals_starts = [x[0] for x in filter_intervals]
    with open(out_path_variants + sample_name + '.variants', 'w') as fout:
        fout.write('# source ' + out_path_vcf + sample_name + '.vcf\n')
        res = set()
        for l in open(out_path_vcf + sample_name + '.vcf').readlines():
            if l.startswith('#'):
                continue
            s = l.strip().split('\t')
            pos = int(s[1])
            ref = s[3]
            alt = s[4]
            if not in_filter_interval(pos, filter_intervals_starts, filter_intervals):
                if len(ref) != 1 or len(alt) != 1 or ref == '*' or alt == '*':
                    res.update(parse_complex_event(pos, ref, alt))
                else:
                    res.add((pos, alt, 'snp'))
        var_list = list(res)
        var_list.sort(key=lambda x: x[0])
        for pos, alt, t in var_list:
            fout.write('%d\t%s\t%s\n' % (pos, alt, t))
    return 1


def parse_variants():
    if not exists(out_path_variants):
        makedirs(out_path_variants)
    filter_intervals = get_intervals_to_filter_out()
    tasks = Parallel(n_jobs=-1)(delayed(parse)(sample_id.strip(), filter_intervals)
                                for sample_id in open(path_to_ids).readlines())
    c = 0
    for task in tasks:
        c += task
    print('%d samples processed' % c)


# def parse_bam(sample_id):
#     samfile = ps.AlignmentFile(bam_path + sample_id + '.bam', "rb")
#     prev = None
#     c = 0
#     for read in samfile.fetch():
#         if prev is None:
#             prev = read
#             continue
#         if read.reference_start < prev.reference_end:
#             c += 1
#         if read.reference_end > prev.reference_end:
#             prev = read
#     print("%s %d" % (sample_id, c))


def parse_bam(sample_id, ref_seq):
    samfile = ps.AlignmentFile(bam_path + sample_id + '.bam', "rb")
    pos_to_variants = {}
    pos_to_indels = {}
    for contig in samfile.fetch():
        ref_pos = contig.reference_start
        contig_pos = 0
        for op, l in contig.cigar:
            if op == 0:
                #match
                for i in range(l):
                    var_list = pos_to_variants.get(ref_pos + i)
                    if var_list is None:
                        pos_to_variants[ref_pos + i] = [(contig.query_alignment_sequence[contig_pos + i], contig.query_alignment_length)]
                        pos_to_indels[ref_pos + i] = [(0, contig.query_alignment_length)]
                    else:
                        var_list.append((contig.query_alignment_sequence[contig_pos + i], contig.query_alignment_length))
                        pos_to_indels[ref_pos + i].append((0, contig.query_alignment_length))
                ref_pos += l
                contig_pos += l
            elif op == 1:
                #insert
                indels = pos_to_indels.get(ref_pos)
                
        if contig.reference_start < prev.reference_end:
            c += 1
        if contig.reference_end > prev.reference_end:
            prev = contig
    print("%s %d" % (sample_id, c))


def parse_bams():
    ref_seq = read_h37rv()
    for sample_id in open(path_to_ids).readlines():
        parse_bam(sample_id.strip(), ref_seq)


if __name__ == '__main__':
    # run_variant_calling()
    # parse_variants()
    parse_bams()

