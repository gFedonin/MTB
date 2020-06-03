from bisect import bisect_left
import pysam
from os.path import exists

from Bio import SeqIO
from sklearn.externals.joblib import Parallel, delayed

from core.annotations import path_to_annotations
from core.constants import data_path, dr_genes, upstream_length
from core.data_reading import read_h37rv

path_to_snps = data_path + 'snps/combined_raw_variants_mq40_keep_complex_std_names_filtered/'
# path_to_ids = data_path + 'combined_no_dup_snp_with_pheno.list'
path_to_ids = data_path + 'casali14/casali14.list'
# out_path_aln_snp = data_path + 'combined_mq40_keep_complex_std_names_filtered_with_DR_snp.fasta'
# out_path_aln_indel = data_path + 'combined_mq40_keep_complex_std_names_filtered_with_DR_indels.fasta'
out_path_aln_snp = data_path + 'casali14/casali_snps_no_DR.fasta'
out_path_aln_indel = data_path + 'casali14/casali_indels_no_DR.fasta'
# out_path_snps = data_path + 'combined_mq40_keep_complex_std_names_filtered_with_DR.snp_list'
# out_path_indels = data_path + 'combined_mq40_keep_complex_std_names_filtered_with_DR.indel_list'
out_path_snps = data_path + 'casali14/casali.snp_list'
out_path_indels = data_path + 'casali14/casali.indel_list'
path_to_canetti = data_path + 'canetti.fasta'
path_to_canetti_bam = data_path + 'canetti.bam'

# thread_num = 144
overwrite = False

filter_out_DR_genes = True
filter_out_PGRS = True
filter_out_recombination_hotspots = True
path_to_recombination_hotspots = data_path + 'list_of_recombining_genes.txt'


def get_intervals_to_filter_out():
    """
    Reads annotations and picks coordinates of DR genes.

    :param filter_out_DR_genes: if True adds coordinates of DR genes from list of DR genes to the results
    :param filter_short_repeats: if True adds coordinates of short tandem repeats to the result
    :return: list of tuples (begin, end), representing intervals
    """
    recombining_genes = set(l.strip() for l in open(path_to_recombination_hotspots).readlines())
    coords = []
    with open(path_to_annotations, 'r') as f:
        for line in f.readlines():
            s = line.strip().split('\t')
            if s[2] == 'Gene':
                strand = s[6]
                if filter_out_DR_genes:
                    for gene_name in dr_genes:
                        if gene_name in s[-1]:
                            if strand == '+':
                                coords.append((int(s[3]) - upstream_length, int(s[4])))
                            else:
                                coords.append((int(s[3]), int(s[4]) + upstream_length))
                if filter_out_PGRS:
                    if 'PGRS' in s[-1] or 'PE' in s[-1]:
                        if strand == '+':
                            coords.append((int(s[3]) - upstream_length, int(s[4])))
                        else:
                            coords.append((int(s[3]), int(s[4]) + upstream_length))
                if filter_out_recombination_hotspots:
                    for gene_name in recombining_genes:
                        if gene_name in s[-1]:
                            if strand == '+':
                                coords.append((int(s[3]) - upstream_length, int(s[4])))
                            else:
                                coords.append((int(s[3]), int(s[4]) + upstream_length))
    coords.sort(key=lambda tup: tup[0])
    return coords


def parse_canetti(h37rv):
    canetti = str([r for r in SeqIO.parse(path_to_canetti, 'fasta')][0].seq)
    samfile = pysam.AlignmentFile(path_to_canetti_bam, "rb")
    snps = {}
    indels = set()
    for read in samfile.fetch():
        ref_pos = read.reference_start
        read_pos = read.query_alignment_start
        # print('ref_pos = %d, contig_pos = %d' % (ref_pos, read_pos))
        # print(read.cigarstring)
        # if read.is_reverse:
        #     print('reverse')
        for op, l in read.cigartuples:
            if op == 0:
                #match
                # print(ref_seq[ref_pos: ref_pos + l])
                # print(mut_str[read_pos: read_pos + l])
                for i in range(l):
                    if h37rv[ref_pos] != canetti[read_pos]:
                        snps[ref_pos + 1] = canetti[read_pos]
                    ref_pos += 1
                    read_pos += 1
            elif op == 1:
                #insert
                indels.add(str(ref_pos + 1) + '\t \t' + canetti[read_pos: read_pos + l])
                read_pos += l
            elif op == 2:
                #delete
                # raw_variants.append(Variant(ref_pos, ref_seq[ref_pos: ref_pos + l], ''))
                for i in range(l):
                    indels.add(str(ref_pos + i + 1) + '_' + h37rv[ref_pos + i] + '')
                ref_pos += l
            elif op == 5:
                #hardclip
                read_pos += l
    samfile.close()
    return snps, indels


def read_vars(sample_id, filter_intervals):
    # sample_id = fname[0: fname.rfind('.variants')]
    snps = {}
    indels = set()
    filter_intervals_starts = [x[0] for x in filter_intervals]
    with open(path_to_snps + sample_id + '.variants') as f1:
        for line in f1.readlines():
            if line[0] == '#':
                continue
            s = line.strip().split('\t')
            pos = int(s[0])
            if len(s[-1]) == 1 and len(s[-2]) == 1:
                inside_filtered_interval = False
                i = bisect_left(filter_intervals_starts, pos)
                if i == 0:
                    if pos == filter_intervals_starts[0]:
                        inside_filtered_interval = True
                else:
                    if i < len(filter_intervals) and filter_intervals_starts[i] == pos:
                        inside_filtered_interval = True
                    elif pos <= filter_intervals[i - 1][1]:
                        inside_filtered_interval = True
                if not inside_filtered_interval:
                    snps[pos] = s[2]
            else:
                # indel
                start = pos
                end = pos + len(s[1])
                inside_filtered_interval = False
                i = bisect_left(filter_intervals_starts, start)
                if i == 0:
                    if start == filter_intervals_starts[0]:
                        inside_filtered_interval = True
                else:
                    if i < len(filter_intervals) and filter_intervals_starts[i] <= end:
                        inside_filtered_interval = True
                    elif start <= filter_intervals[i - 1][1]:
                        inside_filtered_interval = True
                if not inside_filtered_interval:
                    indels.add(str(pos) + '\t' + s[1] + '\t' + s[2])
    return sample_id, snps, indels


def read_snps_no_filter(sample_id):
    # sample_id = fname[0: fname.rfind('.variants')]
    snps = {}
    indels = set()
    with open(path_to_snps + sample_id + '.variants') as f1:
        for line in f1.readlines():
            if line[0] == '#':
                continue
            s = line.strip().split('\t')
            pos = int(s[0])
            if len(s[-1]) == 1 and len(s[-2]) == 1:
                snps[pos] = s[2]
            else:
                indels.add(str(pos) + '\t' + s[1] + '\t' + s[2])
    return sample_id, snps, indels


def main():
    sample_to_snps = {}
    sample_to_indels = {}
    all_snp_pos = set()
    all_indels = set()

    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    # fnames = [fname for fname in listdir(path_to_snps) if fname.endswith('.variants')]

    h37rv = read_h37rv()
    snps_canetti, indels_canetti = parse_canetti(h37rv)

    if filter_out_DR_genes:
        filter_intervals = get_intervals_to_filter_out()
        tasks = Parallel(n_jobs=-1)(delayed(read_vars)(sample_id, filter_intervals) for sample_id in sample_ids
                                    if exists(path_to_snps + sample_id + '.variants'))
    else:
        tasks = Parallel(n_jobs=-1)(delayed(read_snps_no_filter)(sample_id) for sample_id in sample_ids
                                    if exists(path_to_snps + sample_id + '.variants'))


    for sample_id, snps, indels in tasks:
        sample_to_snps[sample_id] = snps
        sample_to_indels[sample_id] = indels
        for pos, alt in snps.items():
            all_snp_pos.add(pos)
        all_indels.update(indels)
    all_snps = list(all_snp_pos)
    all_snps.sort()
    all_indels = list(all_indels)
    all_indels.sort(key=lambda x: int(x.split()[0]))

    if overwrite or not exists(out_path_snps):
        with open(out_path_snps, 'w') as f:
            for pos in all_snps:
                f.write(str(pos))
                f.write('\n')
    if overwrite or not exists(out_path_indels):
        with open(out_path_indels, 'w') as f:
            for indel in all_indels:
                f.write(indel)
                f.write('\n')

    if overwrite or not exists(out_path_aln_snp):
        with open(out_path_aln_snp, 'w') as f_snp:
            f_snp.write('>H37Rv\n')
            for snp_pos in all_snps:
                f_snp.write(h37rv[snp_pos - 1])
            f_snp.write('\n')
            f_snp.write('>canetti\n')
            for snp_pos in all_snps:
                snp_letter = snps_canetti.get(snp_pos)
                if snp_letter is not None:
                    f_snp.write(snp_letter)
                else:
                    f_snp.write(h37rv[snp_pos - 1])
            f_snp.write('\n')
            for sample_id, snps in sample_to_snps.items():
                f_snp.write('>' + sample_id + '\n')
                for snp_pos in all_snps:
                    snp_letter = snps.get(snp_pos)
                    if snp_letter is not None:
                        f_snp.write(snp_letter)
                    else:
                        f_snp.write(h37rv[snp_pos - 1])
                f_snp.write('\n')

    if overwrite or not exists(out_path_aln_indel):
        with open(out_path_aln_indel, 'w') as f_indel:
            f_indel.write('>H37Rv\n')
            for i in range(len(all_indels)):
                f_indel.write('0')
            f_indel.write('\n')
            f_indel.write('>canetti\n')
            for indel in all_indels:
                if indel in indels_canetti:
                    f_indel.write('1')
                else:
                    f_indel.write('0')
            f_indel.write('\n')
            for sample_id, snps in sample_to_snps.items():
                indels = sample_to_indels[sample_id]
                f_indel.write('>' + sample_id + '\n')
                for indel in all_indels:
                    if indel in indels:
                        f_indel.write('1')
                    else:
                        f_indel.write('0')
                f_indel.write('\n')


if __name__ == '__main__':
    main()
