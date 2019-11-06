import gzip
from bisect import bisect_left
from os import listdir, makedirs
from os.path import exists
from sklearn.externals.joblib import Parallel, delayed

from core.annotations import path_to_annotations
from core.constants import data_path, dr_genes, upstream_length

# path_to_vcf = data_path + 'snps/gatk_before_cortex/vcf_split_gatk/'
# out_path = data_path + 'snps/gatk_before_cortex/raw_variants_fixed_no_rep_gatk/'
# path_to_vcf = data_path + 'snps/gatk_before_cortex/vcf/'
# out_path = data_path + 'snps/gatk_before_cortex/raw_variants_no_gvcf_mc10/'
path_to_vcf = '/export/data/kchukreev/data/mutations_files/cortex_vcfs/results/'
out_path = data_path + 'snps/cortex/raw_variants/'


min_alt_frac = 0.75
min_dp = 10

thread_num = 32

filter_annotation_repeats = False
filter_short_repeats = False
path_to_short_tandem_repeats = data_path + 'h37rv.fasta.2.7.7.80.10.20.10.dat'
filter_out_DR_genes = False


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


def parse(fname, filter_intervals):
    filter_intervals_starts = [x[0] for x in filter_intervals]
    i = fname.find('_')
    if i == -1:
        i = fname.find('.')
    sample_name = fname[0:i]
    with open(out_path + sample_name + '.variants', 'w') as fout:
        if fname.endswith('gz'):
            f = gzip.open(path_to_vcf + fname, 'rt')
            lines = f
        else:
            f = open(path_to_vcf + fname)
            lines = f.readlines()
        fout.write('# source ' + path_to_vcf + fname + '\n')
        res = set()
        for l in lines:
            if l.startswith('#'):
                continue
            s = l.strip().split('\t')
            pos = int(s[1])
            ref = s[3]
            alt = s[4]
            stats = s[-1].split(':')
            read_counts = [int(c) for c in stats[1].split(',')]
            total_depth = sum(read_counts)
            if total_depth < min_dp:
                continue
            if read_counts[-1]/total_depth >= min_alt_frac and read_counts[-1] >= min_dp:
                if not in_filter_interval(pos, filter_intervals_starts, filter_intervals):
                    if len(ref) != 1:
                        for i in range(1, len(ref)):
                            res.add((pos + i, '-', 'del'))
                    elif len(alt) != 1:
                        res.add((pos, alt[1:], 'ins'))
                    else:
                        res.add((pos, alt, 'snp'))
        f.close()
        var_list = list(res)
        var_list.sort(key=lambda x: x[0])
        for pos, alt, t in var_list:
            fout.write('%d\t%s\t%s\n' % (pos, alt, t))
    return 1


if __name__ == "__main__":
    if not exists(out_path):
        makedirs(out_path)
    filter_intervals = get_intervals_to_filter_out()
    names = [fname for fname in listdir(path_to_vcf) if fname.endswith('vcf.gz') or fname.endswith('.vcf')]
    tasks = Parallel(n_jobs=thread_num, batch_size=len(names)//thread_num + 1)(delayed(parse)(fname, filter_intervals)
                                                                               for fname in names)
    c = 0
    for task in tasks:
        c += task
    print('%d samples processed' % c)
    # vars = parse_complex_event(0, 'TTATGTCACCAACTTCGCC', 'CTATGTCACCAACTTCGCC')
    # for pos, alt, type in vars:
    #     print('%d %s %s' % (pos, alt, type))

