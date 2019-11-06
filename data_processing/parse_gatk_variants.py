import gzip
from bisect import bisect_left
from os import listdir, makedirs
from os.path import exists

import numpy as np
from Bio import pairwise2
from sklearn.externals.joblib import Parallel, delayed

from core.annotations import path_to_annotations
from core.constants import data_path, dr_genes, upstream_length

# path_to_vcf = data_path + 'snps/gatk_before_cortex/vcf_split_gatk/'
# out_path = data_path + 'snps/gatk_before_cortex/raw_variants_fixed_no_rep_gatk/'
# path_to_vcf = data_path + 'snps/gatk_before_cortex/vcf/'
# out_path = data_path + 'snps/gatk_before_cortex/raw_variants_no_gvcf_mc10/'
path_to_variants = data_path + 'snps/gatk_bwa_mem_rev/'
path_to_vcf = path_to_variants + 'vcf_split_gatk_fa_no_opt/'
out_path = path_to_variants + 'raw_variants_no_mq/'

min_qd = 2.0
max_fs = 60.0
# min_mq = 40.0
min_mq = 0.0
min_MQRankSum = -12.5
min_ReadPosRankSum = -8.0
max_sor = 3.0

std_mult = 2

min_alt_frac = 0.75
min_dp = 10

thread_num = 32

no_filters = False
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


def check_var(annotations):
    ann_dict = {}
    for a in annotations:
        i = a.find('=')
        ann_dict[a[0:i]] = a[i + 1:]
    qd = ann_dict.get('QD')
    if qd is not None and float(qd) < min_qd:
        return False
    fs = ann_dict.get('FS')
    if fs is not None and float(fs) > max_fs:
        return False
    mq = ann_dict.get('MQ')
    if mq is not None and float(mq) < min_mq:
        return False
    MQRankSum = ann_dict.get('MQRankSum')
    if MQRankSum is not None and float(MQRankSum) < min_MQRankSum:
        return False
    ReadPosRankSum = ann_dict.get('ReadPosRankSum')
    if ReadPosRankSum is not None and float(ReadPosRankSum) < min_ReadPosRankSum:
        return False
    sor = ann_dict.get('SOR')
    if sor is not None and float(sor) > max_sor:
        return False
    return True


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


# def parse(fname):
#     sample_name = fname[0:-7]
#     if not exists(out_path + sample_name + '.variants'):
#         with open(out_path + sample_name + '.variants', 'w') as fout:
#             with gzip.open(path_to_vcf + fname, 'rt') as f:
#                 fout.write('# source ' + path_to_vcf + fname + '\n')
#                 for l in f:
#                     if l.startswith('#'):
#                         continue
#                     s = l.strip().split('\t')
#                     annotations = s[7].split(';')
#                     if check_var(annotations):
#                         fout.write(s[1] + '\t')
#                         n = len(s[3])
#                         variants = s[4].split(',')
#                         if len(variants) == 1:
#                             fout.write(s[4] + '\t')
#                         else:
#                             scores = s[-1].split(':')[-1].split(',')
#                             for i in range(1, len(scores)):
#                                 if scores[i] == '0':
#                                     fout.write(variants[i - 1] + '\t')
#                                     break
#                         fout.write(str(n) + '\n')
#         return 1
#     else:
#         return 0


def parse(fname, filter_intervals):
    filter_intervals_starts = [x[0] for x in filter_intervals]
    i = fname.find('_')
    if i == -1:
        i = fname.find('.')
    sample_name = fname[0:i]
    with open(out_path + sample_name + '.variants', 'w') as fout:
        var_list = []
        ann_dict = {}
        if fname.endswith('gz'):
            f = gzip.open(path_to_vcf + fname, 'rt')
            lines = f
        else:
            f = open(path_to_vcf + fname)
            lines = f.readlines()
        fout.write('# source ' + path_to_vcf + fname + '\n')
        for l in lines:
            if l.startswith('#'):
                continue
            s = l.strip().split('\t')
            annotations = s[7].split(';')
            variants = s[4].split(',')
            stats = s[-1].split(':')
            # stat_names = s[-2].split(':')
            read_counts = [int(c) for c in stats[1].split(',')]
            total_depth = sum(read_counts)
            if total_depth < min_dp:
                continue
            if len(variants) == 1:
                if read_counts[-1]/total_depth >= min_alt_frac and read_counts[-1] >= min_dp:
                    var_list.append((int(s[1]), s[3], s[4], annotations))
            else:
                scores = stats[-1].split(',')
                for i in range(1, len(scores)):
                    if scores[i] == '0':
                        if read_counts[i]/total_depth >= min_alt_frac and read_counts[i] >= min_dp:
                            var_list.append((int(s[1]), s[3], variants[i - 1], annotations))
                        break
            for a in annotations:
                if ',' in a:
                    continue
                i = a.find('=')
                vals = ann_dict.get(a[0:i])
                if vals is None:
                    vals = []
                    ann_dict[a[0:i]] = vals
                vals.append(float(a[i + 1:]))
        f.close()
        res = set()
        for pos, ref, alt, ann in var_list:
            if (no_filters or check_var(ann)) and not in_filter_interval(pos, filter_intervals_starts, filter_intervals):
                if len(ref) != 1 or len(alt) != 1 or ref == '*' or alt == '*':
                    res.update(parse_complex_event(pos, ref, alt))
                else:
                    res.add((pos, alt, 'snp'))
        var_list = list(res)
        var_list.sort(key=lambda x: x[0])
        for pos, alt, t in var_list:
            fout.write('%d\t%s\t%s\n' % (pos, alt, t))
    return 1


def check_var_std(annotations, ann_means, ann_stds):
    ann_dict = {}
    for a in annotations:
        i = a.find('=')
        ann_dict[a[0:i]] = a[i + 1:]
    qd = ann_dict.get('QD')
    if qd is not None and float(qd) < ann_means['QD'] - std_mult*ann_stds['QD']:
        return False
    fs = ann_dict.get('FS')
    if fs is not None and float(fs) > ann_means['FS'] + std_mult*ann_stds['FS']:
        return False
    mq = ann_dict.get('MQ')
    if mq is not None and float(mq) < ann_means['MQ'] - std_mult*ann_stds['MQ']:
        return False
    MQRankSum = ann_dict.get('MQRankSum')
    if MQRankSum is not None and float(MQRankSum) < ann_means['MQRankSum'] - std_mult*ann_stds['MQRankSum']:
        return False
    ReadPosRankSum = ann_dict.get('ReadPosRankSum')
    if ReadPosRankSum is not None and float(ReadPosRankSum) < ann_means['ReadPosRankSum'] - std_mult*ann_stds['ReadPosRankSum']:
        return False
    sor = ann_dict.get('SOR')
    if sor is not None and float(sor) > ann_means['SOR'] + std_mult*ann_stds['SOR']:
        return False
    return True


def parse_complex_event(pos, old, new):
    res = []
    if new == '*':
        for i in range(len(old)):
            res.append((pos + i, '-', 'del'))
        return res
    if old == '*':
        res.append((pos, new, 'ins'))
        return res        
    aln = pairwise2.align.globalms(old, new, 2, -1, -1, -.1)[0]
    aln_old = aln[0]
    aln_new = aln[1]
    insert_pos = 0
    insert_s = 0
    insert = False
    p = 0
    for i in range(len(aln_old)):
        c1 = aln_old[i]
        c2 = aln_new[i]
        if c1 != c2:
            if c1 == '-':
                # ins
                if not insert:
                    insert_pos = p
                    insert_s = i
                    insert = True
            elif c2 == '-':
                #del
                if insert:
                    res.append((pos + insert_pos, aln_new[insert_s: i], 'ins'))
                    insert = False
                res.append((pos + p, c2, 'del'))
                p += 1
            else:
                if insert:
                    res.append((pos + insert_pos, aln_new[insert_s: i], 'ins'))
                    insert = False
                res.append((pos + p, c2, 'snp'))
                p += 1
        else:
            if insert:
                res.append((pos + insert_pos, aln_new[insert_s: i], 'ins'))
                insert = False
            p += 1
    if insert:
        res.append((pos + insert_pos, aln_new[insert_s:], 'ins'))
    return res


def parse_std(fname):
    sample_name = fname[0:-7]
    with open(out_path + sample_name + '.variants', 'w') as fout:
        var_list = []
        ann_dict= {}
        with gzip.open(path_to_vcf + fname, 'rt') as f:
            fout.write('# source ' + path_to_vcf + fname + '\n')
            for l in f:
                if l.startswith('#'):
                    continue
                s = l.strip().split('\t')
                annotations = s[7].split(';')
                variants = s[4].split(',')
                if len(variants) == 1:
                    var_list.append((int(s[1]), s[3], s[4], annotations))
                else:
                    scores = s[-1].split(':')[-1].split(',')
                    for i in range(1, len(scores)):
                        if scores[i] == '0':
                            var_list.append((int(s[1]), s[3], variants[i - 1], annotations))
                            break
                for a in annotations:
                    if ',' in a:
                        continue
                    i = a.find('=')
                    vals = ann_dict.get(a[0:i])
                    if vals is None:
                        vals = []
                        ann_dict[a[0:i]] = vals
                    vals.append(float(a[i + 1:]))
        ann_means = {n: np.mean(v) for n,v in ann_dict.items()}
        ann_stds = {n: np.std(v) for n, v in ann_dict.items()}
        res = set()
        for pos, ref, alt, ann in var_list:
            if check_var_std(ann, ann_means, ann_stds):
                if len(ref) != 1 or len(alt) != 1 or ref == '*' or alt == '*':
                    res.update(parse_complex_event(pos, ref, alt))
                else:
                    res.add((pos, alt, 'snp'))
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

