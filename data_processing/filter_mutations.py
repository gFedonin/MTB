from os import mkdir
from pathlib import Path

from Bio import SeqIO
import numpy as np
from bisect import bisect_left, bisect_right

from os.path import exists
from sklearn.externals.joblib import Parallel, delayed

from src.core.constants import upstream_length, dr_genes

path_to_ids = '../../data/all_with_pheno_and_snp.txt'#'./data/test.list'
path_to_snps = '../../data/snps/realigned_vcfs_indels/'
out_path = '../../data/snps/raw_with_DR_with_indel_with_pheno_and_snp_mc10/'
out_path_snp = out_path + 'raw_snp_pos.csv'
path_to_depths = ''#''./data/coverages_with_percentiles/5/'
path_to_annotations = '../../data/AL123456_rev.gff'
path_to_low_coverage_ids = ''#'./data/low_covered21.txt'

qual_thresold = 40
snp_cov_threshold = 10
sample_cov_threshold = 0.9


filter_by_upper_percentile = False
filter_DR_genes = False
filter_by_snp_sample_coverage = False
filter_by_sample_coverage = False
snp_only = False


def get_genes_coords():
    coords = []
    with open(path_to_annotations, 'r') as f:
        for line in f.readlines():
            s = line.split('\t')
            if s[2] == 'Gene':
                strand = s[6]
                for gene_name in dr_genes:
                    if gene_name in s[-1]:
                        if strand == '+':
                            coords.append((int(s[3]) - upstream_length, int(s[4])))
                        else:
                            coords.append((int(s[3]), int(s[4]) + upstream_length))
            elif s[2] in ('Repeat_region', 'mobile_element'):
                coords.append((int(s[3]), int(s[4])))
    return coords


def get_repeats_coords():
    coords = []
    with open(path_to_annotations, 'r') as f:
        for line in f.readlines():
            s = line.split('\t')
            if s[2] in ('Repeat_region', 'mobile_element'):
                coords.append((int(s[3]), int(s[4])))
    return coords


def read_coverage(sample_id):
    coverage = []
    with open(path_to_depths + sample_id + '.coverage') as f:
        f.readline()
        percentile = float(f.readline().split('\t')[-1])
        f.readline()
        for line in f.readlines():
            s = line.split()
            coverage.append((int(s[0]), int(s[1])))
    coverage.sort(key=lambda tup: tup[0])
    return sample_id, coverage, percentile


def pos_covered(coverage, pos):
    for start, end in coverage:
        if start <= pos <= end:
            return True
    return False


def parse_cigar(cigar):
    start = 0
    res = []
    for i in range(len(cigar)):
        if not cigar[i].isdigit():
            res.append((int(cigar[start:i]), cigar[i]))
            start = i + 1
    return res


def read_snps(sample_id, coords, gene_starts, percentile):
    variants = []
    with open(path_to_snps + sample_id + '_h37rv.vcf', 'r') as f1:
        lines = f1.readlines()
        for line in lines:
            line = line.strip()
            if line[0] != '#':
                tokens = line.split('\t')
                qual = float(tokens[5])
                if qual >= qual_thresold:
                    info = tokens[-3]
                    i = info.index('DP=')
                    j = info.index(';', i + 3)
                    cov = int(info[i + 3: j])
                    if cov >= snp_cov_threshold:
                        if filter_by_upper_percentile and cov > percentile:
                            continue
                        i = info.index('TYPE=')
                        j = info.index(';', i + 5)
                        v_type = info[i + 5: j]
                        pos = int(tokens[1])
                        inside_genes = False
                        if filter_DR_genes:
                            i = bisect_left(gene_starts, pos)
                            if i == 0:
                                if pos == gene_starts[0]:
                                    inside_genes = True
                            else:
                                if pos <= coords[i - 1][1]:
                                    inside_genes = True
                        if not inside_genes:
                            alt = tokens[4]
                            if v_type == 'snp':
                                if len(alt) == 1:
                                    variants.append((pos, alt, v_type))
                                else:
                                    i = info.index('CIGAR=')
                                    j = info.index(';', i + 6)
                                    cigar = parse_cigar(info[i + 6: j])
                                    p = 0
                                    for l, o in cigar:
                                        if o == 'M':
                                            pos += l
                                            p += l
                                        elif o == 'X':
                                            for n in range(l):
                                                variants.append((pos, alt[p], 'snp'))
                                                pos += 1
                                                p += 1
                            else:
                                if snp_only:
                                    continue
                                i = info.index('CIGAR=')
                                j = info.index(';', i + 6)
                                cigar = parse_cigar(info[i + 6: j])
                                p = 0
                                for l, o in cigar:
                                    if o == 'M':
                                        pos += l
                                        p += l
                                    elif o == 'X':
                                        for n in range(l):
                                            variants.append((pos, alt[p], 'snp'))
                                            pos += 1
                                            p += 1
                                    elif o == 'D':
                                        for n in range(l):
                                            variants.append((pos, '-', 'del'))
                                            pos += 1
                                    elif o == 'I':
                                        variants.append((pos, alt[p: p + l], 'ins'))
                                        p += l
    return sample_id, variants


def sample_cov(sample_id, all_snps, coverage):
    # coverage = read_coverage(sample_id)
    cov = 0
    lo = 0
    # print(sample_id + ' ' + str(len(coverage)))
    for start, end in coverage:
        s = bisect_left(all_snps, start, lo=lo)
        e = bisect_right(all_snps, end, lo=s + 1)
        # print(str(start) + ' ' + str(end) + ' ' + str(s) + ' ' + str(e))
        cov += e - s
        lo = e
    # print(sample_id + ' cov ' + str(cov))
    return sample_id, cov / len(all_snps)


def compute_snp_cov(all_snps, coverage):
    snp_cov = np.zeros(len(all_snps), dtype=int)
    lo = 0
    for start, end in coverage:
        s = bisect_left(all_snps, start, lo=lo)
        if s < len(all_snps):
            e = bisect_right(all_snps, end, lo=s + 1)
            # print(str(start) + ' ' + str(end) + ' ' + str(s) + ' ' + str(e))
            snp_cov[s:e] = 1
            lo = e
        else:
            break
    return snp_cov


def main():
    if not exists(out_path):
        mkdir(out_path)
    with open(out_path + 'config.txt', 'w') as f:
        f.write('qual_thresold = ' + str(qual_thresold) + '\n')
        f.write('snp_cov_threshold = ' + str(snp_cov_threshold) + '\n')
        f.write('sample_cov_threshold = ' + str(sample_cov_threshold) + '\n')
        f.write('gene_upstream_threshold = ' + str(upstream_length) + '\n')
        f.write('filter_by_upper_percentile = ' + str(filter_by_upper_percentile) + '\n')
        f.write('filter_DR_genes = ' + str(filter_DR_genes) + '\n')
        f.write('filter_by_snp_sample_coverage = ' + str(filter_by_snp_sample_coverage) + '\n')
        f.write('filter_by_sample_coverage = ' + str(filter_by_sample_coverage) + '\n')
        f.write('snp only = ' + str(snp_only) + '\n')

    sample_to_snps = {}
    all_snp_pos = set()
    sample_to_cov = {}

    sample_ids = [sample_id[:-1] for sample_id in open(path_to_ids, 'r').readlines()]
    sample_num = len(sample_ids)

    if filter_DR_genes:
        filter_coords = get_genes_coords()
    else:
        filter_coords = get_repeats_coords()
    filter_coords.sort(key=lambda tup: tup[0])
    gene_starts = [x[0] for x in filter_coords]

    sample_to_cov_list = {}
    percentiles = {}

    if filter_by_upper_percentile or filter_by_sample_coverage or filter_by_snp_sample_coverage:
        coverages = Parallel(n_jobs=-1)(delayed(read_coverage)(sample_id) for sample_id in sample_ids)
        for sample_id, cov, percentile in coverages:
            sample_to_cov_list[sample_id] = cov
            percentiles[sample_id] = percentile
        print('read sample coverages')

    if filter_by_upper_percentile:
        tasks = Parallel(n_jobs=-1)(delayed(read_snps)(sample_id, filter_coords, gene_starts, percentiles[sample_id]) for sample_id in sample_ids)
    else:
        tasks = Parallel(n_jobs=-1)(
            delayed(read_snps)(sample_id, filter_coords, gene_starts, 0) for sample_id in
            sample_ids)
    for sample_id, snps in tasks:
        sample_to_snps[sample_id] = snps
        for pos, alt, v_type in snps:
            all_snp_pos.add(pos)
    all_snps = list(all_snp_pos)
    all_snps.sort()
    snp_num = len(all_snps)
    print('snp reading done')

    if filter_by_sample_coverage:
        tasks = Parallel(n_jobs=-1)(delayed(sample_cov)(sample_id, all_snps, sample_to_cov_list[sample_id]) for sample_id in sample_ids)
        for sample_id, cov in tasks:
            sample_to_cov[sample_id] = cov
        print('done with sample coverages')

    uncovered_snp_pos = set()
    if filter_by_snp_sample_coverage:
        snp_cov_list = Parallel(n_jobs=-1)(
            delayed(compute_snp_cov)(all_snps, sample_to_cov_list[sample_id]) for sample_id in sample_ids)
        snp_cov = np.zeros(snp_num, dtype=int)
        for snp_cov_arr in snp_cov_list:
            snp_cov += snp_cov_arr
        print('done with snp coverages by samples')

        with open('./snp_cov.txt', 'w')as f:
            for i in range(snp_num):
                f.write(str(all_snps[i]) + ' ' + str(snp_cov[i]) + '\n')

        for i in range(snp_num):
            if snp_cov[i]/sample_num < sample_cov_threshold:
                uncovered_snp_pos.add(all_snps[i])

    with open(out_path_snp, 'w') as f:
        for pos in all_snps:
            if pos not in uncovered_snp_pos:
                f.write(str(pos) + '\n')
    print('printed all snps')

    if filter_by_sample_coverage:
        with open(path_to_low_coverage_ids, 'w') as f:
            for sample_id in sample_ids:
                if sample_to_cov[sample_id] < sample_cov_threshold:
                    f.write(sample_id + '\n')

    for sample_id, snps in sample_to_snps.items():
        if not filter_by_sample_coverage or sample_to_cov[sample_id] >= sample_cov_threshold:
            with open(out_path + sample_id + '.variants', 'w') as f:
                for pos, alt, v_type in snps:
                    f.write(str(pos) + '\t' + alt + '\t' + v_type + '\n')


if __name__ == '__main__':
    main()
