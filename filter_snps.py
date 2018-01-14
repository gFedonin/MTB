from Bio import SeqIO
import numpy as np
from bisect import bisect_left, bisect_right

from sklearn.externals.joblib import Parallel, delayed

path_to_ids = './data/Full_subset_filtered_snp_pheno2.txt'
path_to_snps = '/export/data/kchukreev/data/mutations_files/unformated_vcfs/realigned_vcfs/'
out_path = './data/snps/raw_no_DR/'
out_path_snp = './data/snps/raw_snp_pos_no_DR_filtered2.csv'
path_to_depths = './data/coverages_with_percentiles/5/'
path_to_annotations = './data/AL123456_rev.gff'
path_to_low_coverage_ids = './data/low_covered21.txt'

qual_thresold = 40
snp_cov_threshold = 15
sample_cov_threshold = 0.9

gene_upstream_threshold = 100

filter_by_upper_percentile = True
filter_DR_genes = False


list_genes = ['ahpC', 'eis', 'embA', 'embB', 'embC', 'embR', 'fabG1', 'gid', 'gyrA', 'gyrB', 'inhA', 'iniA', 'iniC',
             'katG', 'manB', 'ndh', 'pncA', 'rmlD', 'rpoB', 'rpsA', 'rpsL', 'rrs', 'tlyA', 'ethR', 'fpbC', 'iniB',
              'kasA', 'ethA', 'fabD', 'efpA', 'thyA', 'panD', 'accD6', 'fbpC', 'nat', 'folC', 'rrl', 'rpoC', 'ribD', 'rplC']


def get_genes_coords():
    coords = []
    with open(path_to_annotations, 'r') as f:
        for line in f.readlines():
            s = line.split('\t')
            if s[2] == 'Gene':
                strand = s[6]
                for gene_name in list_genes:
                    if gene_name in s[-1]:
                        if strand == '+':
                            coords.append((int(s[3]) - gene_upstream_threshold, int(s[4])))
                        else:
                            coords.append((int(s[3]), int(s[4]) + gene_upstream_threshold))
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


def read_snps(sample_id, coords, gene_starts, percentile):
    snps = {}
    with open(path_to_snps + sample_id + '_h37rv.vcf', 'r') as f1:
        lines = f1.readlines()
        for line in lines:
            line = line[:-1]
            if line[0] != '#':
                tokens = line.split('\t')
                qual = float(tokens[5])
                if qual >= qual_thresold:
                    info = tokens[-3]
                    i = info.index('DP=')
                    j = info.index(';', i + 3)
                    cov = int(info[i + 3: j])
                    if cov > snp_cov_threshold:
                        if filter_by_upper_percentile and cov > percentile:
                            continue
                        snp_pos = int(tokens[1])
                        inside_genes = False
                        if filter_DR_genes:
                            i = bisect_left(gene_starts, snp_pos)
                            if i == 0:
                                if snp_pos == gene_starts[0]:
                                    inside_genes = True
                            else:
                                if snp_pos <= coords[i - 1][1]:
                                    inside_genes = True
                            # for coord in coords:
                            #     if coord[0] <= snp_pos <= coord[1]:
                            #         inside_genes = True
                            #         break
                        if not inside_genes:
                            snps[snp_pos] = tokens[4]
    return sample_id, snps


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

    coverages = Parallel(n_jobs=-1)(delayed(read_coverage)(sample_id) for sample_id in sample_ids)
    sample_to_cov_list = {}
    percentiles = {}
    for sample_id, cov, percentile in coverages:
        sample_to_cov_list[sample_id] = cov
        percentiles[sample_id] = percentile
    print('read sample coverages')

    tasks = Parallel(n_jobs=-1)(delayed(read_snps)(sample_id, filter_coords, gene_starts, percentiles[sample_id]) for sample_id in sample_ids)
    for sample_id, snps in tasks:
        sample_to_snps[sample_id] = snps
        for snp_pos in snps.keys():
            all_snp_pos.add(snp_pos)
    all_snps = list(all_snp_pos)
    all_snps.sort()
    snp_num = len(all_snps)
    print('snp reading done')

    tasks = Parallel(n_jobs=-1)(delayed(sample_cov)(sample_id, all_snps, sample_to_cov_list[sample_id]) for sample_id in sample_ids)
    for sample_id, cov in tasks:
        sample_to_cov[sample_id] = cov
    # for sample_id in sample_ids:
    #     sample_id, cov = sample_cov(sample_id, all_snps)
    #     sample_to_cov[sample_id] = cov
    print('done with sample coverages')

    snp_cov_list = Parallel(n_jobs=-1)(
        delayed(compute_snp_cov)(all_snps, sample_to_cov_list[sample_id]) for sample_id in sample_ids)
    snp_cov = np.zeros(snp_num, dtype=int)
    for snp_cov_arr in snp_cov_list:
        snp_cov += snp_cov_arr
    print('done with snp coverages by samples')

    with open('./snp_cov.txt', 'w')as f:
        for i in range(snp_num):
            f.write(str(all_snps[i]) + ' ' + str(snp_cov[i]) + '\n')

    uncovered_snp_pos = set()
    for i in range(snp_num):
        if snp_cov[i]/sample_num < sample_cov_threshold:
            uncovered_snp_pos.add(all_snps[i])

    with open(out_path_snp, 'w') as f:
        for snp_pos in all_snps:
            if snp_pos not in uncovered_snp_pos:
                f.write(str(snp_pos) + '\n')
    print('printed all snps')

    with open(path_to_low_coverage_ids, 'w') as f:
        for sample_id in sample_ids:
            if sample_to_cov[sample_id] < sample_cov_threshold:
                f.write(sample_id + '\n')

    for sample_id, snps in sample_to_snps.items():
        if sample_to_cov[sample_id] >= sample_cov_threshold:
            with open(out_path + sample_id + '.snp', 'w') as f:
                for pos, alt in snps.items():
                    f.write(str(pos) + '\t' + alt + '\n')


if __name__ == '__main__':
    main()
