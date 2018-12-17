from bisect import bisect_left, bisect_right
import numpy as np
from sklearn.externals.joblib import Parallel, delayed

from src.core.constants import upstream_length

path_to_snps = '/export/data/kchukreev/data/mutations_files/unformated_vcfs/realigned_vcfs/'
path_to_annotations = './data/AL123456_rev.gff'
path_to_ids = './data/Full_subset_filtered_snp_pheno2.txt'
path_to_depths = './data/coverages_with_percentiles/5/'

filter_by_upper_percentile = True

qual_thresold = 40
snp_cov_threshold = 15
sample_cov_threshold = 0.9


def read_annotations():
    repeats_and_mobile_elements = []
    cds = []
    with open(path_to_annotations, 'r') as f:
        for line in f.readlines():
            s = line.split('\t')
            type = s[2]
            start = int(s[3])
            end = int(s[4])
            strand = 1 if s[6] == '+' else -1
            if type == 'Gene':
                cds.append((start, end))
                if strand == 1:
                    cds.append((start - upstream_length, start))
                else:
                    cds.append((end, end + upstream_length))
            elif type == 'Repeat_region':
                repeats_and_mobile_elements.append((start, end))
            elif type == 'mobile_element':
                repeats_and_mobile_elements.append((start, end))
    cds.sort(key=lambda tup: tup[0])
    repeats_and_mobile_elements.sort(key=lambda tup: tup[0])
    return cds, repeats_and_mobile_elements


def read_snps(sample_id, repeats, repeats_starts, percentile):
    snps = []
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
                    cov = info[i + 3: j]
                    if int(cov) > snp_cov_threshold:
                        if filter_by_upper_percentile and cov > percentile:
                            continue
                        snp_pos = int(tokens[1])
                        in_repeat = False
                        i = bisect_left(repeats_starts, snp_pos)
                        if i == 0:
                            if snp_pos == repeats[0][0]:
                                in_repeat = True
                        else:
                            rep = repeats[i - 1]
                            if rep[1] >= snp_pos:
                                in_repeat = True
                        if not in_repeat:
                            snps.append((snp_pos, tokens[4]))
    return sample_id, snps


def read_coverage(sample_id):
    coverage = []
    with open(path_to_depths + sample_id + '.coverage') as f:
        f.readline()
        percentile = f.readline().split('\t')[-1]
        f.readline()
        for line in f.readlines():
            s = line.split()
            coverage.append((int(s[0]), int(s[1])))
    coverage.sort(key=lambda tup: tup[0])
    return sample_id, coverage, percentile


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


def count_snp(snps, cds_list, cds_starts, uncovered_snp_pos):
    snps_num = 0
    for j in range(len(snps)):
        pos = snps[j]
        if pos in uncovered_snp_pos:
            continue
        i = bisect_left(cds_starts, pos)
        cds = 0
        if i == 0:
            if cds_list[0][0] == pos:
                cds = 1
        else:
            for k in range(i - 1, -1):
                if cds_list[k][1] >= pos:
                    cds += 1
                else:
                    break
        if cds > 1:
            snps_num += 1
    return snps_num


def main():
    sample_to_snps = {}
    all_snp_pos = set()
    sample_to_cov = {}
    sample_to_cov_list = {}
    percentiles = {}

    sample_ids = [sample_id[:-1] for sample_id in open(path_to_ids, 'r').readlines()]
    sample_num = len(sample_ids)

    cds, repeats_and_mobile_elements = read_annotations()
    repeats_starts = [x[0] for x in repeats_and_mobile_elements]
    cds_starts = [x[0] for x in cds]

    coverages = Parallel(n_jobs=-1)(delayed(read_coverage)(sample_id) for sample_id in sample_ids)
    for sample_id, cov, percentile in coverages:
        sample_to_cov_list[sample_id] = cov
        percentiles[sample_id] = percentile

    tasks = Parallel(n_jobs=-1)(
        delayed(read_snps)(sample_id, repeats_and_mobile_elements, repeats_starts, percentiles[sample_id]) for sample_id in sample_ids)
    for sample_id, snps in tasks:
        sample_to_snps[sample_id] = snps
        for snp_pos, alt in snps:
            all_snp_pos.add(snp_pos)
    all_snps = list(all_snp_pos)
    all_snps.sort()
    snp_num = len(all_snps)
    print('done with snp')

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

    uncovered_snp_pos = set()
    for i in range(snp_num):
        if snp_cov[i]/sample_num < sample_cov_threshold:
            uncovered_snp_pos.add(all_snps[i])

    print(count_snp(all_snps, cds, cds_starts, uncovered_snp_pos))


if __name__ == '__main__':
    main()
