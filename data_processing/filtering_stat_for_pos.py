from multiprocessing.pool import Pool

import numpy as np
from sklearn.externals.joblib import Parallel, delayed

from src.core.annotations import CDSType, read_annotations, localize_all_variants
from src.core.constants import ref_len, data_path, upstream_length

path_to_filtered_vars = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_no_win_qual_mqm_std3_mqm30_filter_samples_first/'

# path_to_ids = data_path + 'all_with_pheno_and_snp.txt'#'test.list'
path_to_ids = path_to_filtered_vars + 'samples_filtered.list'
path_to_snps = data_path + 'snps/realigned_vcfs_indels_all/'
path_to_depths = data_path + 'coverages_shrinked/'#''/export/data/kkuleshov/myc/sra/'



path_to_pos_list = path_to_filtered_vars + 'killed_Walker.list'

out_path = path_to_filtered_vars + 'killed_Walker_counts_survived_samples.csv'

qual_threshold = 40
quality_std_multiplier = 3
mqm_threshold = 30
mqm_std_multiplier = 3
variant_cov_threshold = 10
sample_cov_threshold = 0.9

coverage_std_multiplier = 3  # set to x for coverage_std_threshold = x * std(coverage)
filter_by_window_coverage = False
window_len_gene = 400
window_len_nc = 50
window_min_cov_threshold = 10
window_max_cov_threshold_coef = 3  # set to x for window_max_cov_threshold = x * std(coverage)
window_std_cov_threshold = 3

thread_num = -1


def read_cov_file(sample_id):
    """
    Reads shrinked coverage .cov file for given sample id and computes coverage threshold for given percentile

    :param sample_id: sample id to read coverage
    :return: numpy array of length ref_len with coverages of positions, zero-based, upper percentile coverage threshold
    """
    coverage = np.zeros(ref_len, dtype=int)
    i = 0
    with open(path_to_depths + sample_id + '.cov') as f:
        for line in f.readlines():
            s = line.strip()
            coverage[i] = int(s)
            i += 1
    return coverage


def read_vcf_file(sample_id):
    """
    Reads variants from vcf file.
    :param sample_id: sample id
    :return: sample_id, list of tuples (variant_pos, alt, type, coverage, quality_score) representing variants
    """

    def parse_cigar(cigar):
        start = 0
        res = []
        for i in range(len(cigar)):
            if not cigar[i].isdigit():
                res.append((int(cigar[start:i]), cigar[i]))
                start = i + 1
        return res

    variants = []
    with open(path_to_snps + sample_id + '_h37rv.vcf', 'r') as f1:
        for line in f1.readlines():
            line = line.strip()
            if line[0] != '#':
                tokens = line.split('\t')
                qual = float(tokens[5])
                info = tokens[-3]
                i = info.index('DP=')
                j = info.index(';', i + 3)
                cov = int(info[i + 3: j])
                i = info.index('TYPE=')
                j = info.index(';', i + 5)
                v_type = info[i + 5: j]
                i = info.index('MQM=')
                j = info.index(';', i + 4)
                mqm = float(info[i + 4: j])
                pos = int(tokens[1])
                alt = tokens[4]
                if v_type == 'snp':
                    if len(alt) == 1:
                        variants.append((pos, alt, v_type, cov, qual, mqm))
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
                                    variants.append((pos, alt[p], 'snp', cov, qual, mqm))
                                    pos += 1
                                    p += 1
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
                                variants.append((pos, alt[p], 'snp', cov, qual, mqm))
                                pos += 1
                                p += 1
                        elif o == 'D':
                            for n in range(l):
                                variants.append((pos, '-', 'del', cov, qual, mqm))
                                pos += 1
                        elif o == 'I':
                            variants.append((pos, alt[p: p + l], 'ins', cov, qual, mqm))
                            p += l
    return sample_id, variants


def window_is_covered(pos, coverage, mean_coverage, std_coverage, cds):

    if cds is None or cds.type != CDSType.Gene:
        start = max(0, pos - window_len_nc // 2)
        end = min(pos + window_len_nc // 2, ref_len)
    else:
        start = max(0, pos - window_len_gene // 2)
        end = min(pos + window_len_gene // 2, ref_len)
    window = coverage[start: end]
    mean_cov = np.mean(window)
    std_cov = np.std(window)
    undercovered = mean_cov < window_min_cov_threshold
    overcovered = mean_cov > mean_coverage + window_max_cov_threshold_coef*std_coverage
    bigstdcov = std_cov > window_std_cov_threshold*std_coverage
    return undercovered, overcovered, bigstdcov


def quality(sample_id, target_pos_set, variants, pos_to_cds):

    pos_is_undercovered = set()
    pos_is_overcovered = set()
    win_is_undercovered = set()
    win_is_overcovered = set()
    win_has_big_std_cov = set()
    low_quality = set()
    low_quality_std = set()
    low_mqm = set()
    low_mqm_std = set()
    all_pos_in_sample = set()

    coverage = read_cov_file(sample_id)
    std_coverage = np.std(coverage)
    mean_coverage = np.mean(coverage)
    qualities = np.zeros(len(variants))
    mqm_list = np.zeros(len(variants))
    i = 0
    for pos, alt, type, cov, qual, mqm in variants:
        qualities[i] = qual
        mqm_list[i] = mqm
        i += 1
    mean_quality = np.mean(qualities)
    std_quality = np.std(qualities)
    mean_mqm = np.mean(mqm_list)
    std_mqm = np.std(mqm_list)
    for pos, alt, type, cov, qual, mqm in variants:
        if pos not in target_pos_set:
            continue
        all_pos_in_sample.add(pos)
        if qual < qual_threshold:
            low_quality.add(pos)
        if cov < variant_cov_threshold:
            pos_is_undercovered.add(pos)
        if cov > mean_coverage + coverage_std_multiplier*std_coverage:
            pos_is_overcovered.add(pos)
        if mqm < mqm_threshold:
            low_mqm.add(pos)
        if qual < mean_quality - quality_std_multiplier * std_quality:
            low_quality_std.add(pos)
        if mqm < mean_mqm - mqm_std_multiplier * std_mqm:
            low_mqm_std.add(pos)
        if filter_by_window_coverage:
            cds = pos_to_cds.get(pos)
            win_undercovered, win_overcovered, win_bigstdcov = window_is_covered(pos, coverage, mean_coverage,
                                                                    std_coverage, cds)
            if win_undercovered:
                win_is_undercovered.add(pos)
            if win_overcovered:
                win_is_overcovered.add(pos)
            if win_bigstdcov:
                win_has_big_std_cov.add(pos)

    return sample_id, low_quality, low_quality_std, low_mqm, low_mqm_std, pos_is_undercovered, pos_is_overcovered, \
           win_is_undercovered, win_is_overcovered, win_has_big_std_cov, all_pos_in_sample


if __name__ == '__main__':
    sample_ids = list(line.strip() for line in open(path_to_ids).readlines())
    target_pos_list = list(int(line.strip()) for line in open(path_to_pos_list).readlines())
    target_pos_set = set(target_pos_list)
    sample_to_variants = Parallel(n_jobs=thread_num)(delayed(read_vcf_file)(sample_id) for sample_id in sample_ids)
    cds_list = read_annotations(upstream_length)
    pos_to_cds = localize_all_variants(target_pos_list, cds_list)

    print('starting quality computations')
    # with Pool(thread_num) as pool:
    #     tasks = [pool.apply_async(quality, (sample_id, target_pos_set, variants, pos_to_cds)) for sample_id, variants in sample_to_variants]
    tasks = Parallel(n_jobs=thread_num)(delayed(quality)(sample_id, target_pos_set, variants, pos_to_cds)
                                        for sample_id, variants in sample_to_variants)
    #, batch_size=len(sample_ids)//thread_num + 1
    low_quality_counts = {pos: 0 for pos in target_pos_list}
    low_quality_std_counts = {pos: 0 for pos in target_pos_list}
    low_mqm_counts = {pos: 0 for pos in target_pos_list}
    low_mqm_std_counts = {pos: 0 for pos in target_pos_list}
    pos_is_undercovered_counts = {pos: 0 for pos in target_pos_list}
    pos_is_overcovered_counts = {pos: 0 for pos in target_pos_list}
    pos_counts = {pos: 0 for pos in target_pos_list}
    if filter_by_window_coverage:
        win_is_undercovered_counts = {pos: 0 for pos in target_pos_list}
        win_is_overcovered_counts = {pos: 0 for pos in target_pos_list}
        win_has_big_std_cov_counts = {pos: 0 for pos in target_pos_list}

    for sample_id, low_quality, low_quality_std, low_mqm, low_mqm_std, pos_is_undercovered, pos_is_overcovered, \
           win_is_undercovered, win_is_overcovered, win_has_big_std_cov, all_pos_in_sample in tasks:
    # for task in tasks:
    #     sample_id, low_quality, low_quality_std, low_mqm, low_mqm_std, pos_is_undercovered, pos_is_overcovered, \
    #     win_is_undercovered, win_is_overcovered, win_has_big_std_cov, all_pos_in_sample = task.get(timeout=10000)
        for pos in all_pos_in_sample:
            pos_counts[pos] += 1
        for pos in low_quality:
            low_quality_counts[pos] += 1
        for pos in low_quality_std:
            low_quality_std_counts[pos] += 1
        for pos in low_mqm:
            low_mqm_counts[pos] += 1
        for pos in low_mqm_std:
            low_mqm_std_counts[pos] += 1
        for pos in pos_is_overcovered:
            pos_is_overcovered_counts[pos] += 1
        for pos in pos_is_undercovered:
            pos_is_undercovered_counts[pos] += 1
        if filter_by_window_coverage:
            for pos in win_is_overcovered:
                win_is_overcovered[pos] += 1
            for pos in win_is_undercovered:
                win_is_undercovered[pos] += 1
            for pos in win_has_big_std_cov:
                win_has_big_std_cov_counts[pos] += 1

    with open(out_path, 'w') as f:
        if filter_by_window_coverage:
            f.write('pos\ttotal\tlow_quality\tlow_quality_std\tlow_mqm\tlow_mqm_std\tlow_cov\thigh_cov\twin_is_overcovered\twin_is_undercovered\twin_has_big_std_cov\n')
        else:
            f.write('pos\ttotal\tlow_quality\tlow_quality_std\tlow_mqm\tlow_mqm_std\tlow_cov\thigh_cov\n')
        for pos in target_pos_list:
            s = [str(pos), str(pos_counts[pos]), str(low_quality_counts[pos]), str(low_quality_std_counts[pos])]
            s.append(str(low_mqm_counts[pos]))
            s.append(str(low_mqm_std_counts[pos]))
            s.append(str(pos_is_undercovered_counts[pos]))
            s.append(str(pos_is_overcovered_counts[pos]))
            if filter_by_window_coverage:
                s.append(str(win_is_overcovered_counts[pos]))
                s.append(str(win_is_undercovered_counts[pos]))
                s.append(str(win_has_big_std_cov_counts[pos]))
            f.write('\t'.join(s))
            f.write('\n')
