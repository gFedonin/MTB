from os import makedirs
import numpy as np
from bisect import bisect_left, bisect_right
from os.path import exists
from bitarray import bitarray
from sklearn.externals.joblib import Parallel, delayed

from src.core.annotations import read_annotations, localize_all_variants, CDSType
from src.core.constants import upstream_length, dr_genes, data_path, ref_len

path_to_ids = data_path + 'all_with_pheno_and_snp.txt'#'test.list'
path_to_snps = data_path + 'snps/realigned_vcfs_indels_all/'
out_path = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_no_win_qual_mqm_std3_mqm30_no_highcov_str10/'
out_path_win_cov = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_win_cov_m3_filter_samples_first/'
out_path_variants = out_path + 'filtered_raw_variants_pos.csv'
out_path_variants_with_low_cov = out_path + 'low_covered_raw_variants_pos.csv'
out_path_filtered_samples = out_path + 'samples_filtered.list'
out_path_samplies_with_low_coverage = out_path + 'samples_low_covered.list'
path_to_depths = data_path + 'coverages_shrinked/'#''/export/data/kkuleshov/myc/sra/'
path_to_annotations = data_path + 'AL123456_rev.gff'

filter_short_repeats = True
path_to_short_tandem_repeats = data_path + 'h37rv.fasta.2.7.7.80.10.20.10.dat'
qual_threshold = 40
filter_by_quality_std = True
quality_std_multiplier = 3
filter_by_MQM = True
mqm_threshold = 30
filter_by_MQM_std = True
mqm_std_multiplier = 3
variant_cov_threshold = 10
sample_cov_threshold = 0.9

filter_overcovered_positions = False
coverage_std_multiplier = 3  # set to x for coverage_std_threshold = x * std(coverage)
filter_out_DR_genes = False
snp_only = False
filter_by_window_coverage = False
filter_by_window_genes_only = False
print_coverage_statistics = False
window_len_gene = 400
window_len_nc = 50
window_min_cov_threshold = 10
window_max_cov_threshold_coef = 3  # set to x for window_max_cov_threshold = x * std(coverage)
window_std_cov_threshold = 3

thread_num = 32


def get_intervals_to_filter_out(filter_out_DR_genes, filter_short_repeats):
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
            elif s[2] in ('Repeat_region', 'mobile_element'):
                coords.append((int(s[3]), int(s[4])))
    if filter_short_repeats:
        with open(path_to_short_tandem_repeats) as f:
            for line in f.readlines()[15:]:
                s = line.split()
                coords.append((int(s[0]), int(s[1])))
    coords.sort(key=lambda tup: tup[0])
    return coords


def read_depth_file(sample_id):
    """
    Reads raw coverage .depth file for given sample id and computes coverage threshold for given percentile

    :param sample_id: sample id to read coverage
    :return: numpy array of length ref_len with coverages of positions, zero-based, upper percentile coverage threshold
    """
    coverage = np.zeros(ref_len, dtype=int)
    i = 0
    with open(path_to_depths + sample_id + '.depth') as f:
        for line in f.readlines():
            s = line.strip().split('\t')
            coverage[i] = int(s[2])
            i += 1
    return coverage


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


def window_is_covered(pos, coverage, mean_coverage, std_coverage, cds, filter_by_window_genes_only):
    """
    Check coverage based conditions in window centered on given position
    :param pos: genome coordinate of interest
    :param coverage: array representing sample coverage
    :param mean_coverage: mean(coverage)
    :param std_coverage: std(coverage)
    :param cds: CDS to which this position belongs
    :param filter_by_window_genes_only: if True skip all positions which are not inside protein coding genes
    :return: is window coverage < than window_min_cov_threshold, is window coverage >
            window_max_cov_threshold_coef*std_coverage, is std of window coverage > window_std_cov_threshold*std_coverage
    """
    if cds is None or cds.type != CDSType.Gene:
        if filter_by_window_genes_only:
            return False, False, False
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
    # if mean_cov < window_min_cov_threshold:
    #     return False
    # if mean_cov > mean_coverage + window_max_cov_threshold_coef*std_coverage:
    #     return False
    # if stats.levene(window, coverage)[1] < window_std_difference_test_pval:
    #     return False
    # if std_cov > window_std_cov_threshold*std_coverage:
    #     return False
    return undercovered, overcovered, bigstdcov


def filter_variants(sample_id, all_variants_pos_list, variants, filter_intervals, pos_to_cds, snp_only,
                    filter_overcovered_positions, coverage_std_multiplier, filter_by_window_coverage,
                    filter_by_window_genes_only):
    """
    Filters list of variants by given criteria:
    - variant coverage should be >= snp_cov_threshold
    - if filter_by_coverage_std set to True, variant coverage should be <= coverage_std_multiplier*std_coverage
    - variant should not be in filter_intervals

    :param variants: list of variants to filter
    :param filter_intervals: variant in side these intervals will be filtered out
    :param coverage: sample coverage
    :param snp_only: if True, only variants with type == 'snp' will be kept
    :param filter_overcovered_positions: if True, only variants with cov <= coverage_std_multiplier*std_coverage will be kept
    :param coverage_std_multiplier: defines the upper threshold on coverage
    :param std_coverage: stdDev of sample coverage
    :param filter_by_window_coverage: if True, apply window based filters
    :param filter_by_window_genes_only: if True, apply window based filters only to protein coding genes
    :return: sample_id, list of tuples (pos, alt, type) representing variants for this sample, bitarray of length
    """
    filtered_variants = []
    var_pos_num = len(all_variants_pos_list)
    is_properly_covered = bitarray(var_pos_num)
    pos_is_undercovered = bitarray(var_pos_num)
    pos_is_overcovered = bitarray(var_pos_num)
    win_is_undercovered = bitarray(var_pos_num)
    win_is_overcovered = bitarray(var_pos_num)
    win_has_big_std_cov = bitarray(var_pos_num)
    # is_properly_covered = np.zeros(len(all_variants_pos_list), dtype=int)
    coverage = read_cov_file(sample_id)
    std_coverage = np.std(coverage)
    mean_coverage = np.mean(coverage)
    qualities = np.zeros(len(variants))
    mqm_list = np.zeros(len(variants))
    filter_intervals_starts = [x[0] for x in filter_intervals]
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
        cds = pos_to_cds.get(pos)
        if qual >= qual_threshold and cov >= variant_cov_threshold:
            if filter_overcovered_positions and cov > mean_coverage + coverage_std_multiplier*std_coverage:
                continue
            if filter_by_MQM and mqm < mqm_threshold:
                continue
            if snp_only and type != 'snp':
                continue
            if filter_by_quality_std and qual < mean_quality - quality_std_multiplier*std_quality:
                continue
            if filter_by_MQM_std and mqm < mean_mqm - mqm_std_multiplier*std_mqm:
                continue
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
            if inside_filtered_interval:
                continue
            if filter_by_window_coverage:
                win_undercovered, win_overcovered, win_bigstdcov = window_is_covered(pos, coverage, mean_coverage,
                                                                        std_coverage, cds, filter_by_window_genes_only)
                if (win_overcovered or win_undercovered or win_bigstdcov):
                    continue
            filtered_variants.append((pos, alt, type))
    for i in range(len(all_variants_pos_list)):
        pos = all_variants_pos_list[i]
        cov = coverage[pos - 1]
        cds = pos_to_cds.get(pos)
        is_properly_covered[i] = True
        if cov < variant_cov_threshold:
            is_properly_covered[i] = False
            pos_is_undercovered[i] = True
        if filter_overcovered_positions and cov > mean_coverage + coverage_std_multiplier*std_coverage:
            is_properly_covered[i] = False
            pos_is_overcovered[i] = True
        if filter_by_window_coverage:
            win_is_undercovered[i], win_is_overcovered[i], win_has_big_std_cov[i] = \
                window_is_covered(pos, coverage, mean_coverage, std_coverage, cds, filter_by_window_genes_only)
        if filter_by_window_coverage and (win_is_undercovered[i] or win_is_overcovered[i] or win_has_big_std_cov[i]):
            is_properly_covered[i] = False
        # is_properly_covered[i] = 1
    if print_coverage_statistics:
        with open(out_path_win_cov + sample_id + '.cov', 'w') as f:
            f.write('pos\tis_properly_covered\tpos_is_undercovered\tpos_is_overcovered\twin_is_undercovered\t'
                    'win_is_overcovered\twin_has_big_std_cov\n')
            for i in range(len(all_variants_pos_list)):
                pos = all_variants_pos_list[i]
                f.write(str(pos) + '\t')
                if is_properly_covered[i]:
                # if is_properly_covered[i] == 1:
                    f.write('1\t')
                else:
                    f.write('0\t')
                if pos_is_undercovered[i]:
                    f.write('1\t')
                else:
                    f.write('0\t')
                if pos_is_overcovered[i]:
                    f.write('1\t')
                else:
                    f.write('0\t')
                if win_is_undercovered[i]:
                    f.write('1\t')
                else:
                    f.write('0\t')
                if win_is_overcovered[i]:
                    f.write('1\t')
                else:
                    f.write('0\t')
                if win_has_big_std_cov[i]:
                    f.write('1\n')
                else:
                    f.write('0\n')
    return sample_id, filtered_variants, is_properly_covered


def filter_sample_list(all_variants_pos_list, sample_to_variants, filter_intervals, pos_to_cds, snp_only,
                    filter_overcovered_positions, coverage_std_multiplier, filter_by_window_coverage,
                    filter_by_window_genes_only):
    sample_to_cov = {}
    sample_to_vars = {}
    for sample_id, variants in sample_to_variants.items():
        sample_id, filtered_variants, is_properly_covered = filter_variants(sample_id, all_variants_pos_list,
                    variants, filter_intervals, pos_to_cds, snp_only,
                    filter_overcovered_positions, coverage_std_multiplier, filter_by_window_coverage,
                    filter_by_window_genes_only)
        sample_to_vars[sample_id] = filtered_variants
        sample_to_cov[sample_id] = is_properly_covered
    return sample_to_vars, sample_to_cov


def read_all_data(sample_ids, filter_overcovered_positions, coverage_std_multiplier, snp_only,
              filter_out_DR_genes, filter_by_window_coverage, filter_by_window_genes_only):
    """
    Reads and filters variants based on coverages for given list of sample ids

    :param sample_ids: list of sample ids
    :param filter_overcovered_positions: if True, only variants with cov <= coverage_std_multiplier*std_coverage will be kept
    :param coverage_std_multiplier: defines the upper threshold on coverage
    :param snp_only: if True, only variants with type == 'snp' will be kept
    :param filter_out_DR_genes: if True filters out variants inside the list of DR genes
    :param filter_by_window_coverage: if True, apply window based filters
    :param filter_by_window_genes_only: if True, apply window based filters only to protein coding genes
    :return: sample to variants map, list of all variant positions, sample to coverage map
    """

    filter_intervals = get_intervals_to_filter_out(filter_out_DR_genes, filter_short_repeats)
    sample_to_variants = {}
    tasks = Parallel(n_jobs=thread_num)(delayed(read_vcf_file)(sample_id) for sample_id in sample_ids)
    all_snp_pos = set()
    for sample_id, variants in tasks:
        sample_to_variants[sample_id] = variants
        for pos, alt, type, cov, qual, mqm in variants:
            all_snp_pos.add(pos)
    all_variants_pos_list = list(all_snp_pos)
    all_variants_pos_list.sort()

    cds_list = read_annotations(upstream_length)
    pos_to_cds = localize_all_variants(all_variants_pos_list, cds_list)

    sample_lists = {i: {} for i in range(thread_num)}
    for i in range(len(sample_ids)):
        sample_id = sample_ids[i]
        sample_to_vars = sample_lists[i%thread_num]
        sample_to_vars[sample_id] = sample_to_variants[sample_id]

    # tasks = Parallel(n_jobs=thread_num, batch_size=len(sample_ids)//thread_num + 1)(delayed(filter_variants)(sample_id, all_variants_pos_list, variants,
    #             filter_intervals, pos_to_cds, snp_only, filter_overcovered_positions, coverage_std_multiplier,
    #             filter_by_window_coverage, filter_by_window_genes_only) for sample_id, variants in sample_to_variants.items())
    tasks = Parallel(n_jobs=thread_num)(delayed(filter_sample_list)(all_variants_pos_list, sample_vars,
                filter_intervals, pos_to_cds, snp_only, filter_overcovered_positions, coverage_std_multiplier,
                filter_by_window_coverage, filter_by_window_genes_only) for sample_vars in sample_lists.values())

    sample_to_cov = {}
    sample_to_variants = {}
    all_snp_pos.clear()

    for sample_to_vars, sample_to_c in tasks:
        sample_to_variants.update(sample_to_vars)
        sample_to_cov.update(sample_to_c)

    # for sample_id, filtered_variants, is_properly_covered in tasks:
    #     sample_to_variants[sample_id] = filtered_variants
    #     sample_to_cov[sample_id] = is_properly_covered
    #     for pos, alt, type in filtered_variants:
    #         all_snp_pos.add(pos)
    for sample_id, filtered_variants in sample_to_variants.items():
        for pos, alt, type in filtered_variants:
            all_snp_pos.add(pos)
    print('data reading done')
    print('%d total variants after filtering' % len(all_snp_pos))
    return sample_to_variants, all_variants_pos_list, sample_to_cov, all_snp_pos


def print_config():
    with open(out_path + 'config.txt', 'w') as f:
        f.write('qual_thresold = ' + str(qual_threshold) + '\n')
        f.write('snp_cov_threshold = ' + str(variant_cov_threshold) + '\n')
        f.write('sample_cov_threshold = ' + str(sample_cov_threshold) + '\n')
        f.write('gene_upstream_threshold = ' + str(upstream_length) + '\n')
        f.write('filter_DR_genes = ' + str(filter_out_DR_genes) + '\n')
        f.write('snp only = ' + str(snp_only) + '\n')
        f.write('filter_overcovered_positions = ' + str(filter_overcovered_positions) + '\n')
        f.write('coverage_std_multiplier = ' + str(coverage_std_multiplier) + '\n')
        f.write('filter_by_MQM = ' + str(filter_by_MQM) + '\n')
        f.write('mqm_threshold = ' + str(mqm_threshold) + '\n')
        f.write('filter_by_quality_std = ' + str(filter_by_quality_std) + '\n')
        f.write('quality_std_multiplier = ' + str(quality_std_multiplier) + '\n')
        f.write('filter_by_MQM_std = ' + str(filter_by_MQM_std) + '\n')
        f.write('mqm_std_multiplier = ' + str(mqm_std_multiplier) + '\n')
        f.write('filter_by_window_coverage = ' + str(filter_by_window_coverage) + '\n')
        f.write('window_len_gene = ' + str(window_len_gene) + '\n')
        f.write('window_len_nc = ' + str(window_len_nc) + '\n')
        f.write('window_min_cov_threshold = ' + str(window_min_cov_threshold) + '\n')
        f.write('window_max_cov_threshold_coef = ' + str(window_max_cov_threshold_coef) + '\n')
        f.write('window_std_cov_threshold = ' + str(window_std_cov_threshold) + '\n')
        f.write('filter_by_window_genes_only = ' + str(filter_by_window_genes_only) + '\n')
        f.write('filter_short_repeats = ' + str(filter_short_repeats) + '\n')


def main():
    if not exists(out_path):
        makedirs(out_path)
    if print_coverage_statistics:
        if not exists(out_path_win_cov):
            makedirs(out_path_win_cov)
    print_config()

    sample_ids = [sample_id.strip() for sample_id in open(path_to_ids, 'r').readlines()]
    sample_num = len(sample_ids)

    sample_to_variants, all_variants_pos_list, sample_to_cov, all_filtered_variant_poses = read_all_data(sample_ids,
                                        filter_overcovered_positions, coverage_std_multiplier, snp_only,
                                                        filter_out_DR_genes, filter_by_window_coverage,
                                                                             filter_by_window_genes_only)
    variant_pos_num = len(all_variants_pos_list)

    uncovered_variant_pos = set(pos for pos in all_variants_pos_list if pos not in all_filtered_variant_poses)
    uncovered_samples = set()

    pos_counts = {}
    for sample_id, variants in sample_to_variants.items():
        for pos, alt, v_type in variants:
            c = pos_counts.get(pos)
            if c is None:
                pos_counts[pos] = 1
            else:
                pos_counts[pos] = c + 1

    while len(uncovered_samples) < sample_num and len(uncovered_variant_pos) < variant_pos_num:
        filtered_variant_pos_num = variant_pos_num - len(uncovered_variant_pos)
        if filtered_variant_pos_num == 0:
            break
        uncovered_samples_added = False
        for sample_id, is_properly_covered in sample_to_cov.items():
            if sample_id in uncovered_samples:
                continue
            properly_covered_pos_num = 0
            for i in range(variant_pos_num):
                pos = all_variants_pos_list[i]
                if pos not in uncovered_variant_pos and is_properly_covered[i]:
                # if pos not in uncovered_variant_pos and is_properly_covered[i] == 1:
                    properly_covered_pos_num += 1
            if properly_covered_pos_num/filtered_variant_pos_num < sample_cov_threshold:
                uncovered_samples.add(sample_id)
                for pos, alt, v_type in sample_to_variants[sample_id]:
                    c = pos_counts[pos] - 1
                    if c == 0:
                        uncovered_variant_pos.add(pos)
                        filtered_variant_pos_num -= 1
                    else:
                        pos_counts[pos] = c
                uncovered_samples_added = True
        snp_cov = np.zeros(variant_pos_num, dtype=int)
        for sample_id, is_properly_covered in sample_to_cov.items():
            if sample_id in uncovered_samples:
                continue
            for i in range(variant_pos_num):
                pos = all_variants_pos_list[i]
                if pos not in uncovered_variant_pos and is_properly_covered[i]:
                # if pos not in uncovered_variant_pos and is_properly_covered[i] == 1:
                    snp_cov[i] += 1
        uncovered_pos_added = False
        filtered_sample_num = sample_num - len(uncovered_samples)
        for i in range(variant_pos_num):
            pos = all_variants_pos_list[i]
            if pos not in uncovered_variant_pos and snp_cov[i]/filtered_sample_num < sample_cov_threshold:
                uncovered_variant_pos.add(pos)
                uncovered_pos_added = True
        if not uncovered_pos_added and not uncovered_samples_added:
            break

    with open(out_path_samplies_with_low_coverage, 'w') as f:
        for sample_id in uncovered_samples:
            f.write(sample_id + '\n')
    with open(out_path_filtered_samples, 'w') as f:
        for sample_id in sample_ids:
            if sample_id not in uncovered_samples:
                f.write(sample_id + '\n')
    for sample_id, variants in sample_to_variants.items():
        if sample_id not in uncovered_samples:
            with open(out_path + sample_id + '.variants', 'w') as f:
                for pos, alt, v_type in variants:
                    if pos not in uncovered_variant_pos:
                        f.write(str(pos) + '\t' + alt + '\t' + v_type + '\n')
    with open(out_path_variants, 'w') as f:
        for pos in all_variants_pos_list:
            if pos not in uncovered_variant_pos:
                f.write(str(pos) + '\n')
    with open(out_path_variants_with_low_cov, 'w') as f:
        for pos in all_variants_pos_list:
            if pos in uncovered_variant_pos:
                f.write(str(pos) + '\n')


if __name__ == '__main__':
    main()
