import gzip
from bisect import bisect_left, bisect_right
from os import makedirs
from os.path import exists

import numpy as np
from bitarray import bitarray
from sklearn.externals.joblib import Parallel, delayed

from core.annotations import CDSType, localize_all_variants, read_annotations
from core.constants import data_path, dr_genes, ref_len, upstream_length

path_to_ids = data_path + 'all_with_pheno.txt'#'test.list'
path_to_snps = data_path + 'snps/freebayes_before_cortex_vcf/'
out_path = data_path + 'snps/raw_freebayes_before_cortex_no_win_qual_mqm_std3_mqm30_no_highcov/'
out_path_variants = out_path + 'filtered_raw_variants_pos.csv'
path_to_annotations = data_path + 'AL123456_rev.gff'

qual_threshold = 40
filter_by_quality_std = True
quality_std_multiplier = 3
filter_by_MQM = True
mqm_threshold = 30
filter_by_MQM_std = True
mqm_std_multiplier = 3
variant_cov_threshold = 10

filter_overcovered_positions = False
coverage_std_multiplier = 3  # set to x for coverage_std_threshold = x * std(coverage)
filter_out_DR_genes = False
snp_only = False


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
    coords.sort(key=lambda tup: tup[0])
    return coords


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
    path_to_vcf = path_to_snps + sample_id + '_h37rv.vcf'
    if not exists(path_to_vcf):
        path_to_vcf = path_to_snps + sample_id + '_h37rv.vcf.gz'
        if exists(path_to_vcf):
            f1 = gzip.open(path_to_vcf, 'rt')
        else:
            return sample_id, variants
    else:
        f1 = open(path_to_vcf, 'r')
    for line in f1.readlines():
        line = line.strip()
        if line[0] != '#':
            tokens = line.split('\t')
            genotype_fields = tokens[-1].split(':')
            gl = genotype_fields[-1].split(',')
            var_index = gl.index('0')
            qual = float(tokens[5])
            info = tokens[-3]
            i = info.index('DP=')
            j = info.index(';', i + 3)
            cov = int(genotype_fields[2].split(',')[var_index])
            if cov < variant_cov_threshold:
                continue
            i = info.index('TYPE=')
            j = info.index(';', i + 5)
            v_type = info[i + 5: j]
            i = info.index('MQM=')
            j = info.index(';', i + 4)
            mqm = float(info[i + 4: j].split(',')[var_index - 1])
            pos = int(tokens[1])
            alt = tokens[4].split(',')[var_index - 1]
            if v_type == 'snp':
                if len(alt) == 1:
                    variants.append((pos, alt, v_type, cov, qual, mqm))
                else:
                    i = info.index('CIGAR=')
                    j = info.index(';', i + 6)
                    cigar = parse_cigar(info[i + 6: j].split(',')[var_index - 1])
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
                cigar = parse_cigar(info[i + 6: j].split(',')[var_index - 1])
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
    f1.close()
    return sample_id, variants


def filter_variants(sample_id, all_variants_pos_list, variants, filter_intervals, pos_to_cds, snp_only):
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
    # is_properly_covered = np.zeros(len(all_variants_pos_list), dtype=int)
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
            filtered_variants.append((pos, alt, type))
    return sample_id, filtered_variants


def filter_sample_list(all_variants_pos_list, sample_to_variants, filter_intervals, pos_to_cds, snp_only,
                    filter_overcovered_positions, coverage_std_multiplier, filter_by_window_coverage,
                    filter_by_window_genes_only):
    sample_to_vars = {}
    for sample_id, variants in sample_to_variants.items():
        sample_id, filtered_variants = filter_variants(sample_id, all_variants_pos_list,
                    variants, filter_intervals, pos_to_cds, snp_only)
        sample_to_vars[sample_id] = filtered_variants
    return sample_to_vars


def read_all_data(sample_ids, filter_overcovered_positions, coverage_std_multiplier, snp_only,
              filter_out_DR_genes):
    """
    Reads and filters variants based on coverages for given list of sample ids

    :param sample_ids: list of sample ids
    :param filter_overcovered_positions: if True, only variants with cov <= coverage_std_multiplier*std_coverage will be kept
    :param coverage_std_multiplier: defines the upper threshold on coverage
    :param snp_only: if True, only variants with type == 'snp' will be kept
    :param filter_out_DR_genes: if True filters out variants inside the list of DR genes
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

    tasks = Parallel(n_jobs=thread_num)(delayed(filter_sample_list)(all_variants_pos_list, sample_vars,
                filter_intervals, pos_to_cds, snp_only) for sample_vars in sample_lists.values())

    sample_to_variants = {}
    all_snp_pos.clear()

    for sample_to_vars in tasks:
        sample_to_variants.update(sample_to_vars)

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
    return sample_to_variants, all_variants_pos_list, all_snp_pos


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
        f.write('filter_short_repeats = ' + str(filter_short_repeats) + '\n')


def main():
    if not exists(out_path):
        makedirs(out_path)
    print_config()

    sample_ids = [sample_id.strip() for sample_id in open(path_to_ids, 'r').readlines()]

    sample_to_variants, all_variants_pos_list, all_filtered_variant_poses = read_all_data(sample_ids,
                                        filter_overcovered_positions, coverage_std_multiplier, snp_only,
                                                        filter_out_DR_genes)

    for sample_id, variants in sample_to_variants.items():
        with open(out_path + sample_id + '.variants', 'w') as f:
            for pos, alt, v_type in variants:
                f.write(str(pos) + '\t' + alt + '\t' + v_type + '\n')
    with open(out_path_variants, 'w') as f:
        for pos in all_filtered_variant_poses:
            f.write(str(pos) + '\n')


if __name__ == '__main__':
    main()
