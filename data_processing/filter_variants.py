from os import makedirs
from os.path import exists

import numpy as np
from bitarray import bitarray
from sklearn.externals.joblib import Parallel, delayed

from core.constants import data_path, ref_len
from core.data_reading import read_all_variants

path_to_ids = data_path + 'big_filtered.list'
# path_to_ids = data_path + 'all_with_pheno.txt'
# path_to_snps = data_path + 'snps/gatk_before_cortex/raw_variants_fixed_no_rep_gatk/'
# path_to_snps = data_path + 'snps/gatk_before_cortex/raw_variants_mq40_keep_complex/'
# path_to_snps = data_path + 'snps/pilon/raw_variants/'
path_to_snps = data_path + 'snps/combined_raw_variants_mq40_keep_complex_std_names/'

# out_path = data_path + 'snps/gatk_before_cortex/raw_variants_mq40_keep_complex_filtered/'
# out_path = data_path + 'snps/pilon/raw_variants_mq40_filtered/'
out_path = data_path + 'snps/combined_raw_variants_mq40_keep_complex_std_names_filtered/'

out_path_variants = out_path + 'filtered_raw_variants_pos.csv'
out_path_variants_with_low_cov = out_path + 'low_covered_raw_variants_pos.csv'
out_path_filtered_samples = out_path + 'samples_filtered.list'
out_path_samplies_with_low_coverage = out_path + 'samples_low_covered.list'

# path_to_depths = data_path + 'coverages_shrinked/'
path_to_depths = data_path + 'combined_coverages/'

variant_cov_threshold = 10
sample_cov_threshold = 0.9
thread_num = 32


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


def compute_variants_coverage(sample_id, all_variants_pos_list):
    var_pos_num = len(all_variants_pos_list)
    is_properly_covered = bitarray(var_pos_num)
    coverage = read_cov_file(sample_id)
    for i in range(len(all_variants_pos_list)):
        pos = all_variants_pos_list[i]
        cov = coverage[pos - 1]
        is_properly_covered[i] = cov >= variant_cov_threshold
    return sample_id, is_properly_covered


def read_all_data(sample_ids):

    sample_to_variants = read_all_variants(path_to_snps, sample_ids)
    all_snp_pos = set()
    sample_to_vars = {}
    for sample_id, variants in sample_to_variants.items():
        vars = []
        for var in variants:
            s = var.split('\t')
            pos = int(s[0])
            all_snp_pos.add(pos)
            vars.append((pos, s[1], s[2]))
        sample_to_vars[sample_id] = vars
    all_variants_pos_list = list(all_snp_pos)
    all_variants_pos_list.sort()
    print('variant reading done')


    tasks = Parallel(n_jobs=thread_num)(delayed(compute_variants_coverage)(sample_id, all_variants_pos_list)
                                        for sample_id in sample_ids)
    print('coverages done')
    sample_to_cov = {}

    for sample_id, is_properly_covered in tasks:
        sample_to_cov[sample_id] = is_properly_covered

    print('data reading done')
    print('%d total variants after filtering' % len(all_snp_pos))
    return sample_to_vars, all_variants_pos_list, sample_to_cov, all_snp_pos


def main():
    if not exists(out_path):
        makedirs(out_path)

    sample_ids = [sample_id.strip() for sample_id in open(path_to_ids, 'r').readlines()]
    sample_num = len(sample_ids)

    sample_to_variants, all_variants_pos_list, sample_to_cov, all_filtered_variant_poses = read_all_data(sample_ids)
    variant_pos_num = len(all_variants_pos_list)
    print('reading is done')

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
        print('low covered variants %d' % len(uncovered_variant_pos))
        print('low covered samples %d' % len(uncovered_samples))
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


