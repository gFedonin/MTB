from os import mkdir

import numpy as np
from os.path import exists

from sklearn.externals.joblib import Parallel, delayed

from src.core.constants import ref_len, data_path

# path_to_coverages = '/export/data/kkuleshov/myc/sra/'
path_to_coverages = data_path + 'coverages/'
path_to_ids = data_path + 'all_with_pheno_and_snp.txt'
# out_path = data_path + 'coverages_with_percentiles/t5p5/'
out_path = data_path + 'coverages_shrinked/'
coverage_threshold = 5
percentile = 5


def read_coverage_intervals(sample_id):
    coverages = []
    coverage = np.zeros(ref_len, dtype=int)
    with open(path_to_coverages + sample_id + '/' + sample_id + '_h37rv.depth', 'r') as f:
        is_good_interval = False
        start = '1'
        i = 0
        for line in f.readlines():
            s = line.split()
            cov = int(s[2])
            coverage[i] = cov
            i += 1
            if cov >= coverage_threshold:
                if not is_good_interval:
                    is_good_interval = True
                    start = s[1]
            else:
                if is_good_interval:
                    end = s[1]
                    coverages.append((start, end))
                    is_good_interval = False
        if is_good_interval:
            coverages.append((start, str(ref_len)))
    # coverage.sort()
    lower_bound,  = np.percentile(coverage, percentile, overwrite_input=True)
    upper_bound = np.percentile(coverage, 100 - percentile, overwrite_input=True)
    median = np.percentile(coverage, 50, overwrite_input=True)
    return sample_id, lower_bound, median, upper_bound, coverages


def convert_old(sample_ids):
    tasks = Parallel(n_jobs=-1)(delayed(read_coverage_intervals)(sample_id) for sample_id in sample_ids)
    for task in tasks:
        sample_id, lower_bound, median, upper_bound, coverages = task
        with open(out_path + sample_id + '.coverage', 'w') as f:
            # lower_bound, median, upper_bound, coverages = read_coverage(sample_id)
            f.write('coverage_threshold\tpercentile,%\tlower_percentile_val\tmedian\tupper_pesentile_val\n')
            f.write(str(coverage_threshold) + '\t' + str(percentile) + '\t' + str(lower_bound) + '\t' + str(median) + '\t' + str(upper_bound) + '\n')
            f.write('list of 1-based coords of intervals with coverage >= threshold\n')
            for s, e in coverages:
                f.write(s + '\t' + e + '\n')


def read_depth_file_and_shrink(sample_id):
    with open(path_to_coverages + sample_id + '.depth') as fin:
        with open(out_path + sample_id + '.cov', 'w') as fout:
            for line in fin.readlines():
                s = line.strip().split('\t')
                fout.write(s[2] + '\n')


if __name__ == '__main__':
    if not exists(out_path):
        mkdir(out_path)
    sample_ids = [sample_id.strip() for sample_id in open(path_to_ids, 'r').readlines()]
    for sample_id in sample_ids:
        read_depth_file_and_shrink(sample_id)
