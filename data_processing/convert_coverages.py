from os import mkdir, listdir

import numpy as np
from os.path import exists

from sklearn.externals.joblib import Parallel, delayed

from core.constants import ref_len, data_path

# data_set = 'walker18'
# data_set = 'coll18'
data_set = 'missing'
# path_to_coverages = '/export/data/kkuleshov/myc/sra/'
# path_to_coverages = data_path + 'coverages_rm_dup/'
path_to_coverages = data_path + data_set + '/coverages/'
# path_to_ids = data_path + 'all_with_pheno.txt'
path_to_ids = data_path + data_set + '/' + data_set + '_trimmed.list'
# path_to_ids = data_path + data_set + '/' + data_set + '_new.samples'
# out_path = data_path + 'coverages_with_percentiles/t5p5/'
# out_path = data_path + 'coverages_shrinked_rm_dup/'
out_path = data_path + data_set + '/coverages_shrinked/'
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
    if not exists(out_path + sample_id + '.cov'):
        with open(path_to_coverages + sample_id + '.depth') as fin:
            with open(out_path + sample_id + '.cov', 'w') as fout:
                for line in fin.readlines():
                    s = line.strip().split('\t')
                    fout.write(s[2] + '\n')
        return 1


# path_to_mapping = data_path + data_set + '/err_to_samea.csv'
path_to_mapping = data_path + data_set + '/' + data_set + '_id_mapping.csv'


def read_depth_file_shrink_and_rename(sample_id, samea_id):
    if not exists(out_path + samea_id + '.cov'):
        with open(path_to_coverages + sample_id + '.depth') as fin:
            with open(out_path + samea_id + '.cov', 'w') as fout:
                for line in fin.readlines():
                    s = line.strip().split('\t')
                    fout.write(s[2] + '\n')
        return 1


def shrink_and_rename_all():
    if not exists(out_path):
        mkdir(out_path)
    err_to_samea = {}
    for l in open(path_to_mapping).readlines():
        s = l.strip().split('\t')
        err_to_samea[s[0]] = s[1]
    sample_ids = [fname[:fname.index('.')] for fname in listdir(path_to_coverages)]
    tasks = Parallel(n_jobs=-1)(delayed(read_depth_file_shrink_and_rename)(sample_id, err_to_samea[sample_id])
                                for sample_id in sample_ids if not exists(out_path + err_to_samea[sample_id] + '.cov'))
    c = 0
    for task in tasks:
        c += task
    print(str(c))
    # for fname in listdir(path_to_coverages):
    #     sample_id = fname[:fname.index('.')]
    #     read_depth_file_shrink_and_rename(sample_id, err_to_samea[sample_id])


def shrink_all():
    if not exists(out_path):
        mkdir(out_path)
    sample_ids = [fname[:fname.index('.')] for fname in listdir(path_to_coverages)]
    tasks = Parallel(n_jobs=-1)(delayed(read_depth_file_and_shrink)(sample_id)
                                for sample_id in sample_ids if not exists(out_path + sample_id + '.cov'))
    c = 0
    for task in tasks:
        c += task
    print(str(c))
    # for sample_id in sample_ids:
    #     read_depth_file_and_shrink(sample_id)


def find_duplicates():
    err_to_samea = {}
    for l in open(path_to_mapping).readlines():
        s = l.strip().split('\t')
        err_to_samea[s[0]] = s[1]
    sample_ids = [sample_id.strip() for sample_id in open(path_to_ids, 'r').readlines()]
    samea_to_ids = {}
    for sample_id in sample_ids:
        for info, samea in err_to_samea:
            if sample_id in info:
                ids = samea_to_ids.get(samea)
                if ids is None:
                    ids = [sample_id]
                    samea_to_ids[samea] = ids
                else:
                    ids.append(sample_id)
                break
    for samea, ids in samea_to_ids.items():
        if len(ids) > 1:
            print(samea + '\t' + ','.join(ids))


if __name__ == '__main__':
    # shrink_all()
    shrink_and_rename_all()
    # find_duplicates()