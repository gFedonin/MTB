import numpy as np

path_to_coverages = '/export/data/kkuleshov/myc/sra/'
path_to_ids = './data/Full_subset.txt'
out_path = './data/coverages_with_percentiles/5/'
coverage_threshold = 15
persentile = 5
ref_len = 4411532


def read_coverage(sample_id):
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
    coverage.sort()
    lower_bound = np.percentile(coverage, persentile)
    upper_bound = np.percentile(coverage, 100 - persentile)
    median = np.percentile(coverage, 50)
    return lower_bound, median, upper_bound, coverages


if __name__ == '__main__':
    sample_ids = [sample_id[:-1] for sample_id in open(path_to_ids, 'r').readlines()]
    for sample_id in sample_ids:
        with open(out_path + sample_id + '.coverage', 'w') as f:
            lower_bound, median, upper_bound, coverages = read_coverage(sample_id)
            f.write('coverage_threshold\tpercentile,%\tlower_percentile_val\tmedian\tupper_pesentile_val\n')
            f.write(str(coverage_threshold) + '\t' + str(persentile) + '\t' + str(lower_bound) + '\t' + str(median) + '\t' + str(upper_bound) + '\n')
            f.write('list of 1-based coords of intervals with coverage >= threshold\n')
            for s, e in coverages:
                f.write(s + '\t' + e + '\n')
