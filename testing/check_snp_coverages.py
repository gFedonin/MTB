from bisect import bisect_left

from sklearn.externals.joblib import Parallel, delayed

path_to_ids = './data/Full_subset_filtered.txt'
path_to_snps = '/export/data/kchukreev/data/mutations_files/unformated_vcfs/realigned_vcfs/'
path_to_depths = './data/coverages_with_percentiles/5/'
path_to_annotations = './data/AL123456_rev.gff'

out_path = './samples_coverage_anomalies.txt'

qual_thresold = 40
read_len = 150
cov_ratio_threshold = 2
snp_cov_threshold = 15


def get_repeats_coords():
    coords = []
    with open(path_to_annotations, 'r') as f:
        for line in f.readlines():
            s = line.split('\t')
            if s[2] in ('Repeat_region', 'mobile_element'):
                coords.append((int(s[3]), int(s[4])))
    return coords


def read_snps(sample_id, repeats, repeats_starts):
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
                    cov = float(info[i + 3: j])
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
                        snps.append((snp_pos, cov))
    ratio_count = 0
    threshold_count = 0
    triplet_count = 0
    problem_triplets = 0
    for i in range(len(snps)):
        pos_i, cov_i = snps[i]
        for j in range(i + 1, len(snps)):
            pos_j, cov_j = snps[j]
            if pos_j - pos_i < 3:
                triplet_count += 1
                if (cov_i >= snp_cov_threshold) != (cov_j >= snp_cov_threshold):
                    problem_triplets += 1
            if pos_j - pos_i <= read_len:
                if cov_i/cov_j > cov_ratio_threshold or cov_j/cov_i > cov_ratio_threshold:
                    ratio_count += 1
                if (cov_i >= snp_cov_threshold) != (cov_j >= snp_cov_threshold):
                    threshold_count += 1
            else:
                break
    return sample_id, str(len(snps)), str(ratio_count), str(threshold_count), str(triplet_count), str(problem_triplets)


def main():
    sample_ids = [sample_id[:-1] for sample_id in open(path_to_ids, 'r').readlines()]

    filter_coords = get_repeats_coords()
    filter_coords.sort(key=lambda tup: tup[0])
    gene_starts = [x[0] for x in filter_coords]

    tasks = Parallel(n_jobs=-1)(delayed(read_snps)(sample_id, filter_coords, gene_starts) for sample_id in sample_ids)
    with open(out_path, 'w') as f:
        f.write('sample_id\tsnp num\tpairs dist < 150 and cov1/cov2 > 2\tpairs dist < 150 and cov1 >= 15 cov2 < 15' +
                '\tpairs dist < 3\tpairs dist < 3 and cov1 >= 15 cov2 < 15\n')
        for sample_id, snp_num, ratio_count, threshold_count, triplet_count, problem_triplets in tasks:
            f.write(sample_id + '\t' + snp_num + '\t' + ratio_count + '\t' + threshold_count + '\t' + triplet_count +
                    '\t' + problem_triplets + '\n')


if __name__ == '__main__':
    main()
