from math import floor

from core.constants import data_path

# path_to_ids = data_path + 'all_with_pheno.txt'
# path_to_ids = data_path + 'missing.list'
# path_to_ids = data_path + 'list5.txt'
# path_to_ids = data_path + 'coll18/coll5.txt'
# path_to_ids = data_path + 'snps/gatk_before_cortex/combined_ids.list'
# path_to_ids = data_path + 'walker18/missing_err.list'
path_to_ids = data_path + 'snps/combined_raw_variants_mq40_keep_complex_std_names_filtered/samples_filtered.list'
prefix = 'combined'
# path_to_ids = data_path + 'skeza_zero.list'
# path_to_ids = data_path + 'coll18_supp.samples'
# split_array = [0.25, 0.25, 0.25, 0.25]
# split_array = [1, 1]
split_array = [1, 1, 1, 1, 4]
# split_array = [1, 1, 1, 1, 2]

if __name__ == '__main__':
    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    sample_num = len(sample_ids)
    sum_w = sum(split_array)
    counts = [floor(w*sample_num/sum_w) for w in split_array]
    i = 0
    for j in range(len(counts)):
        path = path_to_ids[:path_to_ids.rfind('/') + 1]
        with open(path + prefix + str(j + 1) + '.list', 'w') as f:
            for k in range(counts[j]):
                f.write(sample_ids[i] + '\n')
                i += 1
            if j == len(counts) - 1:
                while i < len(sample_ids):
                    f.write(sample_ids[i] + '\n')
                    i += 1
