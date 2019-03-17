from src.core.constants import data_path
from src.core.data_reading import read_all_variants

path_to_variants = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_no_win_qual_mqm_std3_mqm30_filter_samples_first/'
path_to_sample_ids = path_to_variants + 'samples_filtered.list'
path_to_pos_list = path_to_variants + 'filtered_raw_variants_pos.csv'

out_path = path_to_variants + 'filtered_raw_variants_pos_fixed.csv'


if __name__ == '__main__':
    pos_to_counts = {l.strip(): 0 for l in open(path_to_pos_list).readlines()}
    sample_ids = [l.strip() for l in open(path_to_sample_ids).readlines()]
    sample_to_variants = read_all_variants(path_to_variants, sample_ids)
    for sample_id, variants in sample_to_variants.items():
        for var in variants:
            s = var.split('\t')
            # if s[0] not in pos_to_counts.keys():
            #     print(sample_id)
            pos_to_counts[s[0]] += 1
    with open(out_path, 'w') as fout:
        for pos, c in pos_to_counts.items():
            if c > 0:
                fout.write(pos + '\n')
            else:
                print(pos)
