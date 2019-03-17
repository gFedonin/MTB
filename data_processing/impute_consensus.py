from sklearn.externals.joblib import Parallel, delayed

from src.core.constants import data_path

path_to_sample_variants = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_no_win_m3_genes_filter_samples_first/'
path_to_win_cov = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_win_cov_m3_genes/'
path_to_samples = path_to_sample_variants + 'samples_filtered.list'
path_to_variants = path_to_sample_variants + 'filtered_raw_variants_pos.csv'

out_path = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_win_m3_genes_filter_samples_first'


def read_variants(sample_id):
    variants = []
    uncovered_pos = []
    with open(path_to_sample_variants + sample_id + '.variants', 'r') as f:
        for line in f.readlines():
            variants.append(line.strip())
    with open(path_to_win_cov + sample_id + '.cov', 'r') as f:
        for line in f.readlines()[1:]:
            s = line.strip().split('\t')
            if s[1] == '0':
                uncovered_pos.append(s[0])
    return sample_id, variants, uncovered_pos


if __name__ == '__main__':
    sample_ids = [l.strip() for l in open(path_to_samples).readlines()]
    pos_list = [l.strip() for l in open(path_to_variants).readlines()]
    tasks = Parallel(n_jobs=-1)(delayed(read_variants)(sample_id) for sample_id in sample_ids)
    sample_to_variants = {}
    sample_to_uncovered = {}
    for sample_id, variants, uncovered_pos in tasks:
        sample_to_variants[sample_id] = variants
        sample_to_uncovered[sample_id] = uncovered_pos
    