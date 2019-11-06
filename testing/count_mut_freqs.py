from sklearn.externals.joblib import Parallel, delayed

from src.core.constants import data_path
from src.core.data_reading import read_variants

# path_to_var = data_path + 'snps/gatk_before_cortex/raw_variants_std2/'
path_to_var = data_path + 'snps/gatk_before_cortex/raw_variants_std3/'
path_to_ids = data_path + 'all_with_pheno.txt'
# out_path = '../../res/ml_log_mc3_gatk_before_std2/var.counts'
out_path = '../../res/ml_log_mc3_gatk_before_std3/var.counts'

def count_unique():
    sample_ids = [name.strip() for name in open(path_to_ids).readlines()]
    print(str(len(sample_ids)) + ' samples')

    tasks = Parallel(n_jobs=-1)(delayed(read_variants)(path_to_var, sample_id)
                     for sample_id in sample_ids)

    mut_counts = {}
    for name, l in tasks:
        for snp in l:
            c = mut_counts.get(snp)
            if c is None:
                mut_counts[snp] = 1
            else:
                mut_counts[snp] = c + 1
    with open(out_path, 'w') as f:
        for sample_id, l in tasks:
            unique_snps = 0
            for snp in l:
                if mut_counts[snp] == 1:
                    unique_snps += 1
            f.write('%s\t%d\t%d\n' % (sample_id, unique_snps, len(l)))


def print_all_counts():
    sample_ids = [name.strip() for name in open(path_to_ids).readlines()]
    print(str(len(sample_ids)) + ' samples')

    tasks = Parallel(n_jobs=-1)(delayed(read_variants)(path_to_var, sample_id)
                     for sample_id in sample_ids)

    mut_counts = {}
    for sample_id, l in tasks:
        for snp in l:
            c = mut_counts.get(snp)
            if c is None:
                mut_counts[snp] = 1
            else:
                mut_counts[snp] = c + 1
    counts = [(snp, c) for snp, c in mut_counts.items()]
    counts.sort(key=lambda x: int(x[0].split()[0]))
    with open(out_path, 'w') as f:
        for snp, c in counts:
            f.write('%s\t%d\n' % (snp, c))


if __name__ == '__main__':
    # count_unique()
    print_all_counts()
