from sklearn.externals.joblib import Parallel, delayed

from src.core.data_reading import read_variants

path_to_var = '../../data/snps/raw_with_DR_with_indel_with_pheno_and_snp_mc10/'
path_to_ids = '../../data/dr_covered_with_pheno_and_snp.txt'
out_path = '../../data/dr_covered_with_pheno_and_snp_unique.counts'

if __name__ == '__main__':
    parallel = Parallel(n_jobs=-1)

    with open(path_to_ids, 'r') as f:
        sample_ids = [name.strip() for name in f.readlines()]
    print(str(len(sample_ids)) + ' samples')

    tasks = parallel(delayed(read_variants)(path_to_var, sample_id) for sample_id in sample_ids)

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