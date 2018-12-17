from sklearn.externals.joblib import Parallel, delayed

from src.core.constants import data_path
from src.core.data_reading import read_variants

path_to_var = data_path + 'snps/annotated_with_DR_with_indel_with_pheno_and_snp_mc10/'
path_to_ids = data_path + 'dr_covered_with_pheno_and_snp.txt'
out_path = '../../res/gene_mut.counts'


if __name__ == '__main__':

    parallel = Parallel(n_jobs=-1)

    with open(path_to_ids, 'r') as f:
        sample_ids = [name.strip() for name in f.readlines()]
    print(str(len(sample_ids)) + ' samples')

    tasks = parallel(delayed(read_variants)(path_to_var, sample_id) for sample_id in sample_ids)

    mut_counts = {}
    for name, l in tasks:
        genes_mut = set()
        for snp in l:
            s = snp.split('\t')
            genes_mut.add(s[1])
        for gene in genes_mut:
            c = mut_counts.get(gene)
            if c is None:
                mut_counts[gene] = 1
            else:
                mut_counts[gene] = c + 1
    with open(out_path, 'w') as f:
        for gene, c in mut_counts.items():
            f.write('%s\t%d\n' % (gene, c))