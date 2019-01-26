from sklearn.externals.joblib import Parallel, delayed

from src.core.data_reading import read_variants

path_to_var = '../../data/snps/annotated_with_DR_with_indel_with_pheno_and_snp_mc10/'
path_to_ids = '../../data/dr_covered_with_pheno_and_snp_new.txt'

out_path = '../../res/gene_mut_counts.csv'


def read_all_variants(sample_ids):
    sample_to_variants = {}
    all_snps_set = set()
    tasks = Parallel(n_jobs=-1)(delayed(read_variants)(path_to_var, sample_id) for sample_id in sample_ids)
    for name, l in tasks:
        sample_to_variants[name] = l
        for snp in l:
            all_snps_set.add(snp)
    print('samples\' variants reading done')
    return sample_to_variants, all_snps_set


def print_mut_counts(sample_to_variants):
    snp_to_count = {}
    for sample_id, snp_list in sample_to_variants.items():
        for snp in snp_list:
            c = snp_to_count.get(snp)
            if c is None:
                snp_to_count[snp] = 1
            else:
                snp_to_count[snp] = c + 1
    genes_to_counts = {}
    genes_with_upstreams_to_counts = {}
    genes_with_upstreams_mut_to_counts = {}
    for snp, count in snp_to_count.items():
        s = snp.split('\t')
        if s[0] == 'Gene':
            c = genes_to_counts.get(s[1])
            if c is None:
                genes_to_counts[s[1]] = count
            else:
                genes_to_counts[s[1]] = c + count
            c = genes_with_upstreams_to_counts.get(s[1])
            if c is None:
                genes_with_upstreams_to_counts[s[1]] = count
            else:
                genes_with_upstreams_to_counts[s[1]] = c + count
            mut_to_count = genes_with_upstreams_mut_to_counts.get(s[1])
            if mut_to_count is None:
                mut_to_count = {}
                genes_with_upstreams_mut_to_counts[s[1]] = mut_to_count
            mut_to_count[snp] = count
        elif s[0] == 'upstream':
            c = genes_with_upstreams_to_counts.get(s[1])
            if c is None:
                genes_with_upstreams_to_counts[s[1]] = count
            else:
                genes_with_upstreams_to_counts[s[1]] = c + count
            mut_to_count = genes_with_upstreams_mut_to_counts.get(s[1])
            if mut_to_count is None:
                mut_to_count = {}
                genes_with_upstreams_mut_to_counts[s[1]] = mut_to_count
            mut_to_count[snp] = count
    with open(out_path, 'w') as f:
        f.write('gene_name\tmut_count_in_gene\tmut_count_in_gene_with_upstream\tunique_mutations_count\tmost_frequent_mut_count_in_gene_with_upstream\n')
        for gene, c in genes_to_counts.items():
            mut_to_count = genes_with_upstreams_mut_to_counts[gene]
            max_count = max(mut_to_count.values())
            unique = len(mut_to_count)
            f.write('%s\t%d\t%d\t%d\t%d\n' % (gene, c, genes_with_upstreams_to_counts[gene], unique, max_count))


def main():
    with open(path_to_ids, 'r') as f:
        sample_ids = [name.strip() for name in f.readlines()]
    print(str(len(sample_ids)) + ' samples')
    sample_to_variants, all_snps_set = read_all_variants(sample_ids)
    print_mut_counts(sample_to_variants)


if __name__ == '__main__':
    main()