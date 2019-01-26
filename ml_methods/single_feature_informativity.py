import os
from os.path import exists

from sklearn.externals.joblib import Parallel, delayed

from src.core.annotations import read_annotations, CDSType
from src.core.constants import upstream_length
from src.core.data_reading import read_pheno, read_snp_list, read_variants
from src.data_processing.print_mutations_extended import process_variants

path_to_pheno = '../../data/pheno_mc5/'
path_to_var = '../../data/snps/annotated_with_DR_with_indel_with_pheno_and_snp_mc10/'
path_to_ids = '../../data/dr_covered_with_pheno_and_snp_new.txt'

out_path = '../../res/feature_stats/'

use_extended_features = False
merge_all_mut_in_pos = False
merge_all_mut_in_gene = False
filter_non_cds = False
upstream_indel_breakes_gene = True
snp_count_threshold = 3

filter_by_var_list = True
path_to_var_list = path_to_var + 'dr_genes_var_list.csv'


def read_all_pheno():
    drug_to_pheno = {}
    parallel = Parallel(n_jobs=-1)
    for (dirpath, dirnames, filenames) in os.walk(path_to_pheno):
        tasks = parallel(delayed(read_pheno)(path_to_pheno, drug) for drug in [filename[:-6] for filename in filenames])
        for drug, sample_id_to_pheno in tasks:
            drug_to_pheno[drug] = sample_id_to_pheno
    print(str(len(drug_to_pheno.keys())) + ' drugs')
    return drug_to_pheno


def read_all_variants(sample_ids):
    sample_to_variants = {}
    all_snps_set = set()
    if filter_by_var_list:
        full_snp_list = read_snp_list(path_to_var_list)
        tasks = Parallel(n_jobs=-1)(delayed(read_variants)(path_to_var, sample_id, full_snp_list) for sample_id in sample_ids)
    else:
        tasks = Parallel(n_jobs=-1)(delayed(read_variants)(path_to_var, sample_id) for sample_id in sample_ids)
    if use_extended_features:
        cds_list = read_annotations(upstream_length)
        name_to_type = {}
        for cds in cds_list:
            if cds.type != CDSType.upstream:
                name_to_type[cds.name] = cds.type
        tasks1 = Parallel(n_jobs=-1)(delayed(process_variants)(sample_id, snp_list, name_to_type, merge_all_mut_in_pos,
                                               merge_all_mut_in_gene, filter_non_cds, upstream_indel_breakes_gene)
                          for sample_id, snp_list in tasks)
        tasks = tasks1
    for name, l in tasks:
        sample_to_variants[name] = l
        for snp in l:
            all_snps_set.add(snp)
    print('samples\' variants reading done')
    return sample_to_variants, all_snps_set


def informativity(drug, drug_to_pheno, sample_to_variants):
    sample_id_to_pheno = drug_to_pheno[drug]
    variant_to_counts = {}
    for sample_id, pheno in sample_id_to_pheno:
        snp_list = sample_to_variants[sample_id]
        for snp in snp_list:
            counts = variant_to_counts.get(snp)
            if counts is None:
                counts = [0, 0]
                variant_to_counts[snp] = counts
            counts[int(pheno)] += 1
    with open(out_path + drug + '.csv', 'w') as f:
        f.write('type\tname\tpos\tref\talt\ttype\tS\tR\n')
        for var, counts in variant_to_counts.items():
            f.write('%s\t%d\t%d\n' % (var, counts[0], counts[1]))


def main():

    if not exists(out_path):
        os.makedirs(out_path)

    drug_to_pheno = read_all_pheno()

    with open(path_to_ids, 'r') as f:
        sample_ids = [name.strip() for name in f.readlines()]
    print(str(len(sample_ids)) + ' samples')

    sample_to_variants, all_snps_set = read_all_variants(sample_ids)
    for drug in drug_to_pheno.keys():
        informativity(drug, drug_to_pheno, sample_to_variants)


if __name__ == '__main__':
    main()