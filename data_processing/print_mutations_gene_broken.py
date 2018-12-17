import os
from os.path import exists

from sklearn.externals.joblib import Parallel, delayed

from src.core.annotations import read_annotations, CDSType
from src.core.constants import upstream_length
from src.core.data_reading import read_snp_list, read_variants


path_to_var = '../../data/snps/annotated_with_DR_with_indel_with_pheno_and_snp_mc10/'
path_to_ids = '../../data/dr_covered_with_pheno_and_snp.txt'

out_path = '../../data/snps/gene_broken_mc10/'
all_var_out_path = out_path + 'all_var_list.csv'

upstream_indel_breakes_gene = True

thread_num = 32


def process_variants(sample_id, var_list, name_to_type):
    broken_genes = set()
    for var in var_list:
        s = var.split('\t')
        if s[0] == 'upstream' and name_to_type[s[1]] == CDSType.Gene:
            if upstream_indel_breakes_gene and s[5] in ('del', 'ins'):
                broken_genes.add(s[1])
            continue
        if s[0] == 'Gene' and (s[4] == '*' or s[3] == '*' or s[5] == 'FS'):
            broken_genes.add(s[1])
            continue
    res_list = []
    for gene_name in broken_genes:
        res_list.append(gene_name)
    res_list.sort()
    return sample_id, res_list


def main():
    if not exists(out_path):
        os.makedirs(out_path)
    parallel = Parallel(n_jobs=-1)

    with open(path_to_ids, 'r') as f:
        sample_ids = [name.strip() for name in f.readlines()]
    print(str(len(sample_ids)) + ' samples')

    cds_list = read_annotations(upstream_length)
    name_to_type = {}
    for cds in cds_list:
        if cds.type != CDSType.upstream:
            name_to_type[cds.name] = cds.type

    sample_to_variants = {}
    all_snps_set = set()
    tasks = parallel(delayed(read_variants)(path_to_var, sample_id) for sample_id in sample_ids)
    tasks = parallel(delayed(process_variants)(sample_id, snp_list, name_to_type)
                     for sample_id, snp_list in tasks)
    for name, l in tasks:
        sample_to_variants[name] = l
        with open(out_path + name + '.broken_genes', 'w') as f:
            for var in l:
                f.write(var)
                f.write('\n')
                all_snps_set.add(var)
    with open(all_var_out_path, 'w') as f:
        for var in all_snps_set:
            f.write(var)
            f.write('\n')


if __name__ == '__main__':
    main()
