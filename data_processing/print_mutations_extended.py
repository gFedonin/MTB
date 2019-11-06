import os
from os.path import exists

from sklearn.externals.joblib import Parallel, delayed

from core.annotations import read_annotations, CDSType
from core.constants import upstream_length
from core.data_reading import read_snp_list, read_variants


path_to_var = '../../data/snps/annotated_with_DR_with_indel_with_pheno_and_snp_mc10/'
path_to_ids = '../../data/dr_covered_with_pheno_and_snp.txt'

out_path = '../../data/snps/annotated_with_indel_mc10_extended/'
all_var_out_path = out_path + 'all_var_list.csv'

merge_all_mut_in_pos = True
merge_all_mut_in_gene = True
filter_non_cds = False
upstream_indel_breakes_gene = True

thread_num = 32

filter_by_var_list = False
path_to_var_list = path_to_var + 'dr_genes_var_list.csv'


def process_variants(sample_id, var_list, name_to_type, merge_all_mut_in_pos, merge_all_mut_in_gene, filter_non_cds,
                     upstream_indel_breakes_gene):
    res = set()
    broken_genes = set()
    changed_genes = set()
    upstreams = []

    for var in var_list:
        s = var.split('\t')
        if filter_non_cds and s[0] == 'non_cds':
            continue
        if s[0] == 'upstream':
            upstreams.append(s)
            continue
        if s[0] == 'Gene' and (s[4] == '*' or s[3] == '*' or s[5] == 'FS'):
            broken_genes.add(s[1])
            continue
        res.add(var)
        if s[0] != 'non_cds' and merge_all_mut_in_gene:
            changed_genes.add(s[1])
        if merge_all_mut_in_pos:
            if s[5] != 'ins':
                res.add('\t'.join(s[:-3]))
    for s in upstreams:
        res.add('\t'.join(s))
        if upstream_indel_breakes_gene and s[5] in ('del', 'ins'):
            broken_genes.add(s[1])
        else:
            changed_genes.add(s[1])

    res_list = []
    for var in res:
        s = var.split('\t')
        if s[1] not in broken_genes:
            res_list.append(var)
    for gene_name in broken_genes:
        res_list.append(name_to_type[gene_name].name + '\t' + gene_name + '\tbroken')
    for gene_name in changed_genes:
        if gene_name not in broken_genes:
            res_list.append(name_to_type[gene_name].name + '\t' + gene_name + '\tchanged')
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
    if filter_by_var_list:
        full_snp_list = read_snp_list(path_to_var_list)
        tasks = parallel(delayed(read_variants)(path_to_var, sample_id, full_snp_list) for sample_id in sample_ids)
    else:
        tasks = parallel(delayed(read_variants)(path_to_var, sample_id) for sample_id in sample_ids)
    tasks = parallel(delayed(process_variants)(sample_id, snp_list, name_to_type, merge_all_mut_in_pos,
                                               merge_all_mut_in_gene, filter_non_cds, upstream_indel_breakes_gene)
                     for sample_id, snp_list in tasks)
    for name, l in tasks:
        sample_to_variants[name] = l
        with open(out_path + name + '.variants', 'w') as f:
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
