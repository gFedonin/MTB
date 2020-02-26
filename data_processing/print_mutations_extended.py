import os
from os.path import exists

from sklearn.externals.joblib import Parallel, delayed

from core.annotations import read_annotations, CDSType
from core.constants import upstream_length, data_path
from core.data_reading import read_snp_list


# path_to_var = data_path + 'snps/gatk_and_pilon/intersection/'
# path_to_var = data_path + 'snps/gatk_and_pilon/unification/'
# path_to_var = data_path + 'snps/gatk_before_cortex/annotated_pg_NWds10_mq40_keep_complex_filtered_test/'
# path_to_var = data_path + 'snps/pilon/annotated_pg_NWds10_filtered_test/'
path_to_var = data_path + 'snps/annotated_fixed_no_rep_long_del_pg_NWds10_combined/'
# path_to_ids = data_path + 'all_with_pheno.txt'
path_to_ids = data_path + 'snps/combined_raw_variants_mq40_keep_complex_std_names_filtered/samples_filtered.list'

# out_path = data_path + 'snps/gatk_and_pilon/intersection_ext/'
# out_path = data_path + 'snps/gatk_and_pilon/unification_ext/'
# out_path = data_path + 'snps/gatk_before_cortex/annotated_pg_NWds10_mq40_keep_complex_filtered_ext/'
# out_path = data_path + 'snps/pilon/annotated_pg_NWds10_filtered_ext/'
# out_path = data_path + 'snps/pilon/annotated_pg_NWds10_filtered_ext/'
out_path = data_path + 'snps/annotated_fixed_no_rep_long_del_pg_NWds10_combined_extended/'
all_var_out_path = out_path + 'all_var_list.csv'

use_extended_features = True
merge_all_mut_in_pos = True
merge_all_mut_in_gene = True
filter_non_cds = True
drop_pseudogenes = True
drop_hypothetical = False
keep_proteomic_validated_only = False
keep_only_not_hypotetical_or_proteomic_validated = False
upstream_indel_breakes_gene = True

thread_num = 144

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
        if len(s) < 6:
            print(sample_id)
            exit(-1)
        if s[0] == 'Gene' and (s[5] == 'FS'):# s[4] == '*' or s[3] == '*' or
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


def read_variants(path_to_snp, sample_id, name_to_cds, filter_set=None, keep_type=True):
    snp_list = []
    with open(path_to_snp + sample_id + '.variants', 'r') as f:
        for line in f.readlines():
            l = line.strip()
            if not keep_type:
                i = l.rfind('\t')
                l = l[:i]
            if filter_set is not None:
                if l in filter_set:
                    s = l.split('\t')
                    cds = name_to_cds.get(s[1])
                    if filter_non_cds and s[0] == 'non_cds':
                        continue
                    if drop_pseudogenes and cds is not None and cds.is_pseudogene:
                        continue
                    if drop_hypothetical and cds is not None and cds.is_hypothetical:
                        continue
                    if keep_only_not_hypotetical_or_proteomic_validated and cds is not None:
                        if cds.is_hypothetical and not cds.exists_in_proteom:
                            continue
                    snp_list.append(l)
            else:
                s = l.split('\t')
                cds = name_to_cds.get(s[1])
                if filter_non_cds and s[0] == 'non_cds':
                    continue
                if drop_pseudogenes and cds is not None and cds.is_pseudogene:
                    continue
                if drop_hypothetical and cds is not None and cds.is_hypothetical:
                    continue
                if keep_only_not_hypotetical_or_proteomic_validated and cds is not None:
                    if cds.is_hypothetical and not cds.exists_in_proteom:
                        continue
                snp_list.append(l)
    return sample_id, snp_list


def read_all_variants(sample_ids):
    cds_list = read_annotations(upstream_length, filter_by_gene_len=False)
    name_to_cds = {}
    for cds in cds_list:
        name_to_cds[cds.name] = cds
    sample_to_variants = {}
    if filter_by_var_list:
        full_snp_list = read_snp_list(path_to_var_list)
        tasks = Parallel(n_jobs=thread_num, batch_size=len(sample_ids)//thread_num + 1)(delayed(read_variants)(path_to_var, sample_id, name_to_cds, full_snp_list)
                                    for sample_id in sample_ids if exists(path_to_var + sample_id + '.variants'))
    else:
        tasks = Parallel(n_jobs=thread_num, batch_size=len(sample_ids)//thread_num + 1)(delayed(read_variants)(path_to_var, sample_id, name_to_cds)
                                    for sample_id in sample_ids if exists(path_to_var + sample_id + '.variants'))
    for name, snp_list in tasks:
        sample_to_variants[name] = snp_list
    print('samples\' variants reading done')
    return sample_to_variants


def main():
    if not exists(out_path):
        os.makedirs(out_path)
    parallel = Parallel(n_jobs=-1)

    with open(path_to_ids, 'r') as f:
        sample_ids = [name.strip() for name in f.readlines()]
    print(str(len(sample_ids)) + ' samples')

    cds_list = read_annotations(upstream_length, filter_by_gene_len=False)
    name_to_type = {}
    for cds in cds_list:
        if cds.type != CDSType.upstream:
            name_to_type[cds.name] = cds.type

    sample_to_variants = read_all_variants(sample_ids)
    all_snps_set = set()
    tasks = parallel(delayed(process_variants)(sample_id, snp_list, name_to_type, merge_all_mut_in_pos,
                                               merge_all_mut_in_gene, filter_non_cds, upstream_indel_breakes_gene)
                     for sample_id, snp_list in sample_to_variants.items())
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
