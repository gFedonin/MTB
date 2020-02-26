from bisect import bisect_left
from os import listdir, makedirs
from os.path import exists

from sklearn.externals.joblib import Parallel, delayed

from core.constants import data_path, ref_len, upstream_length
from core.data_reading import read_variants, read_pheno
from core.annotations import read_annotations, localize_all_variants

path_to_ids = data_path + 'all_with_pheno.txt'
# path_to_ids = data_path + 'snps/gatk_pilon_intersection.list'
# path_to_ids = data_path + 'snps/gatk_pilon_old_intersection.list'
# path_to_ids = data_path + 'debug2.list'
# path_to_variants1 = data_path + 'snps/freebayes_before_cortex/raw_no_win_qual_mqm_std3_mqm30_no_highcov_ld/'
# name1 = 'freebayes_before'
# path_to_variants1 = data_path + 'snps/gatk_after_cortex/raw_variants_std2_ld/'
# name1 = 'gatk_after_std2_ld'
# path_to_variants1 = data_path + 'snps/raw_no_win_qual_mqm_std3_mqm30_no_highcov_bam_filtered_ld/'
# name1 = 'freebayes_old_bam_filtered'
# path_to_variants1 = data_path + 'snps/raw_no_win_qual_mqm_std3_mqm30_no_highcov_ld/'
# name1 = 'freebayes_old'
# path_to_variants1 = data_path + 'snps/raw_no_win_qual_mqm_std3_mqm30_no_highcov_no_repeats_no_cov_ld/'
# name1 = 'freebayes_old_no_rep'
# path_to_variants1 = data_path + 'snps/skesa_mummer_raw_ld_mum4/'
# name1 = 'skesa_mummer_raw_ld_mum4'
# path_to_variants1 = data_path + 'snps/gatk_before_cortex/raw_variants_ld/'
# name1 = 'gatk_before_ld'
# path_to_variants1 = data_path + \
#     'snps/gatk_before_cortex/raw_variants_no_gvcf_mc10_ld/'
# name1 = 'gatk_no_gvcf_mc10_ld'
# path_to_variants1 = data_path + 'snps/gatk_minimap/raw_variants_ld/'
# name1 = 'gatk_minimap'
# path_to_variants1 = data_path + 'snps/cortex/raw_variants_ld/'
# name1 = 'cortex'
# path_to_variants1 = data_path + 'snps/skesa_minimap_mapq_raw_ld/'
# name1 = 'skesa_minimap'
# path_to_variants1 = data_path + 'snps/skesa_bwa_mapq_raw_ld/'
# name1 = 'skesa_bwa'
# path_to_variants1 = data_path + 'snps/pilon/raw_variants_filtered/'
path_to_variants1 = data_path + 'snps/pilon/annotated_pg_NWds10_filtered_test/'
name1 = 'pilon_filtered'
# path_to_variants1 = data_path + 'snps/annotated_with_DR_with_indel_with_pheno_and_snp_no_win_mqm_std3_mqm30_long_del_pg_filter_samples_first/'
# name1 = 'old_freebayes_ann'
# path_to_variants2 = data_path + 'snps/gatk_before_cortex/raw_variants_fixed_no_rep_gatk_ld/'
path_to_variants2 = data_path + 'snps/gatk_before_cortex/raw_variants_mq40_keep_complex_filtered/'
# path_to_variants2 = data_path + 'snps/gatk_before_cortex/annotated_pg_NWds10_mq40_keep_complex_filtered_test/'
# name2 = 'gatk_before_ld_fixed_no_rep_gatk'
name2 = 'gatk_before_mq40_keep_complex_filtered'
out_path = '../../res/testing/pilon_vs_gatk_annotated/'
# out_path = '../../res/testing/old_freebayes_ann_vs_gatk_test/'


path_to_short_tandem_repeats = data_path + 'h37rv.fasta.2.7.7.80.10.20.10.dat'
path_to_genes = data_path + 'all_walker_genes_all_drugs.list'
path_to_drug_to_genes = data_path + 'all_walker_genes.list'
gene_set_name = 'walker_genes'
path_to_pheno = data_path + 'pheno_mc5_mega/'


def find_intersections():
    parallel = Parallel(n_jobs=-1)

    with open(path_to_ids, 'r') as f:
        sample_ids = [name.strip() for name in f.readlines()]
    print(str(len(sample_ids)) + ' samples')

    variants1 = parallel(delayed(read_variants)(path_to_variants1, sample_id)
                         for sample_id in sample_ids)
    variants2 = parallel(delayed(read_variants)(path_to_variants2, sample_id)
                         for sample_id in sample_ids)
    all_variants = {}
    for sample_id, var_list in variants1:
        for var in var_list:
            counts = all_variants.get(var)
            if counts is None:
                all_variants[var] = [1, 0]
            else:
                counts[0] += 1
    for sample_id, var_list in variants2:
        for var in var_list:
            counts = all_variants.get(var)
            if counts is None:
                all_variants[var] = [0, 1]
            else:
                counts[1] += 1
    var_list = list(all_variants.keys())
    var_list.sort()
    with open(out_path + name1 + '_vs_' + name2 + '.csv', 'w') as f:
        f.write('pos\talt\ttype\t' + name1 + '\t' + name2 + '\n')
        for var in var_list:
            counts = all_variants[var]
            if counts[0] != counts[1]:
                f.write(var)
                f.write('\t%d\t%d\n' % (counts[0], counts[1]))
    with open(out_path + name1 + '_vs_' + name2 + '.' + name1 + '_unique', 'w') as f:
        for var in var_list:
            counts = all_variants[var]
            if counts[1] == 0:
                f.write(var)
                f.write('\t%d\n' % counts[0])
    with open(out_path + name1 + '_vs_' + name2 + '.' + name2 + '_unique', 'w') as f:
        for var in var_list:
            counts = all_variants[var]
            if counts[0] == 0:
                f.write(var)
                f.write('\t%d\n' % counts[1])


def find_intersections_annotated():
    parallel = Parallel(n_jobs=-1)

    with open(path_to_ids, 'r') as f:
        sample_ids = [name.strip() for name in f.readlines()]
    print(str(len(sample_ids)) + ' samples')

    variants1 = parallel(delayed(read_variants)(path_to_variants1, sample_id)
                         for sample_id in sample_ids)
    variants2 = parallel(delayed(read_variants)(path_to_variants2, sample_id)
                         for sample_id in sample_ids)
    all_variants = {}
    for sample_id, var_list in variants1:
        for var in var_list:
            counts = all_variants.get(var)
            if counts is None:
                all_variants[var] = [1, 0]
            else:
                counts[0] += 1
    for sample_id, var_list in variants2:
        for var in var_list:
            counts = all_variants.get(var)
            if counts is None:
                all_variants[var] = [0, 1]
            else:
                counts[1] += 1
    var_list = list(all_variants.keys())
    var_list.sort()
    with open(out_path + name1 + '_vs_' + name2 + '.csv', 'w') as f:
        f.write('variant\t' + name1 + '\t' + name2 + '\n')
        for var in var_list:
            counts = all_variants[var]
            if counts[0] != counts[1]:
                f.write(var)
                f.write('\t%d\t%d\n' % (counts[0], counts[1]))
    with open(out_path + name1 + '_vs_' + name2 + '.' + name1 + '_unique', 'w') as f:
        for var in var_list:
            counts = all_variants[var]
            if counts[1] == 0:
                f.write(var)
                f.write('\t%d\n' % counts[0])
    with open(out_path + name1 + '_vs_' + name2 + '.' + name2 + '_unique', 'w') as f:
        for var in var_list:
            counts = all_variants[var]
            if counts[0] == 0:
                f.write(var)
                f.write('\t%d\n' % counts[1])



def find_intersections_with_gene_set():
    parallel = Parallel(n_jobs=-1)

    with open(path_to_ids, 'r') as f:
        sample_ids = [name.strip() for name in f.readlines()]
    print(str(len(sample_ids)) + ' samples')

    variants1 = parallel(delayed(read_variants)(path_to_variants1, sample_id)
                         for sample_id in sample_ids)
    variants2 = parallel(delayed(read_variants)(path_to_variants2, sample_id)
                         for sample_id in sample_ids)
    all_variants_counts = {}
    for sample_id, var_list in variants1:
        for var in var_list:
            counts = all_variants_counts.get(var)
            if counts is None:
                all_variants_counts[var] = [1, 0]
            else:
                counts[0] += 1
    for sample_id, var_list in variants2:
        for var in var_list:
            counts = all_variants_counts.get(var)
            if counts is None:
                all_variants_counts[var] = [0, 1]
            else:
                counts[1] += 1
    var_list = list(all_variants_counts.keys())
    var_list.sort()

    pos_list = [int(var.split('\t')[0]) for var in var_list]
    cds_list = read_annotations(upstream_length)
    pos_to_cds = localize_all_variants(pos_list, cds_list)
    gene_set = set(l.strip() for l in open(path_to_genes).readlines())

    with open(out_path + name1 + '_vs_' + name2 + '_' + gene_set_name + '.csv', 'w') as f:
        f.write('pos\talt\ttype\tgene\t' + name1 + '\t' + name2 + '\n')
        for var in var_list:
            pos = int(var.split('\t')[0])
            cds = pos_to_cds.get(pos)
            if cds is None or cds.name not in gene_set:
                continue
            counts = all_variants_counts[var]
            f.write(var)
            f.write('\t%s\t%d\t%d\n' % (cds.name, counts[0], counts[1]))
    with open(out_path + name1 + '_vs_' + name2 + '_' + gene_set_name + '.' + name1 + '_unique', 'w') as f:
        for var in var_list:
            pos = int(var.split('\t')[0])
            cds = pos_to_cds.get(pos)
            if cds is None or cds.name not in gene_set:
                continue
            counts = all_variants_counts[var]
            if counts[1] == 0:
                f.write(var)
                f.write('\t%s\t%d\n' % (cds.name, counts[0]))
    with open(out_path + name1 + '_vs_' + name2 + '_' + gene_set_name + '.' + name2 + '_unique', 'w') as f:
        for var in var_list:
            pos = int(var.split('\t')[0])
            cds = pos_to_cds.get(pos)
            if cds is None or cds.name not in gene_set:
                continue
            counts = all_variants_counts[var]
            if counts[0] == 0:
                f.write(var)
                f.write('\t%s\t%d\n' % (cds.name, counts[1]))


def find_intersections_with_gene_set_annotated():
    parallel = Parallel(n_jobs=-1)

    with open(path_to_ids, 'r') as f:
        sample_ids = [name.strip() for name in f.readlines()]
    print(str(len(sample_ids)) + ' samples')

    variants1 = parallel(delayed(read_variants)(path_to_variants1, sample_id)
                         for sample_id in sample_ids)
    variants2 = parallel(delayed(read_variants)(path_to_variants2, sample_id)
                         for sample_id in sample_ids)
    all_variants_counts = {}
    for sample_id, var_list in variants1:
        for var in var_list:
            counts = all_variants_counts.get(var)
            if counts is None:
                all_variants_counts[var] = [1, 0]
            else:
                counts[0] += 1
    for sample_id, var_list in variants2:
        for var in var_list:
            counts = all_variants_counts.get(var)
            if counts is None:
                all_variants_counts[var] = [0, 1]
            else:
                counts[1] += 1
    var_list = list(all_variants_counts.keys())
    var_list.sort()

    gene_set = set(l.strip() for l in open(path_to_genes).readlines())

    with open(out_path + name1 + '_vs_' + name2 + '_' + gene_set_name + '_ann.csv', 'w') as f:
        f.write('name\tpos\tref\talt\ttype\t' + name1 + '\t' + name2 + '\n')
        for var in var_list:
            s = var.split('\t')
            if s[0] != 'Gene':
                continue
            if s[1] not in gene_set:
                continue
            counts = all_variants_counts[var]
            f.write('\t'.join(s[1:]))
            f.write('\t%d\t%d\n' % (counts[0], counts[1]))
    with open(out_path + name1 + '_vs_' + name2 + '_' + gene_set_name + '.' + name1 + '_unique_ann', 'w') as f:
        for var in var_list:
            s = var.split('\t')
            if s[0] != 'Gene':
                continue
            if s[1] not in gene_set:
                continue
            counts = all_variants_counts[var]
            if counts[1] == 0:
                f.write('\t'.join(s[1:]))
                f.write('\t%d\n' % counts[0])
    with open(out_path + name1 + '_vs_' + name2 + '_' + gene_set_name + '.' + name2 + '_unique_ann', 'w') as f:
        for var in var_list:
            s = var.split('\t')
            if s[0] != 'Gene':
                continue
            if s[1] not in gene_set:
                continue
            counts = all_variants_counts[var]
            if counts[0] == 0:
                f.write('\t'.join(s[1:]))
                f.write('\t%d\n' % counts[1])


def find_intersections_with_PGRS():
    parallel = Parallel(n_jobs=-1)

    with open(path_to_ids, 'r') as f:
        sample_ids = [name.strip() for name in f.readlines()]
    print(str(len(sample_ids)) + ' samples')

    variants1 = parallel(delayed(read_variants)(path_to_variants1, sample_id)
                         for sample_id in sample_ids)
    variants2 = parallel(delayed(read_variants)(path_to_variants2, sample_id)
                         for sample_id in sample_ids)
    all_variants_counts = {}
    for sample_id, var_list in variants1:
        for var in var_list:
            counts = all_variants_counts.get(var)
            if counts is None:
                all_variants_counts[var] = [1, 0]
            else:
                counts[0] += 1
    for sample_id, var_list in variants2:
        for var in var_list:
            counts = all_variants_counts.get(var)
            if counts is None:
                all_variants_counts[var] = [0, 1]
            else:
                counts[1] += 1
    var_list = list(all_variants_counts.keys())
    var_list.sort()

    pos_list = [int(var.split('\t')[0]) for var in var_list]
    cds_list = read_annotations(upstream_length)
    pos_to_cds = localize_all_variants(pos_list, cds_list)

    with open(out_path + name1 + '_vs_' + name2 + '_PPE_PGRS' + '.csv', 'w') as f:
        f.write('pos\talt\ttype\tgene\t' + name1 + '\t' + name2 + '\n')
        for var in var_list:
            pos = int(var.split('\t')[0])
            cds = pos_to_cds.get(pos)
            if cds is None or ('PPE' not in cds.name and 'PGRS' not in cds.name):
                continue
            counts = all_variants_counts[var]
            f.write(var)
            f.write('\t%s\t%d\t%d\n' % (cds.name, counts[0], counts[1]))
    with open(out_path + name1 + '_vs_' + name2 + '_PPE_PGRS' + '.' + name1 + '_unique', 'w') as f:
        for var in var_list:
            pos = int(var.split('\t')[0])
            cds = pos_to_cds.get(pos)
            if cds is None or ('PPE' not in cds.name and 'PGRS' not in cds.name):
                continue
            counts = all_variants_counts[var]
            if counts[1] == 0:
                f.write(var)
                f.write('\t%s\t%d\n' % (cds.name, counts[0]))
    with open(out_path + name1 + '_vs_' + name2 + '_PPE_PGRS' + '.' + name2 + '_unique', 'w') as f:
        for var in var_list:
            pos = int(var.split('\t')[0])
            cds = pos_to_cds.get(pos)
            if cds is None or ('PPE' not in cds.name and 'PGRS' not in cds.name):
                continue
            counts = all_variants_counts[var]
            if counts[0] == 0:
                f.write(var)
                f.write('\t%s\t%d\n' % (cds.name, counts[1]))


def find_intersections_with_gene_set_resistant():
    parallel = Parallel(n_jobs=-1)

    with open(path_to_ids, 'r') as f_res:
        sample_ids = [name.strip() for name in f_res.readlines()]
    print(str(len(sample_ids)) + ' samples')

    variants1 = parallel(delayed(read_variants)(path_to_variants1, sample_id)
                         for sample_id in sample_ids)
    variants2 = parallel(delayed(read_variants)(path_to_variants2, sample_id)
                         for sample_id in sample_ids)
    drug_to_gene_set = {}
    for l in open(path_to_drug_to_genes).readlines():
        s = l.strip().split('\t')
        genes = drug_to_gene_set.get(s[0])
        if genes is None:
            genes = set()
            drug_to_gene_set[s[0]] = genes
        genes.add(s[1])
    i = 0
    for fname in listdir(path_to_pheno):
        drug = fname[:-6]
        gene_set = drug_to_gene_set.get(drug)
        if gene_set is None:
            continue
        print('%d genes for %s' % (len(gene_set), drug))
        for gene in gene_set:
            print(gene)
        drug, sample_id_to_pheno = read_pheno(path_to_pheno, drug)
        sample_to_pheno = {sample_id: pheno for sample_id, pheno in sample_id_to_pheno}
        all_variants_counts_resistant = {}
        all_variants_counts_susceptible = {}
        for sample_id, var_list in variants1:
            pheno = sample_to_pheno.get(sample_id)
            if pheno is None:
                continue
            if pheno == 1:
                for var in var_list:
                    counts = all_variants_counts_resistant.get(var)
                    if counts is None:
                        all_variants_counts_resistant[var] = [1, 0]
                    else:
                        counts[0] += 1
            else:
                for var in var_list:
                    counts = all_variants_counts_susceptible.get(var)
                    if counts is None:
                        all_variants_counts_susceptible[var] = [1, 0]
                    else:
                        counts[0] += 1
        for sample_id, var_list in variants2:
            pheno = sample_to_pheno.get(sample_id)
            if pheno is None:
                continue
            if pheno == 1:
                for var in var_list:
                    counts = all_variants_counts_resistant.get(var)
                    if counts is None:
                        all_variants_counts_resistant[var] = [0, 1]
                    else:
                        counts[1] += 1
            else:
                for var in var_list:
                    counts = all_variants_counts_susceptible.get(var)
                    if counts is None:
                        all_variants_counts_susceptible[var] = [0, 1]
                    else:
                        counts[1] += 1
        all_vars_set = set()
        all_vars_set.update(all_variants_counts_resistant.keys())
        all_vars_set.update(all_variants_counts_susceptible.keys())
        var_list = list(all_vars_set)
        var_list.sort()

        pos_list = [int(var.split('\t')[0]) for var in var_list]
        cds_list = read_annotations(upstream_length)
        pos_to_cds = localize_all_variants(pos_list, cds_list)
        if i == 0:
            f_res = open(out_path + name1 + '_vs_' + name2 + '_' + gene_set_name + '_resistant.csv', 'w')
            f_res.write('pos\talt\ttype\tgene\t' + name1 + '\t' + name2 + '\n')
            f_sus = open(out_path + name1 + '_vs_' + name2 +
                         '_' + gene_set_name + '_susceptible.csv', 'w')
            f_sus.write('pos\talt\ttype\tgene\t' + name1 + '\t' + name2 + '\n')
        else:
            f_res = open(out_path + name1 + '_vs_' + name2 + '_' + gene_set_name + '_resistant.csv', 'a')
            f_sus = open(out_path + name1 + '_vs_' + name2 +
                         '_' + gene_set_name + '_susceptible.csv', 'a')
        f_res.write(drug + '\n')
        f_sus.write(drug + '\n')
        for var in var_list:
            pos = int(var.split('\t')[0])
            cds = pos_to_cds.get(pos)
            if cds is None or cds.name not in gene_set:
                continue
            counts = all_variants_counts_resistant.get(var)
            if counts is not None and counts[0] != counts[1]:
                f_res.write(var)
                f_res.write('\t%s\t%d\t%d\n' % (cds.name, counts[0], counts[1]))
            counts = all_variants_counts_susceptible.get(var)
            if counts is not None and counts[0] != counts[1]:
                f_sus.write(var)
                f_sus.write('\t%s\t%d\t%d\n' %
                            (cds.name, counts[0], counts[1]))
        f_res.close()
        f_sus.close()
        if i == 0:
            f_res = open(out_path + name1 + '_vs_' + name2 + '_' +
                         gene_set_name + '_resistant.' + name1 + '_unique', 'w')
            f_res.write('pos\talt\ttype\tgene\tcount\n')
            f_sus = open(out_path + name1 + '_vs_' + name2 + '_' +
                         gene_set_name + '_susceptible.' + name1 + '_unique', 'w')
            f_sus.write('pos\talt\ttype\tgene\tcount\n')
        else:
            f_res = open(out_path + name1 + '_vs_' + name2 + '_' + gene_set_name + '_resistant.' + name1 + '_unique', 'a')
            f_sus = open(out_path + name1 + '_vs_' + name2 + '_' + gene_set_name + '_susceptible.' + name1 + '_unique', 'a')
        f_res.write(drug + '\n')
        f_sus.write(drug + '\n')
        for var in var_list:
            pos = int(var.split('\t')[0])
            cds = pos_to_cds.get(pos)
            if cds is None or cds.name not in gene_set:
                continue
            counts = all_variants_counts_resistant.get(var)
            if counts is not None and counts[1] == 0:
                f_res.write(var)
                f_res.write('\t%s\t%d\n' % (cds.name, counts[0]))
            counts = all_variants_counts_susceptible.get(var)
            if counts is not None and counts[1] == 0:
                f_sus.write(var)
                f_sus.write('\t%s\t%d\n' % (cds.name, counts[0]))
        f_res.close()
        f_sus.close()
        if i == 0:
            f_res = open(out_path + name1 + '_vs_' + name2 + '_' + gene_set_name + '_resistant.' + name2 + '_unique', 'w')
            f_res.write('pos\talt\ttype\tgene\tcount\n')
            f_sus = open(out_path + name1 + '_vs_' + name2 + '_' +
                         gene_set_name + '_susceptible.' + name2 + '_unique', 'w')
            f_sus.write('pos\talt\ttype\tgene\tcount\n')
        else:
            f_res = open(out_path + name1 + '_vs_' + name2 + '_' + gene_set_name + '_resistant.' + name2 + '_unique', 'a')
            f_sus = open(out_path + name1 + '_vs_' + name2 + '_' +
                         gene_set_name + '_susceptible.' + name2 + '_unique', 'a')
        f_res.write(drug + '\n')
        f_sus.write(drug + '\n')
        for var in var_list:
            pos = int(var.split('\t')[0])
            cds = pos_to_cds.get(pos)
            if cds is None or cds.name not in gene_set:
                continue
            counts = all_variants_counts_resistant.get(var)
            if counts is not None and counts[0] == 0:
                f_res.write(var)
                f_res.write('\t%s\t%d\n' % (cds.name, counts[1]))
            counts = all_variants_counts_susceptible.get(var)
            if counts is not None and counts[0] == 0:
                f_sus.write(var)
                f_sus.write('\t%s\t%d\n' % (cds.name, counts[1]))
        f_res.close()
        f_sus.close()
        i += 1


def find_intersections_with_gene_set_resistant_annotated():
    parallel = Parallel(n_jobs=-1)

    with open(path_to_ids, 'r') as f_res:
        sample_ids = [name.strip() for name in f_res.readlines()]
    print(str(len(sample_ids)) + ' samples')

    variants1 = parallel(delayed(read_variants)(path_to_variants1, sample_id)
                         for sample_id in sample_ids)
    variants2 = parallel(delayed(read_variants)(path_to_variants2, sample_id)
                         for sample_id in sample_ids)
    drug_to_gene_set = {}
    for l in open(path_to_drug_to_genes).readlines():
        s = l.strip().split('\t')
        genes = drug_to_gene_set.get(s[0])
        if genes is None:
            genes = set()
            drug_to_gene_set[s[0]] = genes
        genes.add(s[1])
    i = 0
    for fname in listdir(path_to_pheno):
        drug = fname[:-6]
        gene_set = drug_to_gene_set.get(drug)
        if gene_set is None:
            continue
        print('%d genes for %s' % (len(gene_set), drug))
        for gene in gene_set:
            print(gene)
        drug, sample_id_to_pheno = read_pheno(path_to_pheno, drug)
        sample_to_pheno = {sample_id: pheno for sample_id, pheno in sample_id_to_pheno}
        all_variants_counts_resistant = {}
        all_variants_counts_susceptible = {}
        for sample_id, var_list in variants1:
            pheno = sample_to_pheno.get(sample_id)
            if pheno is None:
                continue
            if pheno == 1:
                for var in var_list:
                    counts = all_variants_counts_resistant.get(var)
                    if counts is None:
                        all_variants_counts_resistant[var] = [1, 0]
                    else:
                        counts[0] += 1
            else:
                for var in var_list:
                    counts = all_variants_counts_susceptible.get(var)
                    if counts is None:
                        all_variants_counts_susceptible[var] = [1, 0]
                    else:
                        counts[0] += 1
        for sample_id, var_list in variants2:
            pheno = sample_to_pheno.get(sample_id)
            if pheno is None:
                continue
            if pheno == 1:
                for var in var_list:
                    counts = all_variants_counts_resistant.get(var)
                    if counts is None:
                        all_variants_counts_resistant[var] = [0, 1]
                    else:
                        counts[1] += 1
            else:
                for var in var_list:
                    counts = all_variants_counts_susceptible.get(var)
                    if counts is None:
                        all_variants_counts_susceptible[var] = [0, 1]
                    else:
                        counts[1] += 1
        all_vars_set = set()
        all_vars_set.update(all_variants_counts_resistant.keys())
        all_vars_set.update(all_variants_counts_susceptible.keys())
        var_list = list(all_vars_set)
        var_list.sort()

        if i == 0:
            f_res = open(out_path + name1 + '_vs_' + name2 + '_' + gene_set_name + '_resistant_ann.csv', 'w')
            f_res.write('gene\tpos\tref\talt\ttype\t' + name1 + '\t' + name2 + '\n')
            f_sus = open(out_path + name1 + '_vs_' + name2 +
                         '_' + gene_set_name + '_susceptible_ann.csv', 'w')
            f_sus.write('gene\tpos\tref\talt\ttype\t' + name1 + '\t' + name2 + '\n')
        else:
            f_res = open(out_path + name1 + '_vs_' + name2 + '_' + gene_set_name + '_resistant_ann.csv', 'a')
            f_sus = open(out_path + name1 + '_vs_' + name2 +
                         '_' + gene_set_name + '_susceptible_ann.csv', 'a')
        f_res.write(drug + '\n')
        f_sus.write(drug + '\n')
        for var in var_list:
            s = var.split('\t')
            if s[0] != 'Gene':
                continue
            if s[1] not in gene_set:
                continue
            counts = all_variants_counts_resistant.get(var)
            if counts is not None and counts[0] != counts[1]:
                f_res.write('\t'.join(s[1:]))
                f_res.write('\t%d\t%d\n' % (counts[0], counts[1]))
            counts = all_variants_counts_susceptible.get(var)
            if counts is not None and counts[0] != counts[1]:
                f_sus.write('\t'.join(s[1:]))
                f_sus.write('\t%d\t%d\n' %
                            (counts[0], counts[1]))
        f_res.close()
        f_sus.close()
        if i == 0:
            f_res = open(out_path + name1 + '_vs_' + name2 + '_' +
                         gene_set_name + '_resistant_ann.' + name1 + '_unique', 'w')
            f_res.write('gene\tpos\tref\talt\ttype\tcount\n')
            f_sus = open(out_path + name1 + '_vs_' + name2 + '_' +
                         gene_set_name + '_susceptible_ann.' + name1 + '_unique', 'w')
            f_sus.write('gene\tpos\tref\talt\ttype\tcount\n')
        else:
            f_res = open(out_path + name1 + '_vs_' + name2 + '_' + gene_set_name + '_resistant_ann.' + name1 + '_unique', 'a')
            f_sus = open(out_path + name1 + '_vs_' + name2 + '_' + gene_set_name + '_susceptible_ann.' + name1 + '_unique', 'a')
        f_res.write(drug + '\n')
        f_sus.write(drug + '\n')
        for var in var_list:
            s = var.split('\t')
            if s[0] != 'Gene':
                continue
            if s[1] not in gene_set:
                continue
            counts = all_variants_counts_resistant.get(var)
            if counts is not None and counts[1] == 0:
                f_res.write('\t'.join(s[1:]))
                f_res.write('\t%d\n' % counts[0])
            counts = all_variants_counts_susceptible.get(var)
            if counts is not None and counts[1] == 0:
                f_sus.write('\t'.join(s[1:]))
                f_sus.write('\t%d\n' % (counts[0]))
        f_res.close()
        f_sus.close()
        if i == 0:
            f_res = open(out_path + name1 + '_vs_' + name2 + '_' + gene_set_name + '_resistant_ann.' + name2 + '_unique', 'w')
            f_res.write('gene\tpos\tref\talt\ttype\tcount\n')
            f_sus = open(out_path + name1 + '_vs_' + name2 + '_' +
                         gene_set_name + '_susceptible_ann.' + name2 + '_unique', 'w')
            f_sus.write('gene\tpos\tref\talt\ttype\tcount\n')
        else:
            f_res = open(out_path + name1 + '_vs_' + name2 + '_' + gene_set_name + '_resistant_ann.' + name2 + '_unique', 'a')
            f_sus = open(out_path + name1 + '_vs_' + name2 + '_' +
                         gene_set_name + '_susceptible_ann.' + name2 + '_unique', 'a')
        f_res.write(drug + '\n')
        f_sus.write(drug + '\n')
        for var in var_list:
            s = var.split('\t')
            if s[0] != 'Gene':
                continue
            if s[1] not in gene_set:
                continue
            counts = all_variants_counts_resistant.get(var)
            if counts is not None and counts[0] == 0:
                f_res.write('\t'.join(s[1:]))
                f_res.write('\t%d\n' % counts[1])
            counts = all_variants_counts_susceptible.get(var)
            if counts is not None and counts[0] == 0:
                f_sus.write('\t'.join(s[1:]))
                f_sus.write('\t%d\n' % counts[1])
        f_res.close()
        f_sus.close()
        i += 1


def compute_enrichment(path, coords):
    starts = [c[0] for c in coords]
    total = 0
    in_intervals = 0
    for l in open(path).readlines():
        pos = int(l.strip().split('\t')[0].split('_')[0])
        i = bisect_left(starts, pos)
        if i == 0:
            if starts[i] == pos:
                in_intervals += 1
        elif i < len(starts):
            if starts[i] == pos:
                in_intervals += 1
            else:
                s, e = coords[i - 1]
                if pos < e:
                    in_intervals += 1
        total += 1
    if total > 0:
        return in_intervals/total
    else:
        return 0


def compute_interval_enrichment():
    coords = []
    l = 0
    with open(path_to_short_tandem_repeats) as f:
        for line in f.readlines()[15:]:
            s = line.split()
            st = int(s[0])
            end = int(s[1]) + 1
            coords.append((st, end))
            l += end - st
    coords.sort(key=lambda tup: tup[0])
    out_path_unique1 = out_path + name1 + '_vs_' + name2 + '.' + name1 + '_unique'
    enrich1 = compute_enrichment(out_path_unique1, coords)
    out_path_unique2 = out_path + name1 + '_vs_' + name2 + '.' + name2 + '_unique'
    enrich2 = compute_enrichment(out_path_unique2, coords)
    print('repeat from ref: %1.2f %s unique repeat rate: %1.2f %s unique repeat rate: %1.2f' % 
    (l/ref_len, name1, enrich1, name2, enrich2))


def indel_distribution(path_to_variants, dist_path):
    parallel = Parallel(n_jobs=-1)

    with open(path_to_ids, 'r') as f:
        sample_ids = [name.strip() for name in f.readlines()]
    print(str(len(sample_ids)) + ' samples')

    variants = parallel(delayed(read_variants)(path_to_variants, sample_id)
                         for sample_id in sample_ids)
    all_vars = set()
    for sample_id, var_list in variants:
        all_vars.update(var_list)
            
    ins_len_to_count = {}
    del_len_to_count = {}
    snp_num = 0
    ins_num = 0
    del_num = 0
    total = len(all_vars)
    for v in all_vars:
        s = v.split('\t')
        if s[-1] == 'snp':
            snp_num += 1
        elif s[-1] == 'ins':
            ins_num += 1
            l = len(s[1])
            c = ins_len_to_count.get(l)
            if c is None:
                ins_len_to_count[l] = 1
            else:
                ins_len_to_count[l] = c + 1
        else:
            del_num += 1
            l = int(s[1])
            c = del_len_to_count.get(l)
            if c is None:
                del_len_to_count[l] = 1
            else:
                del_len_to_count[l] = c + 1
    with open(dist_path, 'w') as f:
        f.write('snp_num=%d(%1.2f) ins_num=%d(%1.2f) del_num=%d(%1.2f)\n' % 
        (snp_num, snp_num/total, ins_num, ins_num/total, del_num, del_num/total))
        f.write('inserts:\n')
        ins_len_to_count_list = [(l, c) for l, c in ins_len_to_count.items()]
        ins_len_to_count_list.sort(key=lambda x: x[0])
        for l, c in ins_len_to_count_list:
            f.write('%d %d (%1.2f)\n' % (l, c, c/ins_num))
        f.write('dels:\n')
        del_len_to_count_list = [(l, c) for l, c in del_len_to_count.items()]
        del_len_to_count_list.sort(key=lambda x: x[0])
        for l, c in del_len_to_count_list:
            f.write('%d %d (%1.2f)\n' % (l, c, c/ins_num))


def indel_distribution_intervals(path_to_variants, dist_path):
    parallel = Parallel(n_jobs=-1)

    with open(path_to_ids, 'r') as f:
        sample_ids = [name.strip() for name in f.readlines()]
    print(str(len(sample_ids)) + ' samples')

    variants = parallel(delayed(read_variants)(path_to_variants, sample_id)
                        for sample_id in sample_ids)
    all_vars = set()
    for sample_id, var_list in variants:
        all_vars.update(var_list)

    coords = []
    with open(path_to_short_tandem_repeats) as f:
        for line in f.readlines()[15:]:
            s = line.split()
            st = int(s[0])
            end = int(s[1]) + 1
            coords.append((st, end))
    coords.sort(key=lambda tup: tup[0])
    starts = [c[0] for c in coords]

    ins_len_to_count = {}
    del_len_to_count = {}
    snp_num = 0
    ins_num = 0
    del_num = 0
    total = len(all_vars)
    for v in all_vars:
        s = v.split('\t')
        pos = int(s[0])
        in_interval = False
        i = bisect_left(starts, pos)
        if i == 0:
            if starts[i] == pos:
                in_interval = True
        elif i < len(starts):
            if starts[i] == pos:
                in_interval = True
            else:
                if pos < coords[i - 1][1]:
                    in_interval = True
        if not in_interval:
            continue
        if s[-1] == 'snp':
            snp_num += 1
        elif s[-1] == 'ins':
            ins_num += 1
            l = len(s[1])
            c = ins_len_to_count.get(l)
            if c is None:
                ins_len_to_count[l] = 1
            else:
                ins_len_to_count[l] = c + 1
        else:
            del_num += 1
            l = int(s[1])
            c = del_len_to_count.get(l)
            if c is None:
                del_len_to_count[l] = 1
            else:
                del_len_to_count[l] = c + 1
    with open(dist_path, 'w') as f:
        f.write('snp_num=%d(%1.2f) ins_num=%d(%1.2f) del_num=%d(%1.2f)\n' %
                (snp_num, snp_num / total, ins_num, ins_num / total, del_num, del_num / total))
        f.write('inserts:\n')
        ins_len_to_count_list = [(l, c) for l, c in ins_len_to_count.items()]
        ins_len_to_count_list.sort(key=lambda x: x[0])
        for l, c in ins_len_to_count_list:
            f.write('%d %d (%1.2f)\n' % (l, c, c / ins_num))
        f.write('dels:\n')
        del_len_to_count_list = [(l, c) for l, c in del_len_to_count.items()]
        del_len_to_count_list.sort(key=lambda x: x[0])
        for l, c in del_len_to_count_list:
            f.write('%d %d (%1.2f)\n' % (l, c, c / ins_num))


def sample_wize(sample_stat_path):
    parallel = Parallel(n_jobs=-1)

    with open(path_to_ids, 'r') as f:
        sample_ids = [name.strip() for name in f.readlines()]
    print(str(len(sample_ids)) + ' samples')

    variants1 = parallel(delayed(read_variants)(path_to_variants1, sample_id)
                         for sample_id in sample_ids)
    variants2 = parallel(delayed(read_variants)(path_to_variants2, sample_id)
                         for sample_id in sample_ids)
    sample_to_vars1 = {sample_id: var_list for sample_id, var_list in variants1}
    with open(sample_stat_path, 'w') as f:
        f.write('sample\t%s\t%s\n' % (name1, name2))
        for sample, var_list2 in variants2:
            var_set1 = set(sample_to_vars1[sample])
            var_set2 = set(var_list2)
            int_len = len(var_set1.intersection(var_set2))
            f.write('%s\t%d\t%d\n' % (sample, len(var_set1) - int_len, len(var_set2) - int_len))


if __name__ == "__main__":
    if not exists(out_path):
        makedirs(out_path)
    # find_intersections()
    find_intersections_annotated()
    # compute_interval_enrichment()
    # find_intersections_with_gene_set()
    # find_intersections_with_PGRS()
    # find_intersections_with_gene_set_resistant()
    # find_intersections_with_gene_set_annotated()
    # find_intersections_with_gene_set_resistant_annotated()
    # indel_distribution(path_to_variants1, out_path + name1 + '.indel_stat')
    # indel_distribution_intervals(path_to_variants1, out_path + name1 + '_repeats.indel_stat')
    # sample_wize(out_path + name1 + '_vs_' + name2 + '.by_sample')
