from os import listdir

from sklearn.externals.joblib import Parallel, delayed
import matplotlib.pyplot as plt
from core.constants import data_path
from core.data_reading import read_variants, read_all_variants

# path_to_var = data_path + 'snps/gatk_before_cortex/raw_variants_std2/'
# path_to_var = data_path + 'snps/gatk_before_cortex/raw_variants_std3/'
# path_to_var = data_path + 'snps/annotated_fixed_no_rep_long_del_pg_NWds10_combined_new/'
from data_processing.print_alignment_gatk import get_intervals_to_filter_out, read_snps

path_to_var = data_path + 'snps/combined_raw_variants_mq40_keep_complex_std_names_filtered/'
# path_to_var = data_path + 'snps/gatk_before_cortex/annotated_pg_NWds10_mq40_keep_complex_filtered/'
# path_to_var = data_path + 'snps/annotated_with_DR_with_indel_with_pheno_and_snp_no_win_qual_mqm_std3_mqm30_no_highcov_long_del_pg_NWds10/'
# path_to_ids = data_path + 'all_with_pheno.txt'
path_to_ids = data_path + 'snps/combined_raw_variants_mq40_keep_complex_std_names_filtered/samples_filtered.list'
# out_path = '../../res/ml_log_mc3_gatk_before_std2/var.counts'
# out_path = '../../res/ml_log_mc3_gatk_before_std3/var.counts'
# out_path = '../../res/base_dataset_var.counts'
# out_path = '../../res/base_dataset_new_var.counts'
# out_path = '../../res/full_dataset_new_var.counts'
out_path = '../../res/full_dataset_raw_5plus.counts'


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


filter_by_count = True
min_count = 5


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
    counts.sort(key=lambda x: x[0].split()[1])
    # counts.sort(key=lambda x: int(x[0].split()[0]))
    with open(out_path, 'w') as f:
        for snp, c in counts:
            if filter_by_count and c < min_count:
                continue
            f.write('%s\t%d\n' % (snp, c))


path_to_mut_counts = '../../res/mut_in_samples001.counts'
print_hist = False
path_to_hist = '../../res/mut_in_samples_hist.png'
freq_threshold = 0.01


def count_mut_in_samples():
    sample_ids = [name.strip() for name in open(path_to_ids).readlines()]
    print(str(len(sample_ids)) + ' samples')

    tasks = Parallel(n_jobs=-1)(delayed(read_variants)(path_to_var, sample_id)
                     for sample_id in sample_ids)
    sample_to_variants = {}
    mut_counts = {}
    for name, l in tasks:
        sample_to_variants[name] = l
        for snp in l:
            c = mut_counts.get(snp)
            if c is None:
                mut_counts[snp] = 1
            else:
                mut_counts[snp] = c + 1
    rare_mut_set = set()
    c_max = freq_threshold * len(sample_ids)
    for snp, c in mut_counts.items():
        if c < c_max:
            rare_mut_set.add(snp)

    filter_intervals = get_intervals_to_filter_out()
    tasks = Parallel(n_jobs=-1)(delayed(read_snps)(sample_id, filter_intervals) for sample_id in sample_ids)
    # sample_to_snps = {}
    # for sample_id, snps in tasks:
    #     var_list = [str(pos) + '\t' + snp for pos, snp in snps.items()]
    #     sample_to_snps[sample_id] = var_list

    total_mut_counts = []
    with open(path_to_mut_counts, 'w') as f:
        f.write('sample\ttotal_variants\trare_variants\tfiltered_snps\n')
        for name, mlist in tasks:
            all_var_list = sample_to_variants[name]
            rnum = 0
            for m in all_var_list:
                if m in rare_mut_set:
                    rnum += 1
            total_mut_counts.append(len(all_var_list))
            f.write('%s\t%d\t%d\t%d\n' % (name, len(all_var_list), rnum, len(mlist)))
    if print_hist:
        fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
        axs.hist(total_mut_counts, bins=20)
        plt.savefig(path_to_hist)


path_to_counts = '../../res/new_dataset_stats/full_dataset_new_var.counts'
dict_path = '../../data/dictionaries_indels/Walker_dictionary_broken_genes.txt'
filtered_counts = '../../res/new_dataset_stats/full_dataset_Walker.counts'


def filter_counts_by_dict():
    walker_muts = {l[:l.rfind('\t')]: '\t0\n' for l in open(dict_path).readlines() if 'snp' in l}
    with open(filtered_counts, 'w') as f:
        for l in open(path_to_counts).readlines():
            mut = l[:l.rfind('\t')]
            if mut in walker_muts:
                walker_muts[mut] = l[l.rfind('\t'):]
        for mut, c in walker_muts.items():
            f.write(mut + c)


path_to_annotations = '../../data/ref_seqs/'
pe_counts_path = '../../data/ref_seqs/pe.counts'
pe_genes_path = '../../data/pe_genes/'


def count_ppe_genes_in_anotatations():
    with open(pe_counts_path, 'w') as f:
        for fname in listdir(path_to_annotations):
            if fname.endswith('.gff3'):
                c = 0
                with open(pe_genes_path + fname, 'w') as f1:
                    for l in open(path_to_annotations + fname).readlines():
                        if 'CDS' not in l:
                            continue
                        s = l.split(';')
                        for si in s:
                            if 'gene' in si:
                                if 'PE' in si or 'PGRS' in si:
                                    c += 1
                                    f1.write(l)
                                    break
                            if 'product' in si:
                                if 'PE FAMILY PROTEIN' in si.upper() or 'PGRS FAMILY PROTEIN' in si.upper() or\
                                        'PE-' in si.upper():
                                    c += 1
                                    f1.write(l)
                                break
                    f.write('%s\t%d\n' % (fname, c))


def find_duplicates():
    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    sample_to_variants = read_all_variants(path_to_var, sample_ids)
    var_sets = [(sample, set(var_list)) for sample, var_list in sample_to_variants.items()]
    num_to_var_sets = {}
    for sample, var_set in var_sets:
        ni = len(var_set)
        set_list = num_to_var_sets.get(ni)
        if set_list is None:
            set_list = [(sample, var_set)]
            num_to_var_sets[ni] = set_list
        else:
            set_list.append((sample, var_set))
    for ni, set_list in num_to_var_sets.items():
        for i in range(len(set_list)):
            sample_i, vars_i = set_list[i]
            for j in range(i + 1, len(set_list)):
                sample_j, vars_j = set_list[j]
                if len(vars_i.intersection(vars_j)) == ni:
                    print(sample_i + ' ' + sample_j)


def find_duplicates_snp_noDR():
    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    filter_intervals = get_intervals_to_filter_out()
    tasks = Parallel(n_jobs=-1)(delayed(read_snps)(sample_id, filter_intervals) for sample_id in sample_ids)
    sample_to_snps = {}
    for sample_id, snps in tasks:
        var_list = [str(pos) + '\t' + snp for pos, snp in snps.items()]
        sample_to_snps[sample_id] = var_list
    var_sets = [(sample, set(var_list)) for sample, var_list in sample_to_snps.items()]
    num_to_var_sets = {}
    for sample, var_set in var_sets:
        ni = len(var_set)
        set_list = num_to_var_sets.get(ni)
        if set_list is None:
            set_list = [(sample, var_set)]
            num_to_var_sets[ni] = set_list
        else:
            set_list.append((sample, var_set))
    for ni, set_list in num_to_var_sets.items():
        for i in range(len(set_list)):
            sample_i, vars_i = set_list[i]
            for j in range(i + 1, len(set_list)):
                sample_j, vars_j = set_list[j]
                if len(vars_i.intersection(vars_j)) == ni:
                    print(sample_i + ' ' + sample_j)


def dist_matrix():
    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    filter_intervals = get_intervals_to_filter_out()
    tasks = Parallel(n_jobs=-1)(delayed(read_snps)(sample_id, filter_intervals) for sample_id in sample_ids)
    sample_to_snps = {}
    for sample_id, snps in tasks:
        var_list = [str(pos) + '\t' + snp for pos, snp in snps.items()]
        sample_to_snps[sample_id] = var_list
    snp_sets = [(sample, set(var_list)) for sample, var_list in sample_to_snps.items()]

    sample_to_variants = read_all_variants(path_to_var, sample_ids)
    var_sets = [(sample, set(var_list)) for sample, var_list in sample_to_variants.items()]
    for i in range(len(var_sets)):
        sample_i, vars_i = var_sets[i]
        sample_i, snps_i = snp_sets[i]
        for j in range(i + 1, len(var_sets)):
            sample_j, vars_j = var_sets[j]
            sample_j, snps_j = snp_sets[j]
            var_sim = len(vars_i.intersection(vars_j))/len(vars_i.union(vars_j))
            snp_sim = len(snps_i.intersection(snps_j))/len(snps_i.union(snps_j))
            print('%s %s %1.2f %1.2f' % (sample_i, sample_j, var_sim, snp_sim))


if __name__ == '__main__':
    # count_unique()
    # print_all_counts()
    count_mut_in_samples()
    # filter_counts_by_dict()
    # count_ppe_genes_in_anotatations()
    # find_duplicates()
    # find_duplicates_snp_noDR()
    # dist_matrix()