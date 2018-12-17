from math import sqrt
from os import makedirs

from os.path import exists

from Bio.SubsMat import MatrixInfo
from scipy.sparse import lil_matrix
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals.joblib import Parallel, delayed
import numpy as np
import os

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score

from src.core.constants import aminoacids, tp, tn, fp, fn
from src.core.data_reading import read_pheno, read_subset, read_variants

path_to_pheno = './data/pheno/'
path_to_snp = './data/snps/annotated_with_DR_with_indel_with_pheno_and_snp_mc10/'
path_to_snp_list = path_to_snp + 'all_var_list.csv'
path_to_ids = './data/dr_covered_with_pheno_and_snp.txt'
path_to_subsets = './data/subsets/'

out_path = './res/fs_models_multi.out'

transformed_out_path = './data/snps/annotated_transformed_genes_with_pheno_and_snp_mc10/'
path_to_transformed_list = transformed_out_path + 'all_var_list.csv'

only_print_features = False
merge_all_mut_in_gene = True
filter_non_cds = True
merge_long_dels = True
merge_all_mut_in_pos = True
snp_count_threshold = 3
use_pairs = False
snp_count_r_threshold = 3
jaccard_index_threshold = 0.9
thread_num = 144


def create_distance_matrix(similarity_matrix):
    dist_matrix = {}
    aa_to_vec = {}
    for aa1 in aminoacids:
        vec = []
        for aa2 in aminoacids:
            if (aa1, aa2) in similarity_matrix.keys():
                vec.append(similarity_matrix[(aa1, aa2)])
            else:
                vec.append(similarity_matrix[(aa2, aa1)])
        aa_to_vec[aa1] = vec
    for aa1 in aminoacids:
        vec1 = aa_to_vec[aa1]
        for aa2 in aminoacids:
            vec2 = aa_to_vec[aa2]
            dist = 0
            for i in range(len(aminoacids)):
                dist += (vec1[i] - vec2[i])**2
            dist_matrix[(aa1, aa2)] = sqrt(dist)
    return dist_matrix


def split_dataset(drug, pheno, sample_to_snps, train_subset):
    snp_train = {}
    snp_test = {}
    pheno_train = {}
    pheno_test = {}
    for sample_id, val in pheno:
        if sample_id in sample_to_snps.keys():
            if sample_id in train_subset:
                snp_train[sample_id] = sample_to_snps[sample_id]
                pheno_train[sample_id] = val
            else:
                snp_test[sample_id] = sample_to_snps[sample_id]
                pheno_test[sample_id] = val
    return drug, snp_train, snp_test, pheno_train, pheno_test


def stat_fixed_multi(X_train, y_train, X_test, y_test, clf):
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    tp_score = tp(y_test, y_pred)
    tn_score = tn(y_test, y_pred)
    fp_score = fp(y_test, y_pred)
    fn_score = fn(y_test, y_pred)
    if tp_score + fp_score > 0:
        ppv = tp_score/(tp_score + fp_score)
    else:
        ppv = 0
    if tn_score + fn_score > 0:
        npv = tn_score/(tn_score + fn_score)
    else:
        npv = 0
    if tp_score + fn_score > 0:
        sensitivity = tp_score/(tp_score + fn_score)
    else:
        sensitivity = 0
    if tn_score + fp_score > 0:
        specificity = tn_score/(tn_score + fp_score)
    else:
        specificity = 0
    return ppv, npv, sensitivity, specificity, f1_score(y_test, y_pred)


def snps_counts(sample_to_snps):
    snp_to_count = {}
    for l in sample_to_snps.values():
        for snp in l:
            s = snp.split('\t')
            if (s[0] == 'Gene' and s[-1] == 'snp') or s[-1] == 'dist':
                f = '\t'.join(s[:-2])
                c = snp_to_count.get(f)
                if c is not None:
                    snp_to_count[f] = c + 1
                else:
                    snp_to_count[f] = 1
            else:
                c = snp_to_count.get(snp)
                if c is not None:
                    snp_to_count[snp] = c + 1
                else:
                    snp_to_count[snp] = 1
    return snp_to_count


def filter_snps(drug, snp_train, snp_test):
    all_snps_set = set()
    sample_num = len(snp_train)
    snp_to_count = snps_counts(snp_train)
    snp_train_filtered = {}
    snp_test_filtered = {}
    for sample_id, snp_list in snp_train.items():
        res = []
        for snp in snp_list:
            s = snp.split('\t')
            if (s[0] == 'Gene' and s[-1] == 'snp') or s[-1] == 'dist':
                f = '\t'.join(s[:-2])
                if snp_count_threshold < snp_to_count[f] < sample_num - snp_count_threshold:
                    res.append(snp)
                    all_snps_set.add(f)
            else:
                if snp_count_threshold < snp_to_count[snp] < sample_num - snp_count_threshold:
                    res.append(snp)
                    all_snps_set.add(snp)
        snp_train_filtered[sample_id] = res
    for sample_id, snp_list in snp_test.items():
        snp_test_filtered[sample_id] = [snp for snp in snp_list if snp in all_snps_set]
    full_snp_list = list(all_snps_set)
    return drug, snp_train_filtered, snp_test_filtered, full_snp_list


def drug_stat_multi(drug, snp_train, snp_test, pheno_train, pheno_test, clf, full_snp_list):
    if len(pheno_train) == 0:
        return drug, len(pheno_train) + len(pheno_test), 0, 0, 0, 0, 0

    sample_ids = list(snp_train.keys())
    sample_num = len(sample_ids)
    snp_num = len(full_snp_list)

    feature_to_index = {}
    for i in range(snp_num):
        feature_to_index[full_snp_list[i]] = i

    X_train = lil_matrix((sample_num, snp_num), dtype=int)
    y_train = []
    for i in range(len(sample_ids)):
        for f in snp_train[sample_ids[i]]:
            s = f.split('\t')
            if s[-1] == 'dist':
                f = '\t'.join(s[:-2])
                X_train[i, feature_to_index[f]] = float(s[-2])
            elif s[0] == 'Gene' and s[-1] == 'snp':
                f = '\t'.join(s[:-2])
                X_train[i, feature_to_index[f]] = float(s[-2])
            else:
                X_train[i, feature_to_index[f]] = 1
        y_train.append(pheno_train[sample_ids[i]])

    sample_ids = list(snp_test.keys())
    sample_num = len(sample_ids)
    X_test = lil_matrix((sample_num, snp_num), dtype=int)
    y_test = []
    for i in range(len(sample_ids)):
        for f in snp_test[sample_ids[i]]:
            s = f.split('\t')
            if s[-1] == 'dist':
                f = '\t'.join(s[:-2])
                X_train[i, feature_to_index[f]] = float(s[-2])
            elif s[0] == 'Gene' and s[-1] == 'snp':
                f = '\t'.join(s[:-2])
                X_train[i, feature_to_index[f]] = float(s[-2])
            else:
                X_test[i, feature_to_index[f]] = 1
        y_test.append(pheno_test[sample_ids[i]])

    ppv, npv, sensitivity, specificity, f1 = stat_fixed_multi(X_train, y_train, X_test, y_test, clf)
    return drug, len(pheno_train) + len(pheno_test), ppv, npv, sensitivity, specificity, f1


def print_stat_fixed_multi(f, drug_to_splits, clf, n_jobs=-1):
    f.write('Drug\tSamples\tPPV\tNPV\tSensitivity\tSpecificity\tF1\n')
    if n_jobs != 1:
        tasks = Parallel(n_jobs=n_jobs)(delayed(drug_stat_multi)(drug, snp_train, snp_test, pheno_train, pheno_test, clf, full_snp_list)
                                        for drug, (snp_train, snp_test, pheno_train, pheno_test, full_snp_list)
                                        in drug_to_splits.items())
        for drug, samples_with_pheno, ppv, npv, sensitivity, specificity, f1 in tasks:
            f.write(drug)
            f.write('\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n' % (samples_with_pheno, ppv, npv, sensitivity, specificity, f1))
    else:
        for drug, (snp_train, snp_test, pheno_train, pheno_test, full_snp_list) in drug_to_splits.items():
            drug, samples_with_pheno, ppv, npv, sensitivity, specificity, f1 = \
                drug_stat_multi(drug, snp_train, snp_test, pheno_train, pheno_test, clf, full_snp_list)
            f.write(drug)
            f.write('\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n' % (samples_with_pheno, ppv, npv, sensitivity, specificity, f1))


def process_snps(sample_id, var_list, distance_matrix, indel_cost):
    res = set()
    gene_to_dist = {}
    deletions = {}
    broken_genes = set()
    for var in var_list:
        s = var.split('\t')
        if filter_non_cds and s[0] == 'non_cds':
            continue
        if s[0] == 'Gene' and (s[4] == '*' or s[3] == '*' or s[5] == 'FS'):
            broken_genes.add(s[0] + '\t' + s[1])
        else:
            if s[0] == 'Gene':
                if s[5] == 'snp':
                    dist = distance_matrix[(s[3], s[4])]
                    res.add(s[0] + '\t' + s[1] + '\t' + s[2] + '\t' + str(dist) + '\tsnp')
                    if merge_all_mut_in_gene:
                        d = gene_to_dist.get(s[0] + '\t' + s[1])
                        if d is None:
                            gene_to_dist[s[0] + '\t' + s[1]] = dist
                        else:
                            gene_to_dist[s[0] + '\t' + s[1]] = np.sqrt(d*d + dist*dist)
                else:
                    res.add(var)
                    if merge_all_mut_in_gene:
                        d = gene_to_dist.get(s[0] + '\t' + s[1])
                        if s[5] == 'ins':
                            c = len(s[4]) * indel_cost*indel_cost
                        else:
                            c = indel_cost*indel_cost
                        if d is None:
                            gene_to_dist[s[0] + '\t' + s[1]] = np.sqrt(c)
                        else:
                            gene_to_dist[s[0] + '\t' + s[1]] = np.sqrt(d*d + c)
                    if merge_long_dels and s[5] == 'del':
                        del_list = deletions.get(s[0] + '\t' + s[1])
                        if del_list is None:
                            del_list = [[int(s[2]), 1]]
                            deletions[s[0] + '\t' + s[1]] = del_list
                        else:
                            pos = int(s[2])
                            for l in del_list:
                                if pos == l[0] + l[1]:
                                    l[1] += 1
                                    break
            else:
                res.add(var)
                if s[0] != 'non_cds' and merge_all_mut_in_gene:
                    d = gene_to_dist.get(s[0] + '\t' + s[1])
                    if s[5] == 'ins':
                        c = len(s[4])
                    else:
                        c = 1
                    if d is None:
                        gene_to_dist[s[0] + '\t' + s[1]] = c
                    else:
                        gene_to_dist[s[0] + '\t' + s[1]] = d + c
                if merge_long_dels and s[5] == 'del':
                    del_list = deletions.get(s[0] + '\t' + s[1])
                    if del_list is None:
                        del_list = [[int(s[2]), 1]]
                        deletions[s[0] + '\t' + s[1]] = del_list
                    else:
                        pos = int(s[2])
                        for l in del_list:
                            if pos == l[0] + l[1]:
                                l[1] += 1
                                break
                if merge_all_mut_in_pos:
                    res.add('\t'.join(s[:-3]))
    for gene_name, del_list in deletions.items():
        for long_del in del_list:
            if long_del[1] > 1:
                res.add(gene_name + '\t' + str(long_del[0]) + '\tlen=\t' + str(long_del[1]) + '\tdel')
    res_list = []
    for var in res:
        s = var.split('\t')
        if s[0] + '\t' + s[1] not in broken_genes:
            res_list.append(var)
    for gene_name in broken_genes:
        res_list.append(gene_name + '\tbroken')
    for gene_name, d in gene_to_dist.items():
        if gene_name not in broken_genes:
            res_list.append(gene_name + '\t' + str(d) + '\tdist')
    res_list.sort()
    return sample_id, res_list


def transpose_snp_matrix(sample_to_snps, pheno):
    snp_to_samples = {}
    snp_to_resistant_count = {}
    for sample_id, snp_list in sample_to_snps.items():
        if pheno[sample_id] == 1:
            is_resistant = True
        else:
            is_resistant = False
        for snp in snp_list:
            sample_list = snp_to_samples.get(snp)
            if sample_list is None:
                snp_to_samples[snp] = [sample_id]
                if is_resistant:
                    snp_to_resistant_count[snp] = 1
                else:
                    snp_to_resistant_count[snp] = 0
            else:
                sample_list.append(sample_id)
                if is_resistant:
                    snp_to_resistant_count[snp] += 1

    return snp_to_samples, snp_to_resistant_count


def gen_pairs_for_snp(snp_pivot_list, sample_to_snps, snp_to_samples, sample_num, pheno, snp_to_resistant_count):
    res = []
    for snp_pivot in snp_pivot_list:
        snp_pair_to_count = {}
        snp_pair_to_r_count = {}
        pivot_r_count = snp_to_resistant_count[snp_pivot]
        if pivot_r_count <= snp_count_r_threshold:
            continue
        pivot_c = len(snp_to_samples[snp_pivot])
        for sample_id in snp_to_samples[snp_pivot]:
            snp_list = sample_to_snps[sample_id]
            if pheno[sample_id] == 1:
                is_resistant = True
            else:
                is_resistant = False
            for snp in snp_list:
                if snp_pivot < snp:
                    pair = (snp_pivot, snp)
                    c = snp_pair_to_count.get(pair)
                    if c is None:
                        snp_pair_to_count[pair] = 1
                        if is_resistant:
                            snp_pair_to_r_count[pair] = 1
                        else:
                            snp_pair_to_r_count[pair] = 0
                    else:
                        snp_pair_to_count[pair] = c + 1
                        if is_resistant:
                            snp_pair_to_r_count[pair] += 1

        for pair, c in snp_pair_to_count.items():
            if snp_count_threshold < c < sample_num - snp_count_threshold:
                snp_c = len(snp_to_samples[pair[1]])
                if pivot_c - c > snp_count_threshold and snp_c - c > snp_count_threshold:
                    if c/(pivot_c + snp_c - c) < jaccard_index_threshold:
                        snp_r_c = snp_to_resistant_count[pair[1]]
                        if snp_r_c <= snp_count_r_threshold:
                            continue
                        pair_r_c = snp_pair_to_r_count[pair]
                        if pair_r_c * pivot_c * snp_c != pivot_r_count * snp_r_c * c:
                            res.append(pair)
    return res


def gen_all_pairs(sample_to_snps, all_snp_list, pheno):
    snp_to_samples, snp_to_resistant_count = transpose_snp_matrix(sample_to_snps, pheno)
    sample_num = len(sample_to_snps)
    # thread_num = os.cpu_count()
    snp_lists = {}
    for i in range(thread_num):
        snp_lists[i] = []
    for i in range(len(all_snp_list)):
        snp_lists[i % thread_num].append(all_snp_list[i])
    tasks = Parallel(n_jobs=thread_num)(delayed(gen_pairs_for_snp)(snp_list, sample_to_snps, snp_to_samples, sample_num,
                                                                   pheno, snp_to_resistant_count)
                                for snp_list in snp_lists.values())
    pairs = set()
    for task in tasks:
        for pair in task:
            pairs.add(pair[0] + '_' + pair[1])
    return pairs


def add_pairs_to_all_snp_lists(sample_to_snps, filtered_pairs):
    for sample_id, snp_list in sample_to_snps.items():
        snp_num = len(snp_list)
        for i in range(snp_num):
            for j in range(i + 1, snp_num):
                pair = snp_list[i] + '_' + snp_list[j]
                if pair in filtered_pairs:
                    snp_list.append(pair)
    return sample_to_snps


def update_snp_lists(sample_to_snps, pairs):
    sample_ids = list(sample_to_snps.keys())
    sample_to_snps_array = {}
    for i in range(thread_num):
        sample_to_snps_array[i] = {}
    for i in range(len(sample_ids)):
        sample_id = sample_ids[i]
        sample_to_snps_array[i % thread_num][sample_id] = sample_to_snps[sample_id]

    tasks1 = Parallel(n_jobs=thread_num)(delayed(add_pairs_to_all_snp_lists)(sample_to_snps_i, pairs)
                      for sample_to_snps_i in sample_to_snps_array.values())

    for task in tasks1:
        for sample_id, snp_list in task.items():
            sample_to_snps[sample_id] = snp_list


def main():
    if not exists(out_path):
        makedirs(out_path)
    drug_to_pheno = {}
    parallel = Parallel(n_jobs=-1)
    for (dirpath, dirnames, filenames) in os.walk(path_to_pheno):
        tasks = parallel(delayed(read_pheno)(path_to_pheno, drug) for drug in [filename[:-6] for filename in filenames])
        for drug, sample_id_to_pheno in tasks:
            drug_to_pheno[drug] = sample_id_to_pheno
    print(str(len(drug_to_pheno.keys())) + ' drugs')

    set_name_to_list = {}
    for (dirpath, dirnames, filenames) in os.walk(path_to_subsets):
        tasks = parallel(delayed(read_subset)(path_to_subsets, name) for name in [filename[:-4] for filename in filenames])
        for name, l in tasks:
            set_name_to_list[name] = l

    with open(path_to_ids, 'r') as f:
        sample_ids = [name.strip() for name in f.readlines()]
    print(str(len(sample_ids)) + ' samples')

    sample_to_variants = {}
    all_snps_set = set()
    tasks = parallel(delayed(read_variants)(path_to_snp, sample_id) for sample_id in sample_ids)
    matrix = create_distance_matrix(MatrixInfo.blosum62)
    indel_cost = max(matrix.values())
    tasks1 = parallel(delayed(process_snps)(sample_id, snp_list, matrix, indel_cost) for sample_id, snp_list in tasks)
    tasks = tasks1
    for name, l in tasks:
        sample_to_variants[name] = l
        for snp in l:
            s = snp.split('\t')
            if s[0] == 'Gene':
                if s[-1] == 'snp' or s[-1] == 'dist':
                    all_snps_set.add('\t'.join(s[:-2]))
                else:
                    all_snps_set.add(snp)
            else:
                if s[-1] == 'dist':
                    all_snps_set.add('\t'.join(s[:-2]))
                else:
                    all_snps_set.add(snp)
    print('sample snp reading done')

    tasks = parallel(delayed(split_dataset)(drug, pheno, sample_to_variants, set_name_to_list['Walker_subset'])
                     for drug, pheno in drug_to_pheno.items())

    drug_to_splits = {}
    for drug, snp_train, snp_test, pheno_train, pheno_test in tasks:
        drug_to_splits[drug] = (snp_train, snp_test, pheno_train, pheno_test)

    tasks1 = parallel(delayed(filter_snps)(drug, snp_train, snp_test)
                      for drug, snp_train, snp_test, pheno_train, pheno_test in tasks)

    for drug, snp_train_filtered, snp_test_filtered, full_snp_list in tasks1:
        snp_train, snp_test, pheno_train, pheno_test = drug_to_splits[drug]
        drug_to_splits[drug] = (snp_train_filtered, snp_test_filtered, pheno_train, pheno_test, full_snp_list)
        print(drug + ': ' + str(len(full_snp_list)) + ' snps after filtering')

    if use_pairs:
        for drug, (snp_train, snp_test, pheno_train, pheno_test, full_snp_list) in drug_to_splits.items():
            pairs = gen_all_pairs(snp_train, full_snp_list, pheno_train)
            print(drug + ' ' + str(len(pairs)) + ' pairs generated')

            if len(pairs) > 0:
                for pair in pairs:
                    full_snp_list.append(pair)
                print(drug + ': ' + str(len(full_snp_list)) + ' snps + pairs after filtering')
                update_snp_lists(snp_train, pairs)
                update_snp_lists(snp_test, pairs)
                print('snp lists extended')

    if only_print_features:
        if not exists(transformed_out_path):
            makedirs(transformed_out_path)
        with open(path_to_transformed_list, 'w') as f:
            for var in all_snps_set:
                f.write(var)
                f.write('\n')
        for sample_id, var_list in sample_to_variants.items():
            with open(transformed_out_path + sample_id + '.variants', 'w') as f:
                for var in var_list:
                    f.write(var)
                    f.write('\n')
    else:
        with open(out_path, 'a') as f:
            # f.write('MI FS k=10000 LR_l1 C1\n\n')
            # clf = MIFSLR(penalty='l1', dual=False, class_weight='balanced', C=1, k=10000)  + pairs
            f.write('LR_l1 C1 transformed extended feature set filtered indels\n\n')
            clf = LogisticRegression(penalty='l1', dual=False, class_weight='balanced', C=1)
            # f.write('RF n=10000\n\n')
            # clf = RandomForestClassifier(n_estimators=10000, n_jobs=-1, class_weight='balanced')
            # f.write('MI GFFS mi_k=10000fs_k=1000 LR_l2 C1\n\n')
            # clf = MIGFFSLR(penalty='l2', dual=False, class_weight='balanced', C=1, mi_k=10000, fs_k=1000, n_jobs=-1)
            print_stat_fixed_multi(f, drug_to_splits, clf, n_jobs=-1)
            f.write('\n')


if __name__ == '__main__':
    main()
