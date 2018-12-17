from scipy.sparse import lil_matrix
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals.joblib import Parallel, delayed
import numpy as np
import os

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score, confusion_matrix, roc_auc_score
from sklearn.neighbors import KNeighborsClassifier

from src.core.annotations import read_annotations, CDSType
from src.core.constants import upstream_length
from src.core.data_reading import read_pheno, read_subset, read_variants, read_snp_list
from src.data_processing.print_mutations_extended import process_variants


path_to_pheno = '../../data/pheno/'
path_to_var = '../../data/snps/annotated_with_DR_with_indel_with_pheno_and_snp_mc10/'
path_to_ids = '../../data/dr_covered_with_pheno_and_snp.txt'
path_to_subsets = '../../data/subsets/'

out_path = '../../res/known_genes_ml.out'

use_extended_features = False
merge_all_mut_in_pos = False
merge_all_mut_in_gene = False
filter_non_cds = False
upstream_indel_breakes_gene = True
snp_count_threshold = 3

use_pairs = False
jaccard_index_threshold = 0.9

thread_num = 32

filter_by_var_list = False
path_to_var_list = path_to_var + 'dr_genes_var_list.csv'


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
    cm = confusion_matrix(y_test, y_pred)
    tp_score = cm[1, 1]
    tn_score = cm[0, 0]
    fp_score = cm[0, 1]
    fn_score = cm[1, 0]
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
    return ppv, npv, sensitivity, specificity, f1_score(y_test, y_pred), roc_auc_score(y_test, y_pred)


def snps_counts(sample_to_snps):
    snp_to_count = {}
    for l in sample_to_snps.values():
        for snp in l:
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
            if snp_count_threshold < snp_to_count[snp] < sample_num - snp_count_threshold:
                res.append(snp)
                all_snps_set.add(snp)
        snp_train_filtered[sample_id] = res
    for sample_id, snp_list in snp_test.items():
        snp_test_filtered[sample_id] = [snp for snp in snp_list if snp in all_snps_set]
    full_snp_list = list(all_snps_set)
    return drug, snp_train_filtered, snp_test_filtered, full_snp_list


def drug_stat_multi(drug, snp_train, snp_test, pheno_train, pheno_test, clf):
    if len(pheno_train) == 0:
        return drug, len(pheno_train) + len(pheno_test), 0, 0, 0, 0, 0, 0

    full_snp_set = set()
    for snp_list in snp_train.values():
        for snp in snp_list:
            full_snp_set.add(snp)
    full_snp_list = list(full_snp_set)

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
            X_train[i, feature_to_index[f]] = 1
        y_train.append(pheno_train[sample_ids[i]])

    sample_ids = list(snp_test.keys())
    sample_num = len(sample_ids)
    X_test = lil_matrix((sample_num, snp_num), dtype=int)
    y_test = []
    for i in range(len(sample_ids)):
        for f in snp_test[sample_ids[i]]:
            X_test[i, feature_to_index[f]] = 1
        y_test.append(pheno_test[sample_ids[i]])

    ppv, npv, sensitivity, specificity, f1, auc = stat_fixed_multi(X_train, y_train, X_test, y_test, clf)
    return drug, len(pheno_train) + len(pheno_test), ppv, npv, sensitivity, specificity, f1, auc


def print_stat_fixed_multi(f, drug_to_splits, clf, n_jobs=-1):
    f.write('Drug\tSamples\tPPV\tNPV\tSensitivity\tSpecificity\tF1\tAUC\n')
    if n_jobs != 1:
        tasks = Parallel(n_jobs=n_jobs)(delayed(drug_stat_multi)(drug, snp_train, snp_test, pheno_train, pheno_test, clf)
                                        for drug, (snp_train, snp_test, pheno_train, pheno_test, full_snp_list)
                                        in drug_to_splits.items())
        for drug, samples_with_pheno, ppv, npv, sensitivity, specificity, f1, auc in tasks:
            f.write(drug)
            f.write('\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n' % (samples_with_pheno, ppv, npv, sensitivity, specificity, f1, auc))
    else:
        for drug, (snp_train, snp_test, pheno_train, pheno_test, full_snp_list) in drug_to_splits.items():
            drug, samples_with_pheno, ppv, npv, sensitivity, specificity, f1, auc = \
                drug_stat_multi(drug, snp_train, snp_test, pheno_train, pheno_test, clf)
            f.write(drug)
            f.write('\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n' % (samples_with_pheno, ppv, npv, sensitivity, specificity, f1, auc))


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
    if filter_by_var_list:
        full_snp_list = read_snp_list(path_to_var_list)
        tasks = parallel(delayed(read_variants)(path_to_var, sample_id, full_snp_list) for sample_id in sample_ids)
    else:
        tasks = parallel(delayed(read_variants)(path_to_var, sample_id) for sample_id in sample_ids)
    if use_extended_features:
        cds_list = read_annotations(upstream_length)
        name_to_type = {}
        for cds in cds_list:
            if cds.type != CDSType.upstream:
                name_to_type[cds.name] = cds.type
        tasks1 = parallel(delayed(process_variants)(sample_id, snp_list, name_to_type, merge_all_mut_in_pos,
                                               merge_all_mut_in_gene, filter_non_cds, upstream_indel_breakes_gene)
                          for sample_id, snp_list in tasks)
        tasks = tasks1
    for name, l in tasks:
        sample_to_variants[name] = l
        for snp in l:
            all_snps_set.add(snp)
    print('samples\' variants reading done')

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
        print(drug + ': ' + str(len(full_snp_list)) + ' variants after filtering')

    if use_pairs:
        for drug, (snp_train, snp_test, pheno_train, pheno_test, full_snp_list) in drug_to_splits.items():
            pairs = gen_all_pairs(snp_train, full_snp_list, pheno_train)
            print(drug + ' ' + str(len(pairs)) + ' pairs generated')

            if len(pairs) > 0:
                for pair in pairs:
                    full_snp_list.append(pair)
                print(drug + ': ' + str(len(full_snp_list)) + ' variants + pairs after filtering')
                update_snp_lists(snp_train, pairs)
                update_snp_lists(snp_test, pairs)
                print('snp lists extended')

    with open(out_path, 'a') as f:
        if filter_by_var_list:
            f.write('known DR genes only\n')
        if filter_non_cds:
            f.write('non cds filtered out\n')
        if use_extended_features:
            f.write('generated features added:\n')
            f.write('gene broken\n')
            if merge_all_mut_in_pos:
                f.write('all mut in pos merged\n')
            if merge_all_mut_in_gene:
                f.write('gene changed\n')
            if upstream_indel_breakes_gene:
                f.write('upstream indel breakes gene\n')
        f.write('snp_count_threshold = %d\n' % snp_count_threshold)

        if use_pairs:
            f.write('pairs added\n')
            f.write('jaccard_index_threshold = %1.2f\n' % jaccard_index_threshold)


        # f.write('LR_l1 C1\n\n')
        # clf = LogisticRegression(penalty='l1', dual=False, class_weight='balanced', C=1)
        f.write('NN k = 1\n\n')
        clf = KNeighborsClassifier(n_neighbors=1, p=1, n_jobs=-1)
        # f.write('RF n=10000\n\n')
        # clf = RandomForestClassifier(n_estimators=10000, n_jobs=-1, class_weight='balanced')
        # f.write('MI GFFS mi_k=10000fs_k=1000 LR_l2 C1\n\n')
        # clf = MIGFFSLR(penalty='l2', dual=False, class_weight='balanced', C=1, mi_k=10000, fs_k=1000, n_jobs=-1)
        print_stat_fixed_multi(f, drug_to_splits, clf, n_jobs=1)
        f.write('\n')



if __name__ == '__main__':
    main()
