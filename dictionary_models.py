import os

import numpy as np

from scipy.sparse import lil_matrix
from sklearn.externals.joblib import Parallel, delayed
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import auc, roc_curve, f1_score
from sklearn.model_selection import cross_val_score, cross_validate

# from lrfslr import LRFSLR
from vocabulary_clf import VocabularyClf
from constants import custom_scoring, tp, tn, fp, fn
from data_reading import read_pheno, read_dict, read_subset, read_variants

path_to_pheno = './data/pheno_mc5/'
path_to_snp = './data/snps/annotated_with_pheno_and_snp_mc5/'
path_to_snp_list = path_to_snp + 'all_snp_list.csv'
path_to_dict = './data/dictionaries1/'
path_to_subsets = './data/subsets/'
path_to_ids = './data/dr_covered_with_pheno_and_snp.txt'
out_path = './res/dict_simple_multi_Walker_train_rest_test.out'

# path_to_selected_features = './res/LRFSLR_features.txt'

scoring = 'roc_auc'
threads = 10


def list_stat_cv_multi(X, y, f_list, feature_to_index, clf):
    f_indices = []
    for f in f_list:
        if f in feature_to_index.keys():
            f_indices.append(feature_to_index[f])
    # print(str(len(f_indices)) + ' features')
    X_slice = X[:, f_indices]
    cv_results = cross_validate(clf, X_slice, y, scoring=custom_scoring, n_jobs=1, return_train_score=False)
    tp = cv_results['test_tp'].sum()
    tn = cv_results['test_tn'].sum()
    fp = cv_results['test_fp'].sum()
    fn = cv_results['test_fn'].sum()
    if tp + fp > 0:
        ppv = tp/(tp + fp)
    else:
        ppv = 0
    if tn + fn > 0:
        npv = tn/(tn + fn)
    else:
        npv = 0
    if tp + fn > 0:
        sensitivity = tp/(tp + fn)
    else:
        sensitivity = 0
    if tn + fp > 0:
        specificity = tn/(tn + fp)
    else:
        specificity = 0
    return ppv, npv, sensitivity, specificity


def list_stat_cv_score(X, y, f_list, feature_to_index, clf):
    f_indices = []
    for f in f_list:
        if f in feature_to_index.keys():
            f_indices.append(feature_to_index[f])
    # print(str(len(f_indices)) + ' features')
    X_slice = X[:, f_indices]
    return cross_val_score(clf, X_slice, y, cv=10, n_jobs=threads, scoring=scoring).mean()


def list_stat_score(X, y, f_list, feature_to_index):
    clf = VocabularyClf()
    f_indices = []
    for f in f_list:
        if f in feature_to_index.keys():
            f_indices.append(feature_to_index[f])
    # print(str(len(f_indices)) + ' features')
    X_slice = X[:, f_indices]
    clf.fit(X_slice, y)
    if scoring == 'roc_auc':
        fpr, tpr, thresholds = roc_curve(y, clf.predict_proba(X_slice)[:, 1])
        return auc(fpr, tpr)
    elif scoring == 'f1':
        return f1_score(y, clf.predict(X_slice))


def list_stat_multi(X, y, f_list, feature_to_index):
    clf = VocabularyClf()
    f_indices = []
    for f in f_list:
        if f in feature_to_index.keys():
            f_indices.append(feature_to_index[f])
    # print(str(len(f_indices)) + ' features')
    X_slice = X[:, f_indices]
    clf.fit(X_slice, y)
    pred = clf.predict(X_slice)
    tp_score = tp(y, pred)
    tn_score = tn(y, pred)
    fp_score = fp(y, pred)
    fn_score = fn(y, pred)

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
    return ppv, npv, sensitivity, specificity


def drug_stat(drug, pheno, dict_name_to_list, X, sample_id_to_index, feature_to_index, subset=None, clf=None):
    sample_indexes = []
    y = []
    if subset is None:
        for sample_id, val in pheno:
            if sample_id in sample_id_to_index.keys():
                sample_indexes.append(sample_id_to_index[sample_id])
                y.append(val)
    else:
        for sample_id, val in pheno:
            if sample_id in subset and sample_id in sample_id_to_index.keys():
                sample_indexes.append(sample_id_to_index[sample_id])
                y.append(val)
    sum = np.sum(y)
    samples_with_pheno = len(sample_indexes)
    dict_name_to_stat = {}
    if samples_with_pheno == 0 or sum == len(y) or sum == 0:
        for dict_name, f_list in dict_name_to_list.items():
            dict_name_to_stat[dict_name] = (0, 0)
    else:
        # print(str(len(sample_indexes)) + ' samples with pheno')
        X_slice = X[sample_indexes, :]
        for dict_name, f_list in dict_name_to_list.items():
            # print(dict_name)
            if drug in f_list.keys():
                features = f_list[drug]
                if clf is not None:
                    dict_name_to_stat[dict_name] = (len(features), list_stat_cv_score(X_slice, y, features, feature_to_index, clf))
                else:
                    dict_name_to_stat[dict_name] = (len(features), list_stat_score(X_slice, y, features, feature_to_index))
            else:
                dict_name_to_stat[dict_name] = (0, 0)
    return samples_with_pheno, dict_name_to_stat


def drug_stat_multi(drug, pheno, dict_name_to_list, X, sample_id_to_index, feature_to_index, subset=None, clf=None):
    sample_indexes = []
    y = []
    if subset is None:
        for sample_id, val in pheno:
            if sample_id in sample_id_to_index.keys():
                sample_indexes.append(sample_id_to_index[sample_id])
                y.append(val)
    else:
        for sample_id, val in pheno:
            if sample_id in subset and sample_id in sample_id_to_index.keys():
                sample_indexes.append(sample_id_to_index[sample_id])
                y.append(val)
    sum = np.sum(y)
    samples_with_pheno = len(sample_indexes)
    dict_name_to_stat = {}
    if samples_with_pheno == 0 or sum == len(y) or sum == 0:
        for dict_name, f_list in dict_name_to_list.items():
            dict_name_to_stat[dict_name] = (0, 0, 0, 0, 0)
    else:
        # print(str(len(sample_indexes)) + ' samples with pheno')
        X_slice = X[sample_indexes, :]
        for dict_name, f_list in dict_name_to_list.items():
            # print(dict_name)
            if drug in f_list.keys():
                features = f_list[drug]
                if clf is not None:
                    # with open(path_to_selected_features, 'a') as f:
                    #     f.write(drug + '\n')
                    #     f.write(';'.join(features) + '\n')
                    ppv, npv, sensitivity, specificity = list_stat_cv_multi(X_slice, y, features, feature_to_index, clf)
                else:
                    ppv, npv, sensitivity, specificity = list_stat_multi(X_slice, y, features, feature_to_index)
                dict_name_to_stat[dict_name] = (len(features), ppv, npv, sensitivity, specificity)
            else:
                dict_name_to_stat[dict_name] = (0, 0, 0, 0, 0)
    return samples_with_pheno, dict_name_to_stat


def print_stat(f, dict_name_to_list, drug_to_pheno, X, sample_id_to_index, feature_to_index, dict_name_to_counts,
               subset=None, clf=None):
    f.write('Drug\tSamples')
    for name in dict_name_to_list.keys():
        f.write('\t' + name)
    f.write('\n')
    for drug, pheno in drug_to_pheno.items():
        f.write(drug)
        samples_with_pheno, dict_name_to_stat = drug_stat(drug, pheno, dict_name_to_list, X, sample_id_to_index,
                                                          feature_to_index, subset, clf)
        f.write('\t' + str(samples_with_pheno))
        for name in dict_name_to_list.keys():
            feature_num, stat = dict_name_to_stat[name]
            try:
                f_count = dict_name_to_counts[name][drug]
            except:
                f_count = 0
            f.write('\t%d/%d(%1.2f)' % (f_count, feature_num, stat))
        f.write('\n')


def print_stat_multi(f, dict_name_to_list, drug_to_pheno, X, sample_id_to_index, feature_to_index, dict_name_to_counts,
                     subset=None, clf=None):
    f.write('Drug\tSamples')
    for name in dict_name_to_list.keys():
        f.write('\t' + name)
    f.write('\n')
    for drug, pheno in drug_to_pheno.items():
        f.write(drug)
        samples_with_pheno, dict_name_to_stat = drug_stat_multi(drug, pheno, dict_name_to_list, X, sample_id_to_index,
                                                                feature_to_index, subset, clf)
        f.write('\t' + str(samples_with_pheno))
        for name in dict_name_to_list.keys():
            feature_num, ppv, npv, sensitivity, specificity = dict_name_to_stat[name]
            try:
                f_count = dict_name_to_counts[name][drug]
            except:
                f_count = 0
            f.write('\t%d;%d;%1.2f;%1.2f;%1.2f;%1.2f' % (feature_num, f_count, ppv, npv, sensitivity, specificity))
        f.write('\n')


def count_snps_presence(feature_to_index, dict_name_to_list):
    dict_name_to_counts = {}
    for name, drug_to_mut_list in dict_name_to_list.items():
        drug_to_count = {}
        for drug, l in drug_to_mut_list.items():
            c = 0
            for f in l:
                if f in feature_to_index.keys():
                    c += 1
            drug_to_count[drug] = c
        dict_name_to_counts[name] = drug_to_count
    return dict_name_to_counts


def list_stat_fixed_multi(X, sample_indexes, y, f_list, feature_to_index, subset_indices, clf):
    f_indices = []
    for f in f_list:
        if f in feature_to_index.keys():
            f_indices.append(feature_to_index[f])
    # print(str(len(f_indices)) + ' features')
    train_indices = []
    test_indices = []
    y_train = []
    y_test = []
    for i in range(len(sample_indexes)):
        s = sample_indexes[i]
        if s in subset_indices:
            train_indices.append(s)
            y_train.append(y[i])
        else:
            test_indices.append(s)
            y_test.append(y[i])
    if len(train_indices) == 0:
        return 0, 0, 0, 0
    X_train = X[train_indices, :][:, f_indices]
    X_test = X[test_indices, :][:, f_indices]
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
    return ppv, npv, sensitivity, specificity


def drug_stat_fixed_multi(drug, pheno, dict_name_to_list, X, sample_id_to_index, feature_to_index, subset, clf):
    sample_indices = []
    y = []
    for sample_id, val in pheno:
        if sample_id in sample_id_to_index.keys():
            sample_indices.append(sample_id_to_index[sample_id])
            y.append(val)
    subset_indices = []
    for sample_id in subset:
        if sample_id in sample_id_to_index.keys():
            subset_indices.append(sample_id_to_index[sample_id])
    sum = np.sum(y)
    samples_with_pheno = len(sample_indices)
    dict_name_to_stat = {}
    if samples_with_pheno == 0 or sum == len(y) or sum == 0:
        for dict_name, f_list in dict_name_to_list.items():
            dict_name_to_stat[dict_name] = (0, 0, 0, 0, 0)
    else:
        # print(str(len(sample_indexes)) + ' samples with pheno')
        for dict_name, f_list in dict_name_to_list.items():
            # print(dict_name)
            if drug in f_list.keys():
                features = f_list[drug]
                ppv, npv, sensitivity, specificity = list_stat_fixed_multi(X, sample_indices, y, features,
                                                                           feature_to_index, subset_indices, clf)
                dict_name_to_stat[dict_name] = (len(features), ppv, npv, sensitivity, specificity)
            else:
                dict_name_to_stat[dict_name] = (0, 0, 0, 0, 0)
    return samples_with_pheno, dict_name_to_stat


def print_stat_fixed_multi(f, dict_name_to_list, drug_to_pheno, X, sample_id_to_index, feature_to_index, dict_name_to_counts,
                     subset=None, clf=None):
    f.write('Drug\tSamples')
    for name in dict_name_to_list.keys():
        f.write('\t' + name)
    f.write('\n')
    for drug, pheno in drug_to_pheno.items():
        f.write(drug)
        samples_with_pheno, dict_name_to_stat = drug_stat_fixed_multi(drug, pheno, dict_name_to_list, X, sample_id_to_index,
                                                                feature_to_index, subset, clf)
        f.write('\t' + str(samples_with_pheno))
        for name in dict_name_to_list.keys():
            feature_num, ppv, npv, sensitivity, specificity = dict_name_to_stat[name]
            try:
                f_count = dict_name_to_counts[name][drug]
            except:
                f_count = 0
            f.write('\t%d;%d;%1.2f;%1.2f;%1.2f;%1.2f' % (feature_num, f_count, ppv, npv, sensitivity, specificity))
        f.write('\n')


def main():

    drug_to_pheno = {}
    for (dirpath, dirnames, filenames) in os.walk(path_to_pheno):
        tasks = Parallel(n_jobs=-1)(delayed(read_pheno)(path_to_pheno, drug) for drug in [filename[:-6] for filename in filenames])
        for drug, sample_id_to_pheno in tasks:
            drug_to_pheno[drug] = sample_id_to_pheno
    print(str(len(drug_to_pheno.keys())) + ' drugs')

    dict_name_to_list = {}
    for (dirpath, dirnames, filenames) in os.walk(path_to_dict):
        tasks = Parallel(n_jobs=-1)(delayed(read_dict)(path_to_dict, name) for name in [filename[:-4] for filename in filenames])
        for name, drug_to_mut_list in tasks:
            dict_name_to_list[name] = drug_to_mut_list

    print(str(len(dict_name_to_list.keys())) + ' dictionaries')

    drug_to_full_list = {}
    for drug in drug_to_pheno.keys():
        full_set = set()
        for name, drug_to_mut_list in dict_name_to_list.items():
            if drug in drug_to_mut_list.keys():
                for f in drug_to_mut_list[drug]:
                    full_set.add(f)
        if len(full_set) > 0:
            drug_to_full_list[drug] = list(full_set)
    dict_name_to_list['All'] = drug_to_full_list
    full_dict_set = set()
    for l in drug_to_full_list.values():
        for f in l:
            full_dict_set.add(f)

    with open(path_to_ids, 'r') as f:
        sample_ids = [name.strip() for name in f.readlines()]

    set_name_to_list = {}
    for (dirpath, dirnames, filenames) in os.walk(path_to_subsets):
        tasks = Parallel(n_jobs=-1)(delayed(read_subset)(path_to_subsets, name) for name in [filename[:-4] for filename in filenames])
        for name, l in tasks:
            set_name_to_list[name] = l

    sample_to_snps = {}
    tasks = Parallel(n_jobs=-1)(delayed(read_variants)(path_to_snp, sample_id, full_dict_set) for sample_id in sample_ids)
    for name, l in tasks:
        sample_to_snps[name] = l

    print(str(len(sample_ids)) + ' samples after filtering')

    # full_snp_list = read_snp_list(path_to_snp_list)
    full_snp_list = list(full_dict_set)
    snp_num = len(full_snp_list)
    print(str(snp_num) + ' snps')

    feature_to_index = {}
    for i in range(snp_num):
        feature_to_index[full_snp_list[i]] = i

    dict_name_to_counts = count_snps_presence(feature_to_index, dict_name_to_list)

    sample_num = len(sample_to_snps)
    sample_id_to_index = {}
    for i in range(len(sample_ids)):
        sample_id_to_index[sample_ids[i]] = i
    X = lil_matrix((sample_num, snp_num), dtype=int)
    for i in range(len(sample_ids)):
        for f in sample_to_snps[sample_ids[i]]:
            if f == '':
                print(sample_ids[i])
                return 0
            X[i, feature_to_index[f]] = 1

    with open(out_path, 'w') as f:
        print_func = print_stat_fixed_multi
        f.write('FeatureNum/PPV/NPV/Sensitivity/Specificity\n\n')
        # f.write('LR_l2 + FS LR_l1 CV\n\n')
        # f.write('Full dataset:\n\n')
        # clf = LRFSLR(dual=False, class_weight='balanced', C=1)
        # f.write('LR_l2 + RFE CV\n\n')
        # f.write('Full dataset:\n\n')
        # clf = RFECV(LogisticRegression(penalty='l2', dual=False, class_weight='balanced'), cv=10)
        # print_func(f, dict_name_to_list, drug_to_pheno, X, sample_id_to_index, feature_to_index, dict_name_to_counts, clf=clf)
        # f.write('FeatureNum/PPV/NPV/Sensitivity/Specificity\n\n')
        # f.write('LinearSVC l2 CV\n\n')
        # f.write('Full dataset:\n\n')
        # clf = LinearSVC(penalty='l2', dual=False, class_weight='balanced', C=10)
        # print_func(f, dict_name_to_list, drug_to_pheno, X, sample_id_to_index, feature_to_index, dict_name_to_counts, clf=clf)
        # f.write('\n')
        # f.write('LinearSVC l1 CV\n\n')
        # f.write('Full dataset:\n\n')
        # clf = LinearSVC(penalty='l1', dual=False, class_weight='balanced', C=10)
        # print_func(f, dict_name_to_list, drug_to_pheno, X, sample_id_to_index, feature_to_index, dict_name_to_counts, clf=clf)
        # f.write('\n')
        # f.write('LogisticRegression l2 CV\n\n')
        # f.write('Full dataset:\n\n')
        # clf = LogisticRegression(penalty="l2", dual=False, class_weight='balanced', C=10)
        # print_func(f, dict_name_to_list, drug_to_pheno, X, sample_id_to_index, feature_to_index, dict_name_to_counts, clf=clf)
        # f.write('\n')
        # f.write('LogisticRegression l1 CV\n\n')
        # f.write('Full dataset:\n\n')
        # clf = LogisticRegression(penalty="l1", dual=False, class_weight='balanced', C=10)
        # print_func(f, dict_name_to_list, drug_to_pheno, X, sample_id_to_index, feature_to_index, dict_name_to_counts, clf=clf)
        # f.write('\n')
        # f.write('RandomForestClassifier n=1000 CV\n\n')
        # f.write('Full dataset:\n\n')
        # clf = RandomForestClassifier(n_estimators=1000, n_jobs=threads, class_weight='balanced')
        # print_func(f, dict_name_to_list, drug_to_pheno, X, sample_id_to_index, feature_to_index, dict_name_to_counts, clf=clf)
        # f.write('\n')
        # f.write('Simple Vocabulary\n\n')
        # f.write('Full dataset:\n\n')
        # print_func(f, dict_name_to_list, drug_to_pheno, X, sample_id_to_index, feature_to_index, dict_name_to_counts)
        # f.write('\n')
        # for name, subset in set_name_to_list.items():
        #     f.write(name + '\n\n')
        #     print_func(f, dict_name_to_list, drug_to_pheno, X, sample_id_to_index, feature_to_index, dict_name_to_counts, subset=subset)
        #     f.write('\n')
        f.write('Simple Dictionary\n\n')
        clf = VocabularyClf()
        f.write('train Walker dataset, test - the rest:\n\n')
        print_func(f, dict_name_to_list, drug_to_pheno, X, sample_id_to_index, feature_to_index, dict_name_to_counts,
                   subset=set_name_to_list['Walker_subset'], clf=clf)
        f.write('\n')
        f.write('LogisticRegression l1 CV\n\n')
        f.write('train Walker dataset, test - the rest:\n\n')
        clf = LogisticRegression(penalty="l1", dual=False, class_weight='balanced')
        print_func(f, dict_name_to_list, drug_to_pheno, X, sample_id_to_index, feature_to_index, dict_name_to_counts,
                   subset=set_name_to_list['Walker_subset'], clf=clf)
        f.write('\n')



if __name__ == '__main__':
    main()
