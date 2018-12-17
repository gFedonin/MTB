from scipy.sparse import lil_matrix
from sklearn.externals.joblib import Parallel, delayed
from sklearn.feature_selection import RFECV
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, make_scorer
from sklearn.model_selection import cross_validate, cross_val_score
import numpy as np
import os
from mlxtend.feature_selection import SequentialFeatureSelector as SFS

from src.core.constants import custom_scoring, tp, tn, fp, fn
from src.core.data_reading import read_pheno, read_subset

path_to_pheno = './data/pheno/'
path_to_ids = './data/dr_covered_with_pheno_and_snp.txt'
path_to_subsets = './data/subsets/'
path_to_tree_features = './data/tree_features_Walker/'

out_path = './res/tree_features_LR_l1_C1_multi.out'

scoring = 'roc_auc'


def stat_cv_multi(X, y, clf):
    cv_results = cross_validate(clf, X, y, scoring=custom_scoring, n_jobs=1, return_train_score=False)
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


def stat_fixed_multi(X, sample_indexes, y, subset_indices, clf):
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
    X_train = X[train_indices, :]
    X_test = X[test_indices, :]
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


def drug_stat(pheno, X, sample_id_to_index, clf):
    sample_indexes = []
    y = []
    for sample_id, val in pheno:
        if sample_id in sample_id_to_index.keys():
            sample_indexes.append(sample_id_to_index[sample_id])
            y.append(val)
    sum = np.sum(y)
    samples_with_pheno = len(sample_indexes)
    if samples_with_pheno == 0 or sum == len(y) or sum == 0:
        return samples_with_pheno, 0
    else:
        # print(str(len(sample_indexes)) + ' samples with pheno')
        X_slice = X[sample_indexes, :]
        return samples_with_pheno, cross_val_score(clf, X_slice, y, cv=10, n_jobs=1, scoring=scoring).mean()


def drug_stat_multi(pheno, X, sample_id_to_index, clf, subset=None):
    sample_indexes = []
    y = []
    for sample_id, val in pheno:
        if sample_id in sample_id_to_index.keys():
            sample_indexes.append(sample_id_to_index[sample_id])
            y.append(val)
    subset_indices = []
    if subset is not None:
        for sample_id in subset:
            if sample_id in sample_id_to_index.keys():
                subset_indices.append(sample_id_to_index[sample_id])
    sum = np.sum(y)
    samples_with_pheno = len(sample_indexes)
    if samples_with_pheno == 0 or sum == len(y) or sum == 0:
        return samples_with_pheno, 0, 0, 0, 0
    else:
        if subset is None:
            X_slice = X[sample_indexes, :]
            ppv, npv, sensitivity, specificity = stat_cv_multi(X_slice, y, clf)
        else:
            ppv, npv, sensitivity, specificity = stat_fixed_multi(X, sample_indexes, y, subset_indices, clf)
        return samples_with_pheno, ppv, npv, sensitivity, specificity



def print_stat_multi(f, drug_to_pheno, X, sample_id_to_index, clf):
    f.write('Drug\tSamples\tPPV\tNPV\tSensitivity\tSpecificity\n')
    for drug, pheno in drug_to_pheno.items():
        f.write(drug)
        samples_with_pheno, ppv, npv, sensitivity, specificity = drug_stat_multi(pheno, X, sample_id_to_index, clf)
        f.write('\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n' % (samples_with_pheno, ppv, npv, sensitivity, specificity))


def print_stat_fixed_multi(f, drug_to_pheno, X, sample_id_to_index, clf, subset):
    f.write('Drug\tSamples\tPPV\tNPV\tSensitivity\tSpecificity\n')
    for drug, pheno in drug_to_pheno.items():
        f.write(drug)
        samples_with_pheno, ppv, npv, sensitivity, specificity = drug_stat_multi(pheno, X, sample_id_to_index, clf, subset)
        f.write('\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n' % (samples_with_pheno, ppv, npv, sensitivity, specificity))


def read_tree_features(filename):
    l = []
    with open(path_to_tree_features + filename, 'r') as f:
        for line in f.readlines():
            # s = line.strip().split(',')
            l.append(line.strip())
    return filename, l



def main():
    drug_to_pheno = {}
    for (dirpath, dirnames, filenames) in os.walk(path_to_pheno):
        tasks = Parallel(n_jobs=-1)(delayed(read_pheno)(path_to_pheno, drug) for drug in [filename[:-6] for filename in filenames])
        for drug, sample_id_to_pheno in tasks:
            drug_to_pheno[drug] = sample_id_to_pheno
    print(str(len(drug_to_pheno.keys())) + ' drugs')

    set_name_to_list = {}
    for (dirpath, dirnames, filenames) in os.walk(path_to_subsets):
        tasks = Parallel(n_jobs=-1)(delayed(read_subset)(path_to_subsets, name) for name in [filename[:-4] for filename in filenames])
        for name, l in tasks:
            set_name_to_list[name] = l

    drug_to_tree_features = {}
    for (dirpath, dirnames, filenames) in os.walk(path_to_tree_features):
        tasks = Parallel(n_jobs=-1)(delayed(read_tree_features)(filename) for filename in filenames)
        for filename, l in tasks:
            name = filename[0:filename.find('.')]
            drug_to_tree_features[name] = l

    feature_to_index = {}
    j = 0
    for l in drug_to_tree_features.values():
        for s in l:
            if s not in feature_to_index.keys():
                feature_to_index[s] = j
                j += 1


    with open(path_to_ids, 'r') as f:
        sample_ids = [name.strip() for name in f.readlines()]
    print(str(len(sample_ids)) + ' samples')

    sample_num = len(sample_ids)
    sample_id_to_index = {}
    for i in range(sample_num):
        sample_id_to_index[sample_ids[i]] = i
    X = lil_matrix((sample_num, len(feature_to_index)), dtype=int)
    for l in drug_to_tree_features.values():
        for s in l:
            samples = s.split(',')
            f = feature_to_index[s]
            for sample_id in samples:
                if sample_id in sample_id_to_index.keys():
                    X[sample_id_to_index[sample_id], f] = 1

    print('sparse matrix created')

    with open(out_path, 'w') as f:
        print_func = print_stat_fixed_multi
        # f.write('LR_l2 + FS LR_l1 C1\n\n')
        # clf = LRFSLR(dual=False, class_weight='balanced', C=1)
        # print_func(f, drug_to_pheno, X, sample_id_to_index, clf=clf)
        # f.write('\n')
        f.write('LR_l1 C1\n\n')
        clf = LogisticRegression(penalty='l1', dual=False, class_weight='balanced', C=1)
        print_func(f, drug_to_pheno, X, sample_id_to_index, clf=clf, subset=set_name_to_list['Walker_subset'])
        f.write('\n')
        # f.write('LR_l2 C1 + SFS\n\n')
        # clf = RFECV(LogisticRegression(penalty='l1', dual=False, class_weight='balanced'), cv=3)
        # clf = SFS(LogisticRegression(dual=False, class_weight='balanced', C=1),
        #           k_features=(1, snp_num),
        #           forward=True,
        #           floating=False,
        #           scoring='accuracy',
        #           cv=3,
        #           n_jobs=-1)
        # print_func(f, drug_to_pheno, X, sample_id_to_index, clf=clf)
        # f.write('\n')


if __name__ == '__main__':
    main()