from os.path import exists

from scipy.sparse import lil_matrix
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals.joblib import Parallel, delayed
import os

from sklearn.linear_model import LogisticRegression

from src.core.annotations import read_annotations
from src.core.constants import data_path, upstream_length
from src.core.data_reading import read_pheno, read_subset, read_snp_list
from src.ml_methods.mtb_predictor import MTBPredictor, MTBPredictorExtended, MTBPredictorPairs

path_to_pheno = data_path + 'pheno_mc5_mega/'
path_to_var = data_path + 'snps/annotated_with_DR_with_indel_with_pheno_and_snp_mc10_long_del_pg/'
path_to_ids = data_path + 'dr_covered_with_pheno_and_snp_new.txt'
path_to_subsets = data_path + 'subsets/'

log_path = '../../res/ml_log_mc1-noncds_long_del_pg_nh_or_pv_ext_pairs/'

use_extended_features = True
merge_all_mut_in_pos = True
merge_all_mut_in_gene = True
filter_non_cds = True
drop_pseudogenes = True
drop_hypothetical = False
keep_proteomic_validated_only = False
keep_only_not_hypotetical_or_proteomic_validated = True
upstream_indel_breakes_gene = True
snp_count_threshold = 1

use_pairs = True
jaccard_index_threshold = 0.9

thread_num = 32

filter_by_var_list = False
path_to_var_list = path_to_var + 'dr_genes_var_list.csv'


def split_dataset(drug, pheno, sample_to_mut, train_subset):
    mut_train = {}
    mut_test = {}
    pheno_train = {}
    pheno_test = {}
    for sample_id, val in pheno:
        if sample_id in sample_to_mut.keys():
            if sample_id in train_subset:
                mut_train[sample_id] = sample_to_mut[sample_id]
                pheno_train[sample_id] = val
            else:
                mut_test[sample_id] = sample_to_mut[sample_id]
                pheno_test[sample_id] = val
    return drug, mut_train, mut_test, pheno_train, pheno_test


def drug_stat_multi(drug, snp_train, snp_test, pheno_train, pheno_test, clf):
    if len(pheno_train) == 0:
        return drug, len(pheno_train), len(pheno_test), 0, 0, 0, 0, 0, 0
    if use_extended_features:
        if use_pairs:
            clf = MTBPredictorPairs(drug, clf, snp_count_threshold, merge_all_mut_in_pos, merge_all_mut_in_gene,
                                    upstream_indel_breakes_gene, jaccard_index_threshold, thread_num,
                                    log_path + drug + '.log')
        else:
            clf = MTBPredictorExtended(drug, clf, snp_count_threshold, merge_all_mut_in_pos, merge_all_mut_in_gene,
                                       upstream_indel_breakes_gene, log_path + drug + '.log')
    else:
        clf = MTBPredictor(drug, clf, snp_count_threshold, log_path + drug + '.log')
    clf.train(snp_train, pheno_train)
    ppv, npv, sensitivity, specificity, f1, auc = clf.compute_stat(snp_test, pheno_test)
    clf.print_selected_features()
    return drug, len(pheno_train), len(pheno_test), ppv, npv, sensitivity, specificity, f1, auc


def print_stat_fixed_multi(f, drug_to_splits, clf, n_jobs=-1):
    f.write('Drug\tTrain\tTest\tPPV\tNPV\tSensitivity\tSpecificity\tF1\tAUC\n')
    if n_jobs != 1:
        tasks = Parallel(n_jobs=n_jobs)(delayed(drug_stat_multi)(drug, snp_train, snp_test, pheno_train, pheno_test, clf)
                                        for drug, (snp_train, snp_test, pheno_train, pheno_test)
                                        in drug_to_splits.items())
        for drug, train_num, test_num, ppv, npv, sensitivity, specificity, f1, auc in tasks:
            f.write(drug)
            f.write('\t%d\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n' % (train_num, test_num, ppv, npv,
                                                                              sensitivity, specificity, f1, auc))
    else:
        for drug, (snp_train, snp_test, pheno_train, pheno_test) in drug_to_splits.items():
            drug, train_num, test_num, ppv, npv, sensitivity, specificity, f1, auc = \
                drug_stat_multi(drug, snp_train, snp_test, pheno_train, pheno_test, clf)
            f.write(drug)
            f.write('\t%d\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n' % (train_num, test_num, ppv, npv,
                                                                              sensitivity, specificity, f1, auc))


def read_all_pheno():
    drug_to_pheno = {}
    for (dirpath, dirnames, filenames) in os.walk(path_to_pheno):
        tasks = Parallel(n_jobs=-1)(delayed(read_pheno)(path_to_pheno, drug) for drug in
                                    [filename[:-6] for filename in filenames])
        for drug, sample_id_to_pheno in tasks:
            drug_to_pheno[drug] = sample_id_to_pheno
    print(str(len(drug_to_pheno.keys())) + ' drugs')
    return drug_to_pheno


def read_all_subsets():
    set_name_to_list = {}
    for (dirpath, dirnames, filenames) in os.walk(path_to_subsets):
        tasks = Parallel(n_jobs=-1)(delayed(read_subset)(path_to_subsets, name) for name in
                                    [filename[:-4] for filename in filenames])
        for name, l in tasks:
            set_name_to_list[name] = l
    return set_name_to_list


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
        tasks = Parallel(n_jobs=-1)(delayed(read_variants)(path_to_var, sample_id, name_to_cds, full_snp_list)
                                    for sample_id in sample_ids)
    else:
        tasks = Parallel(n_jobs=-1)(delayed(read_variants)(path_to_var, sample_id, name_to_cds) for sample_id in sample_ids)
    for name, snp_list in tasks:
        sample_to_variants[name] = snp_list
    print('samples\' variants reading done')
    return sample_to_variants


def main():

    if not exists(log_path):
        os.makedirs(log_path)

    drug_to_pheno = read_all_pheno()

    set_name_to_list = read_all_subsets()

    with open(path_to_ids, 'r') as f:
        sample_ids = [name.strip() for name in f.readlines()]
    print(str(len(sample_ids)) + ' samples')

    sample_to_variants = read_all_variants(sample_ids)

    tasks = Parallel(n_jobs=-1)(delayed(split_dataset)(drug, pheno, sample_to_variants, set_name_to_list['Walker_subset'])
                     for drug, pheno in drug_to_pheno.items())

    drug_to_splits = {}
    for drug, snp_train, snp_test, pheno_train, pheno_test in tasks:
        drug_to_splits[drug] = (snp_train, snp_test, pheno_train, pheno_test)

    with open(log_path + 'ml.stats', 'w') as f:
        if filter_by_var_list:
            f.write('known DR genes only\n')
        if filter_non_cds:
            f.write('mut in non cds are filtered out\n')
        if drop_pseudogenes:
            f.write('mut in pseudogenes are filtered out\n')
        if drop_hypothetical:
            f.write('mut in hypothetical genes are filtered out\n')
        if keep_proteomic_validated_only:
            f.write('mut in not identified in proteomics study genes are filtered out\n')
        if keep_only_not_hypotetical_or_proteomic_validated:
            f.write('mut in hypothetical and not identified in proteomics study genes are filtered out\n')
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

        f.write('LR_l1 C1\n\n')
        clf = LogisticRegression(penalty='l1', dual=False, class_weight='balanced', C=1)
        # f.write('NN k = 1\n\n')
        # clf = KNeighborsClassifier(n_neighbors=1, p=1, n_jobs=-1)
        # f.write('RF n=10000\n\n')
        # clf = RandomForestClassifier(n_estimators=10000, n_jobs=-1, class_weight='balanced')
        # f.write('MI GFFS mi_k=10000fs_k=1000 LR_l2 C1\n\n')
        # clf = MIGFFSLR(penalty='l2', dual=False, class_weight='balanced', C=1, mi_k=10000, fs_k=1000, n_jobs=-1)
        print_stat_fixed_multi(f, drug_to_splits, clf, n_jobs=-1)
        f.write('\n')


if __name__ == '__main__':
    main()
