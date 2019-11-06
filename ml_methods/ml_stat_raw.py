from os.path import exists

from sklearn.ensemble import RandomForestClassifier
from sklearn.externals.joblib import Parallel, delayed
import os

from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC

from core.annotations import read_annotations
from core.constants import data_path, upstream_length
from core.data_reading import read_pheno, read_subset, read_variants
from ml_methods.mtb_predictor import MTBPredictor

path_to_pheno = data_path + 'pheno_mc5_mega/'
# path_to_var = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_no_win_qual_mqm_std3_mqm30_no_highcov_bam_filtered/'
# path_to_var = data_path + 'snps/gatk_before_cortex/raw_variants_ld/'
# path_to_var = data_path + 'snps/freebayes_after_cortex/raw_no_win_qual_mqm_std3_mqm30_no_highcov_ld/'
# path_to_var = data_path + 'snps/skesa_mummer_raw_ld_mum4_q30/'
path_to_var = data_path + 'snps/skesa_minimap_mapq_raw/'
# path_to_var = data_path + 'snps/gatk_before_cortex/raw_variants_no_gvcf_mc10_ld/'
path_to_ids = data_path + 'all_with_pheno.txt'
path_to_subsets = data_path + 'subsets/'

# log_path = '../../res/ml_log_mc3_bam_filtered/'
# log_path = '../../res/ml_svm_mc3_gatk_before/'
# log_path = '../../res/ml_log_mc3_freebayes_after_cortex/'
# log_path = '../../res/ml_log_mc3_skesa_mum4_q30/'
log_path = '../../res/ml_log_mc3_skesa_minimap/'
# log_path = '../../res/ml_log_mc3_gatk_before_no_gvcf_mc10/'

snp_count_threshold = 3

thread_num = 32


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


def read_all_variants(sample_ids):
    sample_to_variants = {}
    tasks = Parallel(n_jobs=thread_num, batch_size=len(sample_ids)//thread_num + 1)(delayed(read_variants)(path_to_var, sample_id) for sample_id in sample_ids)
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

    tasks = Parallel(n_jobs=thread_num)(delayed(split_dataset)(drug, pheno, sample_to_variants, set_name_to_list['Walker_subset'])
                     for drug, pheno in drug_to_pheno.items())

    drug_to_splits = {}
    for drug, snp_train, snp_test, pheno_train, pheno_test in tasks:
        drug_to_splits[drug] = (snp_train, snp_test, pheno_train, pheno_test)

    with open(log_path + 'ml.stats', 'w') as f:
        f.write('snp_count_threshold = %d\n' % snp_count_threshold)

        f.write('LR_l1 C1\n\n')
        clf = LogisticRegression(penalty='l1', dual=False, class_weight='balanced', C=1)
        # clf = LinearSVC(penalty='l2', loss='hinge', class_weight='balanced', max_iter=10000)
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
