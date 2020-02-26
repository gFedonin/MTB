from os.path import exists

from sklearn.ensemble import RandomForestClassifier
from sklearn.externals.joblib import Parallel, delayed
import os

from sklearn.linear_model import LogisticRegression
from sklearn.utils import shuffle

from core.constants import data_path
from core.data_reading import read_pheno, read_subset, read_variants
from ml_methods.bslmm import BSLMM
from ml_methods.mtb_predictor import MTBPredictor

# path_to_pheno = data_path + 'pheno_mc5_mega/'
path_to_pheno = data_path + 'combined_pheno/'
# path_to_var = data_path + 'snps/annotated_with_DR_with_indel_with_pheno_and_snp_no_win_mqm_std3_mqm30_long_del_pg_filter_samples_first/'
# path_to_var = data_path + 'snps/annotated_pg_NWds10_no_win_qual_mqm_std3_mqm30_filter_samples_first/'
# path_to_var = data_path + 'snps/gatk_before_cortex/annotated_long_del_pg_NWds10/'
# path_to_var = data_path + 'snps/skesa_mum4_annotated_long_del_pg_NWds10/'
# path_to_var = data_path + 'snps/pilon/annotated_pg_NWds10_filtered_test/'
# path_to_var = data_path + 'snps/pilon/annotated_pg_NWds10_filtered_ext/'
# path_to_var = data_path + 'snps/gatk_before_cortex/annotated_pg_NWds10_mq40_keep_complex_filtered_test/'
# path_to_var = data_path + 'snps/gatk_before_cortex/annotated_pg_NWds10_mq40_keep_complex_filtered_ext/'
# path_to_var = data_path + 'snps/gatk_and_pilon/unification_ext/'
# path_to_var = data_path + 'snps/gatk_and_pilon/intersection_ext/'
path_to_var = data_path + 'snps/annotated_fixed_no_rep_long_del_pg_NWds10_combined_extended/'
# path_to_ids = data_path + 'snps/intersect.list'#'snps/raw_with_DR_with_indel_with_pheno_and_snp_no_win_m3_filter_samples_first/samples_filtered.list'#'dr_covered_with_pheno_and_snp_new.txt'
# path_to_ids = data_path + 'all_with_pheno.txt'
# path_to_ids = data_path + 'snps/gatk_pilon_old_intersection.list'
path_to_ids = data_path + 'snps/combined_raw_variants_mq40_keep_complex_std_names_filtered/samples_filtered.list'
path_to_subsets = data_path + 'subsets/'

# log_path = '../../res/ml_log_mc1-noncds_long_del_pg_ext_NWds10_no_win_qual_mqm_std3_long_del_pg_filter_samples_first/'
# log_path = '../../res/ml_log_mc1-noncds_long_del_pg_no_win_qual_mqm_std3_filter_samples_first/'
# log_path = '../../res/ml_log_mc3_gatk_before_annotated_long_del_pg_NWds10/'
# log_path = '../../res/ml_log_mc3_skesa_mum4_annotated_long_del_pg_NWds10/'
# log_path = '../../res/ml_log_mc1_pilon_annotated_long_del_pg_NWds10_filtered/'
# log_path = '../../res/ml_log_mc1_gatk_annotated_long_del_pg_NWds10_mq40_keep_complex_filtered/'
# log_path = '../../res/ml_log_mc1_gatk_pilon_unification_cv_rep/'
# log_path = '../../res/ml_log_mc1_gatk_pilon_intersection_cv_rep/'
# log_path = '../../res/lmm_gatk_pilon_unification_cv_rep/'
# log_path = '../../res/lmm_gatk_pilon_intersection_cv_rep/'
# log_path = '../../res/ml_log_mc1_combined/'
log_path = '../../res/lmm_mc1_combined/'
# log_path = '../../res/lmm_gatk_cv_rep/'
# log_path = '../../res/lmm_pilon_cv_rep/'

use_bslmm = True

use_extended_features = True
merge_all_mut_in_pos = True
merge_all_mut_in_gene = True
filter_non_cds = True
drop_pseudogenes = True
drop_hypothetical = False
keep_proteomic_validated_only = False
keep_only_not_hypotetical_or_proteomic_validated = False
upstream_indel_breakes_gene = True
snp_count_threshold = 1

use_pairs = False
jaccard_index_threshold = 0.9

thread_num = 144

filter_by_var_list = False
path_to_var_list = path_to_var + 'dr_genes_var_list.csv'


def split_dataset(pheno, sample_to_mut, cv):
    splits = []
    all_samples = [sample_id for sample_id, val in pheno if sample_id in sample_to_mut]
    shuffle(all_samples, random_state=0)
    for i in range(cv):
        train_i = []
        test_i = []
        for j in range(len(all_samples)):
            sample_id = all_samples[j]
            if j%cv == i:
                test_i.append(sample_id)
            else:
                train_i.append(sample_id)
        splits.append((train_i, test_i))
    return splits


def drug_stat_multi(predictor, snp_train, snp_test, pheno_train, pheno_test, split_id):

    predictor.train(snp_train, pheno_train, split_id)

    ppv, npv, sensitivity, specificity, f1, auc = predictor.compute_stat(snp_test, pheno_test, split_id)

    return predictor.drug, ppv, npv, sensitivity, specificity, f1, auc


def gen_predictors(sample_to_variants, drug_to_pheno, clf, n_jobs=-1):
    drug_to_predictor = {}
    if n_jobs == 1:
        for drug, pheno in drug_to_pheno.items():
            if use_bslmm:
                predictor = BSLMM(drug, pheno, log_path, sample_to_variants, snp_count_threshold)
            else:
                predictor = MTBPredictor(drug, clf, snp_count_threshold, log_path + drug + '.log')
            drug_to_predictor[drug] = predictor
    else:
        if use_bslmm:
            predictors = Parallel(n_jobs=n_jobs)(delayed(BSLMM)(drug, pheno, log_path, sample_to_variants, snp_count_threshold)
                                                 for drug, pheno in drug_to_pheno.items())
        else:
            predictors = Parallel(n_jobs=n_jobs)(delayed(MTBPredictor)(drug, clf, snp_count_threshold,
                                                                       log_path + drug + '.log')
                                                 for drug, pheno in drug_to_pheno.items())
        it = iter(predictors)
        for drug, pheno in drug_to_pheno.items():
            drug_to_predictor[drug] = next(it)
    return drug_to_predictor


def print_stat_fixed_multi(f, drug_to_splits, sample_to_variants, drug_to_pheno, clf, n_jobs=-1):
    drug_to_predictors = gen_predictors(sample_to_variants, drug_to_pheno, clf, n_jobs)
    f.write('Drug\tSample_num\tPPV\tNPV\tSensitivity\tSpecificity\tF1\tAUC\n')
    if n_jobs != 1:
        tasks = []
        for drug, splits in drug_to_splits.items():
            predictor = drug_to_predictors[drug]
            pheno = drug_to_pheno[predictor.drug]
            pheno_dict = {sample_id: ph for sample_id, ph in pheno}
            split_id = 0
            for train, test in splits:
                snp_train = {sample_id: sample_to_variants[sample_id] for sample_id in train}
                snp_test = {sample_id: sample_to_variants[sample_id] for sample_id in test}
                pheno_train = {sample_id: pheno_dict[sample_id] for sample_id in train}
                pheno_test = {sample_id: pheno_dict[sample_id] for sample_id in test}
                tasks.append((predictor, snp_train, snp_test, pheno_train, pheno_test, split_id))
                split_id += 1
        stats = Parallel(n_jobs=n_jobs)(delayed(drug_stat_multi)(predictor, snp_train, snp_test, pheno_train,
                                                                 pheno_test, split_id)
                                        for predictor, snp_train, snp_test, pheno_train, pheno_test, split_id in tasks)
        drug_to_stat = {}
        for drug, pheno in drug_to_pheno.items():
            stat = {}
            drug_to_stat[drug] = stat
            stat['av_ppv'] = 0
            stat['av_npv'] = 0
            stat['av_sens'] = 0
            stat['av_spec'] = 0
            stat['av_f1'] = 0
            stat['av_auc'] = 0
        for drug, ppv, npv, sensitivity, specificity, f1, auc in stats:
            stat = drug_to_stat[drug]
            stat['av_ppv'] += ppv
            stat['av_npv'] += npv
            stat['av_sens'] += sensitivity
            stat['av_spec'] += specificity
            stat['av_f1'] += f1
            stat['av_auc'] += auc
        for drug, stat in drug_to_stat.items():
            splits = drug_to_splits[drug]
            pheno = drug_to_pheno[drug]
            split_num = len(splits)
            for k,v in stat.items():
                stat[k] /= split_num
            f.write(drug)
            f.write('\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n' % (len(pheno), stat['av_ppv'], stat['av_npv'],
                    stat['av_sens'], stat['av_spec'], stat['av_f1'], stat['av_auc']))
    else:
        for drug, splits in drug_to_splits.items():
            predictor = drug_to_predictors[drug]
            pheno = drug_to_pheno[drug]
            split_num = len(splits)
            av_ppv = 0
            av_npv = 0
            av_sens = 0
            av_spec = 0
            av_f1 = 0
            av_auc = 0
            split_id = 0
            for train, test in splits:
                drug, ppv, npv, sensitivity, specificity, f1, auc = \
                    drug_stat_multi(predictor, train, test, sample_to_variants, pheno, split_id)
                av_ppv += ppv
                av_npv += npv
                av_sens += sensitivity
                av_spec += specificity
                av_f1 += f1
                av_auc += auc
                split_id += 1

            f.write(drug)
            f.write('\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n' % (len(pheno), av_ppv / split_num,
                av_npv / split_num, av_sens / split_num, av_spec / split_num, av_f1 / split_num, av_auc / split_num))


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
    tasks = Parallel(n_jobs=thread_num, batch_size=len(sample_ids)//thread_num + 1)(delayed(read_variants)(path_to_var, sample_id)
                                    for sample_id in sample_ids if exists(path_to_var + sample_id + '.variants'))
    for name, snp_list in tasks:
        sample_to_variants[name] = snp_list
    print('samples\' variants reading done')
    return sample_to_variants


def shrink_features(sample_to_variants):
    sample_to_shrinked = {}
    var_set = set()
    for sample, vars in sample_to_variants.items():
        var_set.update(vars)
    var_list = list(var_set)
    var_list.sort()
    var_to_id = {var: i for i, var in enumerate(var_list)}
    for sample, vars in sample_to_variants.items():
        sample_to_shrinked[sample] = [var_to_id[var] for var in vars]
    return sample_to_shrinked, var_list


def main():

    if not exists(log_path):
        os.makedirs(log_path)

    drug_to_pheno = read_all_pheno()

    with open(path_to_ids, 'r') as f:
        sample_ids = [name.strip() for name in f.readlines()]
    print(str(len(sample_ids)) + ' samples')

    sample_to_variants = read_all_variants(sample_ids)
    sample_to_shrinked, var_list = shrink_features(sample_to_variants)

    drug_to_splits = {}
    for drug, pheno in drug_to_pheno.items():
        drug_to_splits[drug] = split_dataset(pheno, sample_to_shrinked, 10)
    print('splitting done')

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
        clf = LogisticRegression(penalty='l1', dual=False, class_weight='balanced', C=1, solver='liblinear')
        # f.write('NN k = 1\n\n')
        # clf = KNeighborsClassifier(n_neighbors=1, p=1, n_jobs=-1)
        # f.write('RF n=10000\n\n')
        # clf = RandomForestClassifier(n_estimators=10000, n_jobs=-1, class_weight='balanced')
        # f.write('MI GFFS mi_k=10000fs_k=1000 LR_l2 C1\n\n')
        # clf = MIGFFSLR(penalty='l2', dual=False, class_weight='balanced', C=1, mi_k=10000, fs_k=1000, n_jobs=-1)
        print_stat_fixed_multi(f, drug_to_splits, sample_to_shrinked, drug_to_pheno, clf, n_jobs=-1)
        f.write('\n')


if __name__ == '__main__':
    main()
