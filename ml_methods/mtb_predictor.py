from bisect import bisect_left
from random import shuffle

import scipy.stats as stat

from scipy.sparse import lil_matrix
from sklearn import clone
from sklearn.metrics import confusion_matrix, roc_auc_score, f1_score
from sklearn.externals.joblib import Parallel, delayed

from core.annotations import read_annotations, CDSType
from core.constants import upstream_length
from data_processing.print_mutations_extended import process_variants
import numpy as np


class MTBPredictor:

    max_p_value = 0.01
    iterNum = 100
    compute_pvalues = False

    def __init__(self, drug, clf, snp_count_threshold, log_path):
        self.drug = drug
        self.snp_count_threshold = snp_count_threshold
        self.clf = clf
        self.log_path = log_path
        self.p_values = []

    def snps_counts(self, sample_to_snps):
        snp_to_count = {}
        for l in sample_to_snps.values():
            for snp in l:
                c = snp_to_count.get(snp)
                if c is not None:
                    snp_to_count[snp] = c + 1
                else:
                    snp_to_count[snp] = 1
        return snp_to_count

    def filter_features(self, mut_train):
        all_snps_set = set()
        sample_num = len(mut_train)
        snp_to_count = self.snps_counts(mut_train)
        mut_train_filtered = {}
        for sample_id, snp_list in mut_train.items():
            res = []
            for snp in snp_list:
                s = snp.split('\t')
                if self.snp_count_threshold < snp_to_count[snp] < sample_num - self.snp_count_threshold:
                    res.append(snp)
                    all_snps_set.add(snp)
            mut_train_filtered[sample_id] = res
        full_snp_list = list(all_snps_set)
        return mut_train_filtered, full_snp_list

    def create_sparse_matrix(self, mut_train, pheno_train):
        mut_train, self.full_snp_list = self.filter_features(mut_train)

        sample_ids = list(mut_train.keys())
        sample_num = len(sample_ids)
        self.snp_num = len(self.full_snp_list)

        self.feature_to_index = {}
        for i in range(self.snp_num):
            self.feature_to_index[self.full_snp_list[i]] = i

        X_train = lil_matrix((sample_num, self.snp_num), dtype=int)
        y_train = []
        for i in range(len(sample_ids)):
            for f in mut_train[sample_ids[i]]:
                X_train[i, self.feature_to_index[f]] = 1
            y_train.append(pheno_train[sample_ids[i]])
        return X_train, y_train

    def pvalues(self, X_train, y_train):
        random_coefs = []
        for i in range(self.iterNum):
            clf = clone(self.clf)
            y_rand = y_train.copy()
            shuffle(y_rand)
            clf.fit(X_train, y_rand)
            abs_coefs = [abs(c) for c in clf.coef_[0]]
            random_coefs.append(max(abs_coefs))
        random_coefs.sort()
        for i in range(self.snp_num):
            j = bisect_left(random_coefs, abs(self.clf.coef_[0][i]))
            self.p_values.append((self.iterNum - j)/self.iterNum)

    def train(self, mut_train, pheno_train):
        X_train, y_train = self.create_sparse_matrix(mut_train, pheno_train)
        self.clf.fit(X_train, y_train)
        if self.compute_pvalues:
            self.pvalues(X_train, y_train)
        # denom = (2.0*(1.0+np.cosh(self.clf.decision_function(X_train))))
        # F_ij = np.dot((X_train/denom[:,None]).T, X_train) ## Fisher Information Matrix
        # Cramer_Rao = np.linalg.inv(F_ij) ## Inverse Information Matrix
        # sigma_estimates = np.array([np.sqrt(Cramer_Rao[i,i]) for i in range(Cramer_Rao.shape[0])]) # sigma for each coefficient
        # self.z_scores = self.clf.coef_[0]/sigma_estimates # z-score for eaach model coefficient
        # self.p_values = [stat.norm.sf(abs(x))*2*self.snp_num for x in self.z_scores] ### two tailed test for p-values

    def compute_stat(self, mut_test, pheno_test):
        sample_ids = list(mut_test.keys())
        sample_num = len(sample_ids)
        X_test = lil_matrix((sample_num, self.snp_num), dtype=int)
        y_test = []
        for i in range(len(sample_ids)):
            for f in mut_test[sample_ids[i]]:
                index = self.feature_to_index.get(f)
                if index is not None:
                    X_test[i, index] = 1
            y_test.append(pheno_test[sample_ids[i]])
        y_pred = self.clf.predict(X_test)
        cm = confusion_matrix(y_test, y_pred)
        tp_score = cm[1, 1]
        tn_score = cm[0, 0]
        fp_score = cm[0, 1]
        fn_score = cm[1, 0]
        if tp_score + fp_score > 0:
            ppv = tp_score / (tp_score + fp_score)
        else:
            ppv = 0
        if tn_score + fn_score > 0:
            npv = tn_score / (tn_score + fn_score)
        else:
            npv = 0
        if tp_score + fn_score > 0:
            sensitivity = tp_score / (tp_score + fn_score)
        else:
            sensitivity = 0
        if tn_score + fp_score > 0:
            specificity = tn_score / (tn_score + fp_score)
        else:
            specificity = 0
        return ppv, npv, sensitivity, specificity, f1_score(y_test, y_pred), roc_auc_score(y_test, y_pred)

    def print_selected_features(self):
        with open(self.log_path, 'w') as f:
            if self.compute_pvalues:
                f.write('feature\tcoef\tp_value\n')
            else:
                f.write('feature\tcoef\n')
            for i in range(self.snp_num):
                # if self.p_values[i] < self.max_p_value:
                f.write(self.full_snp_list[i])
                f.write('\t')
                f.write(str(self.clf.coef_[0][i]))
                if self.compute_pvalues:
                    f.write('\t')
                    f.write(str(self.p_values[i]))
                f.write('\n')


class MTBPredictorExtended(MTBPredictor):

    def __init__(self, drug, clf, snp_count_threshold, merge_all_mut_in_pos, merge_all_mut_in_gene,
                 upstream_indel_breakes_gene, log_path):
        super().__init__(drug, clf, snp_count_threshold, log_path)
        self.merge_all_mut_in_pos = merge_all_mut_in_pos
        self.merge_all_mut_in_gene = merge_all_mut_in_gene
        self.upstream_indel_breakes_gene = upstream_indel_breakes_gene

    def create_sparse_matrix(self, mut_train, pheno_train):
        cds_list = read_annotations(upstream_length, filter_by_gene_len=False)
        self.name_to_type = {}
        for cds in cds_list:
            if cds.type != CDSType.upstream:
                self.name_to_type[cds.name] = cds.type
        for sample_id, snp_list in mut_train.items():
            sample_id, mut_list = process_variants(sample_id, snp_list, self.name_to_type, self.merge_all_mut_in_pos,
                                      self.merge_all_mut_in_gene, False,
                                      self.upstream_indel_breakes_gene)
            mut_train[sample_id] = mut_list
        return super().create_sparse_matrix(mut_train, pheno_train)

    def compute_stat(self, mut_test, pheno_test):
        for sample_id, snp_list in mut_test.items():
            sample_id, mut_list = process_variants(sample_id, snp_list, self.name_to_type, self.merge_all_mut_in_pos,
                                      self.merge_all_mut_in_gene, False, self.upstream_indel_breakes_gene)
            mut_test[sample_id] = mut_list
        return super().compute_stat(mut_test, pheno_test)


class MTBPredictorPairs(MTBPredictorExtended):

    def __init__(self, drug, clf, snp_count_threshold, merge_all_mut_in_pos, merge_all_mut_in_gene,
                 upstream_indel_breakes_gene, jaccard_index_threshold,
                 n_jobs, log_path):
        super().__init__(drug, clf, snp_count_threshold, merge_all_mut_in_pos, merge_all_mut_in_gene,
                         upstream_indel_breakes_gene, log_path)
        self.n_jobs = n_jobs
        self.jaccard_index_threshold = jaccard_index_threshold

    def transpose_snp_matrix(self, sample_to_snps, pheno):
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

    def gen_pairs_for_snp(self, snp_pivot_list, sample_to_snps, snp_to_samples, sample_num, pheno, snp_to_resistant_count):
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
                if self.snp_count_threshold < c < sample_num - self.snp_count_threshold:
                    snp_c = len(snp_to_samples[pair[1]])
                    if pivot_c - c > self.snp_count_threshold and snp_c - c > self.snp_count_threshold:
                        if c / (pivot_c + snp_c - c) < self.jaccard_index_threshold:
                            snp_r_c = snp_to_resistant_count[pair[1]]
                            pair_r_c = snp_pair_to_r_count[pair]
                            if pair_r_c * pivot_c * snp_c != pivot_r_count * snp_r_c * c:
                                res.append(pair)
        return res

    def gen_all_pairs(self, sample_to_snps, all_snp_list, pheno):
        snp_to_samples, snp_to_resistant_count = self.transpose_snp_matrix(sample_to_snps, pheno)
        sample_num = len(sample_to_snps)
        # thread_num = os.cpu_count()
        snp_lists = {}
        for i in range(self.n_jobs):
            snp_lists[i] = []
        for i in range(len(all_snp_list)):
            snp_lists[i % self.n_jobs].append(all_snp_list[i])
        tasks = Parallel(n_jobs=self.n_jobs)(
            delayed(self.gen_pairs_for_snp)(snp_list, sample_to_snps, snp_to_samples, sample_num,
                                       pheno, snp_to_resistant_count)
            for snp_list in snp_lists.values())
        pairs = set()
        for task in tasks:
            for pair in task:
                pairs.add(pair[0] + '_' + pair[1])
        return pairs

    def add_pairs_to_all_snp_lists(self, sample_to_snps, filtered_pairs):
        for sample_id, snp_list in sample_to_snps.items():
            snp_num = len(snp_list)
            for i in range(snp_num):
                for j in range(i + 1, snp_num):
                    pair = snp_list[i] + '_' + snp_list[j]
                    if pair in filtered_pairs:
                        snp_list.append(pair)
        return sample_to_snps

    def update_snp_lists(self, sample_to_snps, pairs):
        sample_ids = list(sample_to_snps.keys())
        sample_to_snps_array = {}
        for i in range(self.n_jobs):
            sample_to_snps_array[i] = {}
        for i in range(len(sample_ids)):
            sample_id = sample_ids[i]
            sample_to_snps_array[i % self.n_jobs][sample_id] = sample_to_snps[sample_id]

        tasks1 = Parallel(n_jobs=self.n_jobs)(delayed(self.add_pairs_to_all_snp_lists)(sample_to_snps_i, pairs)
                                             for sample_to_snps_i in sample_to_snps_array.values())

        for task in tasks1:
            for sample_id, snp_list in task.items():
                sample_to_snps[sample_id] = snp_list

    def create_sparse_matrix(self, mut_train, pheno_train):

        cds_list = read_annotations(upstream_length, filter_by_gene_len=False)
        self.name_to_type = {}
        for cds in cds_list:
            if cds.type != CDSType.upstream:
                self.name_to_type[cds.name] = cds.type
        for sample_id, snp_list in mut_train.items():
            sample_id, mut_list = process_variants(sample_id, snp_list, self.name_to_type, self.merge_all_mut_in_pos,
                                      self.merge_all_mut_in_gene, False,
                                      self.upstream_indel_breakes_gene)
            mut_train[sample_id] = mut_list
        mut_train, full_snp_list = self.filter_features(mut_train)
        self.pairs = self.gen_all_pairs(mut_train, full_snp_list, pheno_train)
        with open(self.log_path, 'a') as f:
            f.write(self.drug + ' ' + str(len(self.pairs)) + ' pairs generated\n')

            if len(self.pairs) > 0:
                for pair in self.pairs:
                    full_snp_list.append(pair)
                f.write(self.drug + ': ' + str(len(full_snp_list)) + ' variants + pairs after filtering\n')
                self.update_snp_lists(mut_train, self.pairs)

        sample_ids = list(mut_train.keys())
        sample_num = len(sample_ids)
        self.snp_num = len(full_snp_list)

        self.feature_to_index = {}
        for i in range(self.snp_num):
            self.feature_to_index[full_snp_list[i]] = i

        X_train = lil_matrix((sample_num, self.snp_num), dtype=int)
        y_train = []
        for i in range(len(sample_ids)):
            for f in mut_train[sample_ids[i]]:
                X_train[i, self.feature_to_index[f]] = 1
            y_train.append(pheno_train[sample_ids[i]])
        return X_train, y_train

    def compute_stat(self, mut_test, pheno_test):
        for sample_id, snp_list in mut_test.items():
            sample_id, mut_list = process_variants(sample_id, snp_list, self.name_to_type, self.merge_all_mut_in_pos,
                                      self.merge_all_mut_in_gene, False,
                                      self.upstream_indel_breakes_gene)
            mut_test[sample_id] = mut_list
        if len(self.pairs) > 0:
            self.update_snp_lists(mut_test, self.pairs)
        return super().compute_stat(mut_test, pheno_test)