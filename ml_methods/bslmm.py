import gzip
from os.path import exists

import pandas as pd
from sklearn.metrics import roc_auc_score, f1_score
from sklearn.metrics import confusion_matrix
import subprocess

path_to_gemma = '/export/home/fedonin/gemma'

class BSLMM:

    def print_genotype_bimbam(self, pheno, sample_to_variants):
        var_set = set()
        samples = []
        self.samples = samples
        for sample, ph in pheno:
            vars = sample_to_variants.get(sample)
            if vars is None:
                continue
            var_set.update(vars)
            samples.append(sample)
        var_list = list(var_set)
        self.var_list = var_list
        var_list.sort()
        samples.sort()
        if not exists(self.temp_path + self.drug + '.bimbam.gz'):
            vars_to_samples = {}
            for sample in samples:
                vars = sample_to_variants[sample]
                for var in vars:
                    sample_set = vars_to_samples.get(var)
                    if sample_set is None:
                        vars_to_samples[var] = {sample}
                    else:
                        sample_set.add(sample)
            with gzip.open(self.temp_path + self.drug + '.bimbam.gz', 'wt') as f:
                for var in var_list:
                    sample_set = vars_to_samples[var]
                    if self.snp_count_threshold < len(sample_set) < len(samples) - self.snp_count_threshold:
                        f.write(str(var) + ',-,-')
                        for sample in samples:
                            if sample in sample_set:
                                f.write(',1')
                            else:
                                f.write(',0')
                        f.write('\n')

    def __init__(self, drug, pheno, temp_path, sample_to_variants, snp_count_threshold):
        self.temp_path = temp_path
        self.drug = drug
        self.snp_count_threshold = snp_count_threshold
        self.print_genotype_bimbam(pheno, sample_to_variants)

    def train(self, mut_train, pheno_train, split_id):
        if not exists(self.temp_path + self.drug + '_' + str(split_id) + '.pheno'):
            with open(self.temp_path + self.drug + '_' + str(split_id) + '.pheno', 'w') as f:
                for sample in self.samples:
                    ph = pheno_train.get(sample)
                    if ph is None:
                        f.write('NA\n')
                    else:
                        f.write(str(ph) + '\n')

        if not exists(self.drug + '_' + str(split_id) + '.param.txt'):
            subprocess.call([path_to_gemma,
                                   "-g", self.temp_path + self.drug + ".bimbam.gz",
                                   "-p", self.temp_path + self.drug + '_' + str(split_id) + ".pheno",
                                   "-bslmm", "1", "-o", self.drug + '_' + str(split_id)],
                                  cwd=self.temp_path)

    def compute_stat(self, snp_test, pheno_test, split_id):
        if not exists(self.temp_path + "output/" + self.drug + '_' + str(split_id) + ".prdt.txt"):
            subprocess.call([path_to_gemma,
                                   "-g", self.temp_path + self.drug + ".bimbam.gz",
                                   "-p", self.temp_path + self.drug + '_' + str(split_id) + ".pheno",
                                   "-epm", "output/" + self.drug + '_' + str(split_id) + ".param.txt",
                                   "-emu", "output/" + self.drug + '_' + str(split_id) + ".log.txt",
                                   "-predict", "1", "-o", self.drug + '_' + str(split_id)],
                                  cwd=self.temp_path)
        sample_to_pred = {}
        f = open(self.temp_path + "output/" + self.drug + '_' + str(split_id) + ".prdt.txt")
        for sample in self.samples:
            s = f.readline().strip()
            if s != 'NA':
                if float(s) > 0.5:
                    sample_to_pred[sample] = 1
                else:
                    sample_to_pred[sample] = 0
        y_pred = []
        y_test = []
        for sample, pheno in pheno_test.items():
            y_test.append(pheno_test[sample])
            y_pred.append(sample_to_pred[sample])
        tn, fp, fn, tp = confusion_matrix(y_test, y_pred, labels=[0,1]).ravel()
        if tp + fp > 0:
            ppv = tp / (tp + fp)
        else:
            ppv = 0
        if tn + fn > 0:
            npv = tn / (tn + fn)
        else:
            npv = 0
        if tp + fn > 0:
            sensitivity = tp / (tp + fn)
        else:
            sensitivity = 0
        if tn + fp > 0:
            specificity = tn / (tn + fp)
        else:
            specificity = 0
        try:
            f1 = f1_score(y_test, y_pred)
        except:
            f1 = 0
        try:
            auc = roc_auc_score(y_test, y_pred)
        except:
            auc = 0
        return ppv, npv, sensitivity, specificity, f1, auc