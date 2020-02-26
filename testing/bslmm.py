import pandas as pd
import numpy as np

from core.constants import data_path

PATH_TO_GENOTYPES = data_path + "snps/gatk_and_pilon/intersection_ext/"
PATH_TO_PHENOTYPES = "pheno_mc5_mega"
path_to_ids = data_path + 'snps/gatk_pilon_old_intersection.list'
import os

dir_ = PATH_TO_GENOTYPES #path to genotype files folder
directory = os.fsencode(dir_)

# obtain set of all features
Set_of_features = set()
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    with open(dir_ + filename) as f:
        s = f.read().splitlines()
        Set_of_features.update(s)
List_of_features = list(Set_of_features)
N = 0 # total number of samples
for file in os.listdir(directory):
    N += 1

#generate naive features to X numpy array
D = {i:j for j, i in enumerate(List_of_features)}
X = np.zeros((N, len(List_of_features)))
k = 0
names = [] #names of samples
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    with open(dir_ + filename) as f:
        features = f.read().splitlines()
        for i in features:
            X[k, D[i]] = 1
        names.append(filename[:-9])
    k += 1

# extract mutations that occured more than 5 times
# and write them to df with columns corresponding to Mean Genotype Data file format
from itertools import compress
s = X.sum(axis = 0)
df = pd.DataFrame(X[:,s  > 5], index = names,
                  columns = list(compress(List_of_features, s > 5)))
cols = [i.replace("\t", "") + ",A,T" for i in df.columns]
df.columns = cols

import warnings
warnings.filterwarnings('ignore')
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from sklearn.metrics import roc_auc_score, make_scorer, accuracy_score
from sklearn.metrics import confusion_matrix
roc = make_scorer(roc_auc_score)
acc = make_scorer(accuracy_score)
conf = make_scorer(confusion_matrix)
from sklearn.model_selection import StratifiedShuffleSplit
import os
import subprocess
dir_ = PATH_TO_PHENOTYPES #path to phenotype files folder

f = open("summary.txt", "w") # file with models metrics

for filename in os.listdir(dir_):
    print(filename[:-6])
    antibio = filename[:-6]
  #  if (filename[:-6] = "Isoniazid"): #You can chose which antibiotics to explore in current job
    #    continue
    f.write(filename[:-6] + "\n")
    resistance = pd.read_csv(dir_ + "/" + filename, delimiter="[ \t]", names = ["name", "resistance"])
    # calculate metrics of logistic regrission on the data
    X = df.loc[resistance.name].values
    y = resistance.resistance.values
    model = LogisticRegression(penalty="l1", max_iter=100, C = 1, solver = "liblinear")
    Split = StratifiedShuffleSplit(n_splits = 1, test_size=0.2, train_size=0.8, random_state=1)
    print("Logistic Regression")
    for train_index, test_index in Split.split(X, y):
        model.fit(X[train_index], y[train_index]) 
        logreg_auc = roc(model, X[test_index], y[test_index])
        logreg_conf = conf(model, X[test_index], y[test_index])
        print("AUC:", roc(model, X[test_index], y[test_index]))
        print("Confusion matrix:")
        print(conf(model, X[test_index], y[test_index]))
    # write Mean Genotype file
    df.loc[resistance.name].T.to_csv("run_bslmm/" + antibio +"01.mgt", header = False, index=True, sep = ",")
    subprocess.call(["sed", "-i", 's/"//g', "run_bslmm/" + antibio +"01.mgt"])
    subprocess.call(["sed", "-i", 's/,/, /g', "run_bslmm/" + antibio +"01.mgt"])
    # write Phenotype file
    y = y.astype("float")
    np.savetxt("run_bslmm/" + antibio +"01.truepheno", y, delimiter = "\t", fmt = '%1.1f')
    y[test_index] = np.nan
    np.savetxt("run_bslmm/" + antibio +"01.pheno", y, delimiter = "\t", fmt = '%1.1f')
    subprocess.call(["sed", "-i", 's/nan/NA/g', "run_bslmm/" + antibio + "01.pheno"])
    #run gemma (fitting model)
    subprocess.check_call(["../gemma",
                           "-g", antibio + "01.mgt", 
                           "-p", antibio + "01.pheno", 
                           "-bslmm", "1",
                           "-o", antibio + "01"],
                          cwd="run_bslmm/")
    #predict phenotypes
    subprocess.check_call(["../gemma",
                           "-g", antibio + "01.mgt",
                           "-p", antibio + "01.pheno",
                           "-epm", "output/" + antibio + "01" + ".param.txt",
                           "-emu", "output/" + antibio + "01" + ".log.txt",
                           "-predict", "1",
                           "-o", antibio + "01"],
                          cwd="run_bslmm/")
    #obtain predicted phenotypes
    true = pd.read_csv("run_bslmm/" + antibio + "01.truepheno", sep = " ", names = ["label"])
    pred = pd.read_csv("run_bslmm/output/" + antibio + "01" +".prdt.txt", sep = " ", names = ["prob"])
    y_true = true.label[pred.prob.notna()]
    y_proba = pred.prob[pred.prob.notna()]
    bslmm_auc = roc_auc_score(y_true = y_true, y_score = y_proba)
    #simple grid search for bslmm threshold
    xx = np.linspace(0., 1., 100)
    thr = -1
    for i in xx:
        labels = y_proba > i
        if thr < accuracy_score(y_true, labels):
            pred_labels = labels
            thr = accuracy_score(y_true, labels)
            bslmm_conf= confusion_matrix(y_true = y_true, y_pred = labels)
    print("BSLMM")
    print("AUC:", bslmm_auc)
    print("Confusion matrix:")
    print(bslmm_conf)
    # write information obtained to summary file
    f.write("Logistic Regression\n")
    f.write("AUC score:" + "\n" + str(logreg_auc) + "\n")
    f.write("Confusion matrix:\n")
    f.write(str(logreg_conf[0, 0]) + " " + str(logreg_conf[0, 1]) + "\n")
    f.write(str(logreg_conf[1, 0]) + " " + str(logreg_conf[1, 1]) + "\n")
    f.write("BSLMM\n")
    f.write("AUC score:" + "\n" + str(bslmm_auc) + "\n")
    f.write("Confusion matrix:\n")
    f.write(str(bslmm_conf[0, 0]) + " " + str(bslmm_conf[0, 1]) + "\n")
    f.write(str(bslmm_conf[1, 0]) + " " + str(bslmm_conf[1, 1]) + "\n")
    f.write("\n")
f.close()
