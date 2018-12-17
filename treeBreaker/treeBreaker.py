import os

from os.path import exists
from sklearn.externals.joblib import Parallel, delayed

path_to_tries = '../data/trees_mc5_Walker/'
path_to_pheno = '../data/pheno_mc5_Walker/'
out_path = '../data/tree_breaker_mc5_Walker/'


def run_tree_braker(drug):
    os.system('./treeBreaker ' + path_to_tries + drug + '.nw ' + path_to_pheno + drug + '.pheno ' + out_path + drug + '.out')
    return 0

if not exists(out_path):
    os.mkdir(out_path)

drugs = []
for (dirpath, dirnames, filenames) in os.walk(path_to_tries):
    for filename in filenames:
        drugs.append(filename.split('.')[0])

tasks = Parallel(n_jobs=-1)(delayed(run_tree_braker)(drug) for drug in drugs)
c = 0
for task in tasks:
    c += task
print(str(c))