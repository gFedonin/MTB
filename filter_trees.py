from os import mkdir
from os.path import exists
import os
from ete3 import Tree

from sklearn.externals.joblib import Parallel, delayed

path_to_trees = './data/trees_mc5/'
path_to_pheno = './data/pheno_mc5/'

path_to_filtered_trees = './data/trees_mc5_Walker/'
path_to_filtered_pheno = './data/pheno_mc5_Walker/'

path_to_subset = './data/subsets/Walker_subset.txt'


def filter_pheno(sample_ids, filename):
    with open(path_to_pheno + filename, 'r') as fin:
        with open(path_to_filtered_pheno + filename, 'w')as fout:
            for line in fin.readlines():
                if line.split('\t')[0] in sample_ids:
                    fout.write(line)
    return 0


def filter_tree(sample_ids, filename):
    t = Tree(path_to_trees + filename)
    l = []
    for s in t.get_leaves():
      if s.name in sample_ids:
          l.append(s.name)
    if len(l) > 0:
        t.prune(l)
        t.write(format=1, outfile=path_to_filtered_trees + filename)
    else:
        print(filename)
    return 0


def main():
    if not exists(path_to_filtered_pheno):
        mkdir(path_to_filtered_pheno)
    if not exists(path_to_filtered_trees):
        mkdir(path_to_filtered_trees)
    sample_ids = set()
    with open(path_to_subset, 'r') as f:
        for line in f.readlines():
            sample_ids.add(line.strip())
    for (dirpath, dirnames, filenames) in os.walk(path_to_pheno):
        tasks = Parallel(n_jobs=-1)(delayed(filter_pheno)(sample_ids, filename) for filename in filenames)
        c = 0
        for task in tasks:
            c += task
    for (dirpath, dirnames, filenames) in os.walk(path_to_trees):
        tasks = Parallel(n_jobs=-1)(delayed(filter_tree)(sample_ids, filename) for filename in filenames)
        c = 0
        for task in tasks:
            c += task


if __name__ == '__main__':
    main()