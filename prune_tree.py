import pandas as pd
from ete3 import Tree
from os.path import exists

from os import makedirs
from sklearn.externals.joblib import Parallel, delayed

path_to_tree = './data/tree_dr_covered_with_pheno_and_snp_mc5_rooted.nw'
path_to_phenotypes = './data/dr_covered_with_pheno_and_snp.csv'
out_tree = './data/trees_mc5/'
out_pheno = './data/pheno_mc5/'


def prune(pheno, drug):
    t = Tree(path_to_tree)
    # t.set_outgroup('canetti')
    sample_to_pheno = pheno[['Organism_name', drug]][pheno[drug].isin(('S', 'R'))]
    sample_to_pheno[drug] = sample_to_pheno[drug].map({'S': 0, 'R': 1})
    t.prune(sample_to_pheno['Organism_name'].tolist())
    t.write(format=1, outfile=out_tree + drug + ".nw")
    sample_to_pheno.to_csv(out_pheno + drug + '.pheno', header=None, index=None, sep='\t')
    return 0


def main():
    if not exists(out_tree):
        makedirs(out_tree)
    if not exists(out_pheno):
        makedirs(out_pheno)
    pheno = pd.read_csv(path_to_phenotypes, sep='\t')
    tasks = Parallel(n_jobs=-1)(delayed(prune)(pheno, drug) for drug in pheno.columns if drug != 'Organism_name')
    c = 0
    for task in tasks:
        c += task
    print(c)


if __name__ == '__main__':
    main()