import pandas as pd
from ete3 import Tree
from sklearn.externals.joblib import Parallel, delayed

path_to_tree = './data/big_tree_filtered2.newick'
path_to_phenotypes = './data/phenotype_filtered.csv'
out_tree = './data/trees/'
out_pheno = './data/pheno/'


def prune(pheno, drug):
    t = Tree(path_to_tree)
    t.set_outgroup('canetti')
    sample_to_pheno = pheno[['Organism_name', drug]][pheno[drug].isin(('S', 'R'))]
    sample_to_pheno[drug] = sample_to_pheno[drug].map({'S': 0, 'R': 1})
    t.prune(sample_to_pheno['Organism_name'].tolist())
    t.write(format=1, outfile=out_tree + drug + ".nw")
    sample_to_pheno.to_csv(out_pheno + drug + '.pheno', header=None, index=None, sep='\t')
    return 0


def main():
    pheno = pd.read_csv(path_to_phenotypes, sep='\t')
    tasks = Parallel(n_jobs=-1)(delayed(prune)(pheno, drug) for drug in pheno.columns if drug != 'Organism_name')
    c = 0
    for task in tasks:
        c += task
    print(c)


if __name__ == '__main__':
    main()