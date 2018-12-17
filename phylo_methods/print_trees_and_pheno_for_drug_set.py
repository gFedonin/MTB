import pandas as pd
from ete3 import Tree
from os.path import exists

from os import makedirs
from sklearn.externals.joblib import Parallel, delayed

from src.core.constants import data_path

path_to_tree = data_path + 'tree_with_pheno_and_snp_mc5_mega_rooted.nw'
path_to_phenotypes = data_path + 'dr_covered_with_pheno_and_snp.csv'
out_tree = data_path + 'trees_mc5_mega_mix/'
out_pheno = data_path + 'pheno_mc5_mega_mix/'

# drug_names = ('Isoniazid', 'Rifampicin', 'Ethambutol', 'Pyrazinamide', 'Streptomycin')
drug_names = ('Isoniazid', 'Rifampicin', 'Ethambutol', 'Pyrazinamide', 'Streptomycin', 'Moxifloxacin', 'Ofloxacin',
              'Amikacin', 'Capreomycin', 'Kanamycin')


def prune(sample_ids, drug):
    t = Tree(path_to_tree)
    # t.set_outgroup('canetti')
    t.prune(sample_ids)
    t.write(format=1, outfile=out_tree + drug + ".nw")
    return 0


def process_pheno():
    all_pheno = pd.read_csv(path_to_phenotypes, sep='\t')
    drug_samples_pairs = []
    for drug in all_pheno.columns[1:]:
        all_pheno[drug] = all_pheno[drug].map({'S': '0', 'R': '1'})

    all_pheno_no_na = all_pheno.dropna(subset=drug_names)
    for drug in drug_names:
        all_pheno_no_na[['Organism_name', drug]].to_csv(out_pheno + drug + '.pheno', sep='\t', index=False, header=False)
        drug_samples_pairs.append((drug, all_pheno_no_na['Organism_name'].tolist()))

    return drug_samples_pairs


def main():
    if not exists(out_tree):
        makedirs(out_tree)
    if not exists(out_pheno):
        makedirs(out_pheno)
    drug_samples_pairs = process_pheno()
    tasks = Parallel(n_jobs=-1)(delayed(prune)(sample_ids, drug) for drug, sample_ids in drug_samples_pairs)
    c = 0
    for task in tasks:
        c += task
    print(c)


if __name__ == '__main__':
    main()
