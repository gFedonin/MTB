import pandas as pd
from ete3 import Tree
from os.path import exists
from os import makedirs

from core.constants import data_path

path_to_tree = data_path + 'tree_with_pheno_and_snp_mc5_mega_rooted.nw'
path_to_phenotypes = data_path + 'dr_covered_with_pheno_and_snp.csv'
out_tree = data_path + 'trees_mc5_mega_mix_rep/'
out_pheno = data_path + 'pheno_mc5_mega_mix_rep/'
out_sample_ids = data_path + '10drugs.sample_list'

# drug_names = ('Isoniazid', 'Rifampicin', 'Ethambutol', 'Pyrazinamide', 'Streptomycin')
drug_names = ('Isoniazid', 'Rifampicin', 'Ethambutol', 'Pyrazinamide', 'Streptomycin', 'Moxifloxacin', 'Ofloxacin',
              'Amikacin', 'Capreomycin', 'Kanamycin')

def process_pheno():
    all_pheno = pd.read_csv(path_to_phenotypes, sep='\t')
    for drug in all_pheno.columns[1:]:
        all_pheno[drug] = all_pheno[drug].map({'S': '0', 'R': '1'})

    all_pheno_no_na = all_pheno.dropna(subset=drug_names)
    for drug in drug_names:
        all_pheno_no_na[['Organism_name', drug]].to_csv(out_pheno + drug + '.pheno', sep='\t', index=False, header=False)

    return all_pheno_no_na['Organism_name'].tolist()


def main():
    if not exists(out_tree):
        makedirs(out_tree)
    if not exists(out_pheno):
        makedirs(out_pheno)
    sample_ids = process_pheno()
    with open(out_sample_ids, 'w') as f:
        for id in sample_ids:
            f.write(id + '\n')
    t = Tree(path_to_tree)
    t.prune(sample_ids)
    for drug in drug_names:
        t.write(format=1, outfile=out_tree + drug + ".nw")


if __name__ == '__main__':
    main()
