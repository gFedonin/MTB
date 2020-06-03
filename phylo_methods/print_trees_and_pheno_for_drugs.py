import pandas as pd
from ete3 import Tree
from os.path import exists

from os import makedirs, listdir
from sklearn.externals.joblib import Parallel, delayed

from core.constants import data_path

# path_to_tree = data_path + 'tree_with_pheno_and_snp_mc5_mega_rooted.nw'
path_to_tree = data_path + '/iqtree_with_indels/combined_mq40_keep_complex_std_names_filtered_no_DR_partitions_rooted.nw'
# path_to_phenotypes = data_path + 'dr_covered_with_pheno_and_snp.csv'
path_to_phenotypes = data_path + 'combined_phenotypes_filtered.csv'
# out_tree = data_path + 'trees_mc5_mega/'
out_tree = data_path + 'trees_combined_with_indels/'
# out_pheno = data_path + 'pheno_mc5_mega/'
out_pheno = data_path + 'pheno_combined/'

no_missing_values = True


def prune(sample_ids, drug):
    t = Tree(path_to_tree)
    samples_in_tree = set()
    for node in t.traverse("postorder"):
        samples_in_tree.add(node.name)
    # t.set_outgroup('canetti')
    sample_ids_filtered = [sid for sid in sample_ids if sid in samples_in_tree]
    t.prune(sample_ids_filtered)
    t.write(format=1, outfile=out_tree + drug + ".nw")
    return 0


def process_pheno():
    all_pheno = pd.read_csv(path_to_phenotypes, sep='\t')
    drug_samples_pairs = []
    for drug in all_pheno.columns[1:]:
        all_pheno[drug] = all_pheno[drug].map({'S': '0', 'R': '1'})
        drug_pheno = all_pheno[all_pheno[drug].isin(('0', '1'))]
        drug_pheno[['Organism_name', drug]].to_csv(out_pheno + drug + '.pheno', sep='\t', index=False, header=False)
        drug_samples_pairs.append((drug, drug_pheno['Organism_name'].tolist()))

    if no_missing_values:
        all_pheno_no_na = all_pheno.dropna(subset=['Isoniazid', 'Rifampicin'])
        mdr = all_pheno_no_na.loc[(all_pheno['Isoniazid'] == '1') & (all_pheno['Rifampicin'] == '1')]
        non_mdr = all_pheno_no_na.loc[(all_pheno['Isoniazid'] == '0') | (all_pheno['Rifampicin'] == '0')]
    else:
        mdr = all_pheno.loc[(all_pheno['Isoniazid'] == '1') & (all_pheno['Rifampicin'] == '1')]
        non_mdr = all_pheno.loc[(all_pheno['Isoniazid'] == '0') | (all_pheno['Rifampicin'] == '0')]
    sample_ids = []
    with open(out_pheno + 'MDR.pheno', 'w') as f:
        for index, row in mdr.iterrows():
            sample_ids.append(row['Organism_name'])
            f.write(row['Organism_name'] + '\t1\n')
        for index, row in non_mdr.iterrows():
            sample_ids.append(row['Organism_name'])
            f.write(row['Organism_name'] + '\t0\n')
    drug_samples_pairs.append(('MDR', sample_ids))

    sample_ids = []
    if no_missing_values:
        all_pheno_no_na = all_pheno.dropna(subset=['Isoniazid', 'Rifampicin', 'Ciprofloxacin', 'Moxifloxacin',
                                                   'Ofloxacin', 'Amikacin', 'Capreomycin', 'Kanamycin', 'Ciprofloxacin',
                                                   'Moxifloxacin', 'Ofloxacin', 'Amikacin', 'Capreomycin', 'Kanamycin'])
        xdr = all_pheno_no_na.loc[(all_pheno['Isoniazid'] == '1') & (all_pheno['Rifampicin'] == '1') &
                      ((mdr['Ciprofloxacin'] == '1') | (mdr['Moxifloxacin'] == '1') | (mdr['Ofloxacin'] == '1')) &
                      ((mdr['Amikacin'] == '1') | (mdr['Capreomycin'] == '1') | (mdr['Kanamycin'] == '1'))]
        non_xdr = all_pheno_no_na.loc[((all_pheno['Ciprofloxacin'] == '0') & (all_pheno['Moxifloxacin'] == '0') &
                                 (all_pheno['Ofloxacin'] == '0')) | ((all_pheno['Amikacin'] == '0') &
                                                                     (all_pheno['Capreomycin'] == '0') & (
                                                                                 all_pheno['Kanamycin'] == '0'))
                                | (all_pheno['Isoniazid'] == '0') | (all_pheno['Rifampicin'] == '0')]
    else:
        xdr = mdr.loc[((mdr['Ciprofloxacin'] == '1') | (mdr['Moxifloxacin'] == '1') | (mdr['Ofloxacin'] == '1')) &
                      ((mdr['Amikacin'] == '1') | (mdr['Capreomycin'] == '1') | (mdr['Kanamycin'] == '1'))]
        non_xdr = all_pheno.loc[((all_pheno['Ciprofloxacin'] == '0') & (all_pheno['Moxifloxacin'] == '0') &
                                 (all_pheno['Ofloxacin'] == '0')) | ((all_pheno['Amikacin'] == '0') &
                                                                     (all_pheno['Capreomycin'] == '0') & (
                                                                                 all_pheno['Kanamycin'] == '0'))
                                | (all_pheno['Isoniazid'] == '0') | (all_pheno['Rifampicin'] == '0')]
    with open(out_pheno + 'XDR.pheno', 'w') as f:
        for index, row in xdr.iterrows():
            sample_ids.append(row['Organism_name'])
            f.write(row['Organism_name'] + '\t1\n')
        for index, row in non_xdr.iterrows():
            sample_ids.append(row['Organism_name'])
            f.write(row['Organism_name'] + '\t0\n')
    drug_samples_pairs.append(('XDR', sample_ids))

    if no_missing_values:
        all_pheno_no_na = all_pheno.dropna(subset=['Isoniazid', 'Rifampicin'])
        non_mdr = all_pheno_no_na.loc[(all_pheno['Isoniazid'] == '0') & (all_pheno['Rifampicin'] == '0')]
    else:
        non_mdr = all_pheno.loc[(all_pheno['Isoniazid'] == '0') & (all_pheno['Rifampicin'] == '0')]
    sample_ids = []
    with open(out_pheno + 'MDR_S.pheno', 'w') as f:
        for index, row in mdr.iterrows():
            sample_ids.append(row['Organism_name'])
            f.write(row['Organism_name'] + '\t1\n')
        for index, row in non_mdr.iterrows():
            sample_ids.append(row['Organism_name'])
            f.write(row['Organism_name'] + '\t0\n')
    drug_samples_pairs.append(('MDR_S', sample_ids))

    if no_missing_values:
        all_pheno_no_na = all_pheno.dropna(subset=['Isoniazid', 'Rifampicin', 'Ciprofloxacin', 'Moxifloxacin',
                                                   'Ofloxacin', 'Amikacin', 'Capreomycin', 'Kanamycin', 'Ciprofloxacin',
                                                   'Moxifloxacin', 'Ofloxacin', 'Amikacin', 'Capreomycin', 'Kanamycin'])
        non_xdr = all_pheno_no_na.loc[((all_pheno['Ciprofloxacin'] == '0') & (all_pheno['Moxifloxacin'] == '0') &
                                 (all_pheno['Ofloxacin'] == '0')) & ((all_pheno['Amikacin'] == '0') &
                                                                     (all_pheno['Capreomycin'] == '0') & (
                                                                                 all_pheno['Kanamycin'] == '0'))
                                & (all_pheno['Isoniazid'] == '0') & (all_pheno['Rifampicin'] == '0')]
    else:
        non_xdr = all_pheno.loc[((all_pheno['Ciprofloxacin'] == '0') & (all_pheno['Moxifloxacin'] == '0') &
                                 (all_pheno['Ofloxacin'] == '0')) & ((all_pheno['Amikacin'] == '0') &
                                                                     (all_pheno['Capreomycin'] == '0') & (
                                                                             all_pheno['Kanamycin'] == '0'))
                                & (all_pheno['Isoniazid'] == '0') & (all_pheno['Rifampicin'] == '0')]
    with open(out_pheno + 'XDR_S.pheno', 'w') as f:
        for index, row in xdr.iterrows():
            sample_ids.append(row['Organism_name'])
            f.write(row['Organism_name'] + '\t1\n')
        for index, row in non_xdr.iterrows():
            sample_ids.append(row['Organism_name'])
            f.write(row['Organism_name'] + '\t0\n')
    drug_samples_pairs.append(('XDR', sample_ids))

    return drug_samples_pairs


def merge_pheno():
    path_to_pheno = data_path + 'combined_pheno/'
    sample_to_pheno = {}
    all_drugs = []
    for fname in listdir(path_to_pheno):
        drug = fname[:-6]
        all_drugs.append(drug)
        for l in open(path_to_pheno + fname).readlines():
            s = l.strip().split('\t')
            if s[1] == '0':
                pheno = 'S'
            else:
                pheno = 'R'
            drug_to_pheno = sample_to_pheno.get(s[0])
            if drug_to_pheno is None:
                sample_to_pheno[s[0]] = {drug: pheno}
            else:
                drug_to_pheno[drug] = pheno
    t = Tree(path_to_tree)
    samples_in_tree = set()
    for node in t.traverse("postorder"):
        samples_in_tree.add(node.name)
    with open(path_to_phenotypes, 'w') as f:
        f.write('Organism_name\t' + '\t'.join(all_drugs) + '\n')
        for sample, drug_to_pheno in sample_to_pheno.items():
            if sample not in samples_in_tree:
                continue
            f.write(sample)
            for drug in all_drugs:
                pheno = drug_to_pheno.get(drug)
                if pheno is None:
                    f.write('\t-')
                else:
                    f.write('\t' + pheno)
            f.write('\n')


def main():
    if not exists(out_tree):
        makedirs(out_tree)
    if not exists(out_pheno):
        makedirs(out_pheno)
    merge_pheno()
    drug_samples_pairs = process_pheno()
    tasks = Parallel(n_jobs=-1)(delayed(prune)(sample_ids, drug) for drug, sample_ids in drug_samples_pairs)
    c = 0
    for task in tasks:
        c += task
    print(c)


if __name__ == '__main__':
    main()
