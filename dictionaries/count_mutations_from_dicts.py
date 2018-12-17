import os

from sklearn.externals.joblib import Parallel, delayed

from src.core.data_reading import read_pheno, read_dict, read_variants

path_to_pheno = '../../data/pheno_mc5_Walker/'
path_to_snp = '../../data/snps/annotated_with_pheno_and_snp_mc5_snp_only/'
path_to_dict = '../../data/dictionaries_new/'
path_to_ids = '../../data/dr_covered_with_pheno_and_snp.txt'
out_path = '../../res/dict_snp_counts_Walker_test.out'


def count_snps_presence(sample_to_snps, drug_to_pheno, dict_name_to_list):
    with open(out_path, 'w') as f:
        for dict_name, drug_to_list in dict_name_to_list.items():
            f.write(dict_name + '\n\n')
            for drug, snp_list in drug_to_list.items():
                mut_to_count = {}
                mut_to_R = {}
                mut_to_S = {}
                if drug not in drug_to_pheno.keys():
                    continue
                sample_to_pheno = drug_to_pheno[drug]
                for mut in snp_list:
                    mut_to_count[mut] = 0
                    mut_to_R[mut] = 0
                    mut_to_S[mut] = 0
                for sample_id, pheno in sample_to_pheno:
                    mut_list = sample_to_snps[sample_id]
                    for mut in mut_list:
                        # s = mut.split('\t')
                        # if len(s) == 6:
                        #     mut = '\t'.join(s[:-1])
                        c = mut_to_count.get(mut)
                        if c is not None:
                            mut_to_count[mut] = c + 1
                            if pheno == 1:
                                mut_to_R[mut] += 1
                            else:
                                mut_to_S[mut] += 1

                f.write(drug + '\n')
                for mut, c in mut_to_count.items():
                    f.write(mut + ' ' + str(c) + ' R: ' + str(mut_to_R[mut]) + ' S: ' + str(mut_to_S[mut]) + '\n')
                f.write('\n')


def main():

    drug_to_pheno = {}
    for (dirpath, dirnames, filenames) in os.walk(path_to_pheno):
        tasks = Parallel(n_jobs=-1)(delayed(read_pheno)(path_to_pheno, drug) for drug in [filename[:-6] for filename in filenames])
        for drug, sample_id_to_pheno in tasks:
            drug_to_pheno[drug] = sample_id_to_pheno
    print(str(len(drug_to_pheno.keys())) + ' drugs')

    dict_name_to_list = {}
    for (dirpath, dirnames, filenames) in os.walk(path_to_dict):
        tasks = Parallel(n_jobs=-1)(delayed(read_dict)(path_to_dict, name) for name in [filename[:-4] for filename in filenames])
        for name, drug_to_mut_list in tasks:
            dict_name_to_list[name] = drug_to_mut_list

    print(str(len(dict_name_to_list.keys())) + ' dictionaries')

    drug_to_full_list = {}
    for drug in drug_to_pheno.keys():
        full_set = set()
        for name, drug_to_mut_list in dict_name_to_list.items():
            if drug in drug_to_mut_list.keys():
                for f in drug_to_mut_list[drug]:
                    full_set.add(f)
        if len(full_set) > 0:
            drug_to_full_list[drug] = list(full_set)

    full_dict_set = set()
    for l in drug_to_full_list.values():
        for f in l:
            full_dict_set.add(f)

    with open(path_to_ids, 'r') as f:
        sample_ids = [name.strip() for name in f.readlines()]

    sample_to_snps = {}
    tasks = Parallel(n_jobs=-1)(delayed(read_variants)(path_to_snp, sample_id, keep_type=False) for sample_id in sample_ids)# , full_dict_set
    for name, l in tasks:
        sample_to_snps[name] = l

    print(str(len(sample_ids)) + ' samples after filtering')

    # full_snp_list = read_snp_list(path_to_snp_list)
    full_snp_list = list(full_dict_set)
    snp_num = len(full_snp_list)
    print(str(snp_num) + ' snps')
    count_snps_presence(sample_to_snps, drug_to_pheno, dict_name_to_list)


if __name__ == '__main__':
    main()

