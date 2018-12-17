import os

from sklearn.externals.joblib import Parallel, delayed

path_to_snp_list = './data/snps/annotated_with_pheno_and_snp/all_snp_list.csv'
path_to_dict = './data/dictionaries1/'

drug_names = {'STR': 'Streptomycin', 'INH': 'Isoniazid', 'EMB': 'Ethambutol', 'ETH': 'Ethionamide',
              'CIP': 'Ciprofloxacin', 'OFX': 'Ofloxacin', 'PZA': 'Pyrazinamide', 'RIF': 'Rifampicin',
              'AMK': 'Amikacin', 'CAP': 'Capreomycin', 'KAN': 'Kanamycin', 'SM': 'Streptomycin',
              'MOX': 'Moxifloxacin', 'PTH': 'Prothionamide',
              'FLQ': ['Ciprofloxacin', 'Moxifloxacin', 'Ofloxacin'], 'AMI':'Amikacin'}


def read_snp_list():
    list = set()
    with open(path_to_snp_list, 'r') as f:
        for line in f.readlines():
            list.add(line.strip())
    return list


def read_snp_set():
    list = {}
    with open(path_to_snp_list, 'r') as f:
        for line in f.readlines():
            s = line.strip().split('\t')
            list['\t'.join(s[0:3])] = (s[3],s[4])
    return list


def read_dict(name):
    mut_list = set()
    with open(path_to_dict + name + '.txt', 'r') as f:
        for line in f.readlines():
            line = line.strip()
            j = line.rfind('\t')
            drug = line[j + 1:]
            if drug in drug_names.keys():
                full_name = drug_names[drug]
                if type(full_name) is list:
                    mut_list.add(line[: j])
                else:
                    mut_list.add(line[: j])
            else:
                print(drug + ' not found')
    return name, mut_list


def main():
    all_snps = read_snp_list()
    dict_name_to_list = {}
    # all_snps = read_snp_set()
    for (dirpath, dirnames, filenames) in os.walk(path_to_dict):
        tasks = Parallel(n_jobs=-1)(delayed(read_dict)(name) for name in [filename[:-4] for filename in filenames])
        for name, durg_to_mut_list in tasks:
            dict_name_to_list[name] = durg_to_mut_list
    for name, mut_list in dict_name_to_list.items():
        print(name)
        for f in mut_list:
            # s = f.split('\t')
            # if '\t'.join(s[0:3]) in all_snps.keys():
            #     ref, alt = all_snps['\t'.join(s[0:3])]
            #     if ref != s[3]:
            #         print(name + ' ' + f)
            if f not in all_snps:
                print(f)


if __name__ == '__main__':
    main()
