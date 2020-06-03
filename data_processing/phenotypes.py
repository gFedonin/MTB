from os import listdir, makedirs
from os.path import exists

from core.constants import data_path

# path_to_pheno_csv = data_path + 'coll18/pheno.csv'
from core.data_reading import read_all_pheno

path_to_pheno_csv = data_path + 'phenotypes.csv'
# out_path = data_path + 'coll18/pheno/'
out_path = data_path + 'combined_pheno/'
path_to_ids = data_path + 'snps/combined_raw_variants_mq40_keep_complex_std_names_filtered/samples_filtered.list'


def split_pheno():
    if not exists(out_path):
        makedirs(out_path)
    sample_ids = set(l.strip() for l in open(path_to_ids).readlines())
    f = open(path_to_pheno_csv)
    s = f.readline().strip().split('\t')
    drugs = s[1:]
    drug_to_pheno = {drug: [] for drug in drugs}
    for l in f.readlines():
        s = l.strip().split('\t')
        if s[0] not in sample_ids:
            continue
        for i in range(1, len(s)):
            if s[i] in ('R', 'S') or s[i] in ('1', '0'):
                if s[i] == 'S':
                    s[i] = '0'
                elif s[i] == 'R':
                    s[i] = '1'
                else:
                    print(l)
                    exit(1)
                drug_to_pheno[drugs[i - 1]].append(s[0] + '\t' + s[i])
    f.close()
    for drug, sample_to_pheno in drug_to_pheno.items():
        with open(out_path + drug[0] + drug[1:].lower() + '.pheno', 'w') as f:
            f.write('\n'.join(sample_to_pheno))


path1 = data_path + 'coll18/pheno/'
path2 = data_path + 'pheno_mc5/'
merged_pheno_path = data_path + 'pheno_base_and_coll18/'


def merge_pheno():
    fnames1 = [fname for fname in listdir(path1)]
    fnames2 = set(fname for fname in listdir(path2))
    for fname in fnames1:
        if fname in fnames2:
            with open(merged_pheno_path + fname, 'w') as f:
                f.write(''.join(open(path1 + fname).readlines()))
                f.write('\n')
                f.write(''.join(open(path2 + fname).readlines()))
        else:
            with open(merged_pheno_path + fname, 'w') as f:
                f.write(''.join(open(path1 + fname).readlines()))
    fnames1 = set(fnames1)
    for fname in fnames2:
        if fname not in fnames1:
            with open(merged_pheno_path + fname, 'w') as f:
                f.write(''.join(open(path2 + fname).readlines()))


path_to_walker_pheno = 'd:\\MyProjects\\Neverov\\MTB\\readme\\data\\Walker2018\\pheno_big.csv'
path_to_id_mapping = 'd:\\MyProjects\\Neverov\\MTB\\readme\\data\\Walker2018\\lib_err_samea.csv'
out_path_pheno = 'd:\\MyProjects\\Neverov\\MTB\\readme\\data\\Walker2018\\pheno_big_filtered.csv'


def parse_walker():
    lib_to_samea = {}
    err_to_samea = {}
    for l in open(path_to_id_mapping).readlines():
        s = l.strip().split('\t')
        lib_to_samea[s[0]] = s[2]
        if ';' in s[1]:
            s[1] = s[1].split(';')[0]
        err_to_samea[s[1]] = s[2]
    with open(out_path_pheno, 'w') as f:
        for l in open(path_to_walker_pheno).readlines():
            s = l.strip().split()
            if s[2] == 'NA':
                if s[1] == 'NA':
                    samea = lib_to_samea.get(s[0])
                    if samea is None:
                        continue
                    s[2] = samea
                else:
                    samea = err_to_samea.get(s[1])
                    if samea is None:
                        continue
                    s[2] = samea
                f.write('\t'.join(s[2:]) + '\n')
            else:
                f.write('\t'.join(s[2:]) + '\n')


path_to_farhat_pheno = 'd:\\MyProjects\\Neverov\\MTB\\readme\\data\\Farhat2019\\pheno_mic.csv'
path_to_thresholds = 'd:\\MyProjects\\Neverov\\MTB\\readme\\data\\Farhat2019\\thresholds.csv'
path_to_filtered_list = 'd:\\MyProjects\\Neverov\\MTB\\readme\\data\\Farhat2019\\farhat19_new.samples'
out_farhat_pheno = 'd:\\MyProjects\\Neverov\\MTB\\readme\\data\\Farhat2019\\pheno.csv'


def convert_farhat():
    thresholds = {}
    for l in open(path_to_thresholds).readlines():
        s = l.strip().split('\t')
        thresholds[s[0]] = float(s[1])
    # samples = set(l.strip() for l in open(path_to_filtered_list).readlines())
    with open(out_farhat_pheno, 'w') as f:
        lines = open(path_to_farhat_pheno).readlines()
        drugs = lines[0].strip().split('\t')[1:]
        f.write(lines[0])
        for l in lines[1:]:
            s = l.strip().split('\t')
            # if s[0] not in samples:
            #     continue
            pheno = [s[0]]
            for i in range(len(drugs)):
                val = s[i + 1]
                if val == 'NA':
                    pheno.append('-')
                    continue
                drug = drugs[i]
                if val.startswith('<='):
                    if float(val[2:]) <= thresholds[drug]:
                        pheno.append('S')
                    else:
                        pheno.append('-')
                elif val.startswith('>'):
                    if thresholds[drug] <= float(val[1:]):
                        pheno.append('R')
                    else:
                        pheno.append('-')
                else:
                    j = val.find('-')
                    if j == -1:
                        print(l)
                    else:
                        lower = float(val[:j])
                        upper = float(val[j + 1:])
                        if upper <= thresholds[drug]:
                            pheno.append('S')
                        elif lower > thresholds[drug]:
                            pheno.append('R')
                        else:
                            pheno.append('-')
            f.write('\t'.join(pheno) + '\n')


path_to_dup = data_path + 'dup_snp_noDR.list'
path_to_pheno = data_path + 'combined_pheno/'
path_to_list_to_filter = data_path + 'combined_no_dup_with_pheno.list'
path_to_undup_list = data_path + 'combined_no_dup_snp_with_pheno.list'
path_to_best_in_clusters = data_path + 'best_in_clusters.list'


def filter_duplicates():
    clusters = []
    centroid = None
    cluster = None
    duplicates = set()
    for l in open(path_to_dup).readlines():
        s = l.strip().split(' ')
        if centroid != s[0]:
            if s[0] not in duplicates:
                centroid = s[0]
                cluster = s
                clusters.append(cluster)
                duplicates.update(s)
        else:
            cluster.append(s[1])
            duplicates.add(s[1])

    drug_to_pheno = read_all_pheno(path_to_pheno)
    drug_to_samples = {}
    drug_to_counts = []
    for drug, sample_to_pheno in drug_to_pheno.items():
        sample_set = {sample for sample, pheno in sample_to_pheno}
        drug_to_samples[drug] = sample_set
        drug_to_counts.append((drug, len(sample_set)))
    drug_to_counts.sort(key=lambda x: x[1])
    with open(path_to_best_in_clusters, 'w') as f:
        for cluster in clusters:
            best = cluster[0]
            best_drugs = {drug for drug, sample_set in drug_to_samples.items() if best in sample_set}
            for candidate in cluster[1:]:
                for drug, count in drug_to_counts:
                    if candidate in drug_to_samples[drug] and drug not in best_drugs:
                        best = candidate
                        best_drugs = {drug for drug, sample_set in drug_to_samples.items() if best in sample_set}
                        break
            duplicates.remove(best)
            f.write(best + '\n')
    with open(path_to_undup_list, 'w') as f:
        for l in open(path_to_list_to_filter).readlines():
            if l.strip() not in duplicates:
                f.write(l)


if __name__ == '__main__':
    # split_pheno()
    # merge_pheno()
    # parse_walker()
    # convert_farhat()
    filter_duplicates()