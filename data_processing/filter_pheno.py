# path_to_phenotypes = '../../data/all_with_pheno.csv'
from core.constants import data_path

path_to_phenotypes = data_path + 'combined_phenotypes.csv'
# out_path = '../../data/dr_covered_with_pheno_and_snp.csv'
out_path = data_path + 'combined_phenotypes_filtered.csv'
# path_to_ids = '../../data/dr_covered_with_pheno_and_snp_new.txt'
path_to_ids = data_path + 'snps/combined_raw_variants_mq40_keep_complex_std_names_filtered/samples_filtered.list'


def filter_by_list():
    with open(path_to_ids, 'r') as f:
        sample_ids = set()
        for line in f.readlines():
            sample_ids.add(line.strip())
    with open(path_to_phenotypes, 'r') as f:
        pheno = f.readlines()
    with open(out_path, 'w') as f:
       for line in pheno:
           if 'Organism_name' in line:
               f.write(line)
           s = line.split(sep='\t')
           if s[0] in sample_ids:
               f.write(line)


def filter_dup():
    lines = open(path_to_phenotypes).readlines()
    sample_counts = {}
    for line in lines[1:]:
        s = line.split('\t')
        if s[0] in sample_counts:
            sample_counts[s[0]] += 1
        else:
            sample_counts[s[0]] = 1
    for sample, c in sample_counts.items():
        if c > 1:
            print(sample)


if __name__ == '__main__':
    filter_dup()
