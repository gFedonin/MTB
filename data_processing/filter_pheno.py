path_to_phenotypes = '../../data/all_with_pheno.csv'
out_path = '../../data/dr_covered_with_pheno_and_snp.csv'
path_to_ids = '../../data/dr_covered_with_pheno_and_snp_new.txt'


if __name__ == '__main__':
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
