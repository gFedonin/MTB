path_to_phenotypes = './data/phenotype_filtered.csv'
out_path = './data/phenotype_filtered1.csv'
path_to_ids = './data/sample_ids_filtered2.txt'

if __name__ == '__main__':
    with open(path_to_ids, 'r') as f:
        sample_ids = set()
        for line in f.readlines():
            sample_ids.add(line[:-1])
    with open(path_to_phenotypes, 'r') as f:
        pheno = f.readlines()
    with open(out_path, 'w') as f:
       for line in pheno:
           if 'Organism_name' in line:
               f.write(line)
           s = line.split(sep='\t')
           if s[0] in sample_ids:
               f.write(line)
