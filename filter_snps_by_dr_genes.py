from constants import dr_genes

path_to_snp = './data/snps/annotated_with_pheno_and_snp_mc5/'
path_to_snp_list = path_to_snp + 'all_snp_list.csv'
out_path = path_to_snp + 'dr_genes_snp_list.csv'

filter_set = set(dr_genes)
with open(path_to_snp_list, 'r') as fin:
    with open(out_path, 'w') as fout:
        for line in fin.readlines():
            s = line.split('\t')
            if s[0] == 'Gene' or s[0] == 'upstream':
                if s[1] in filter_set:
                    fout.write(line)
