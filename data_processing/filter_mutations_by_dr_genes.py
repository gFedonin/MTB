from src.core.constants import dr_genes

path_to_snp = '../../data/snps/annotated_with_DR_with_indel_with_pheno_and_snp_mc10_long_del/'
path_to_snp_list = path_to_snp + 'all_var_list.csv'
out_path = path_to_snp + 'dr_genes_var_list.csv'

filter_set = set(dr_genes)
with open(path_to_snp_list, 'r') as fin:
    with open(out_path, 'w') as fout:
        for line in fin.readlines():
            s = line.split('\t')
            if s[0] == 'Gene' or s[0] == 'upstream':
                if s[1] in filter_set:
                    fout.write(line)
