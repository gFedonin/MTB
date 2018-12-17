from os import makedirs

from os.path import exists

from src.core.constants import data_path

filter_list_path = data_path + 'phylogenetic_markers_converted.csv'
path_to_ids = data_path + 'all_with_pheno_and_snp.txt'
var_path = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_mc10/'
out_path = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_mc10_no_phylo_markers/'

if __name__ == '__main__':
    if not exists(out_path):
        makedirs(out_path)
    filter_set = set()
    with open(filter_list_path) as fin:
        for line in fin.readlines():
            filter_set.add(line.strip())
    with open(path_to_ids) as fin:
        for line in fin.readlines():
            mut_list = []
            c = 0
            with open(var_path + line.strip() + '.variants') as fin2:
                for l in fin2.readlines():
                    if l[:l.rfind('\t')] not in filter_set:
                        mut_list.append(l)
                    else:
                        c += 1
            print(line.strip() + ' ' + str(c))
            with open(out_path + line.strip() + '.variants', 'w') as fout:
                fout.write(''.join(mut_list))
