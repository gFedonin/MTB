import operator
import os
from os.path import exists

from src.core.annotations import read_annotations
from src.core.constants import data_path, upstream_length
from src.core.data_reading import read_h37rv
from src.phylo_methods.association_from_xparr import read_xparr
from src.phylo_methods.convert_coordinates_na_to_aa import localize_vars, read_Walker, print_pos

path_to_xparr = data_path + 'xparr/mc10_mega_snp_genes/'
path_to_gene_list = data_path + 'xparr/mc10_mega_snp_genes/gene_list.txt'
out_path = '../../res/mc10_mega_snp_genes_stat/'


if __name__ == '__main__':
    if not exists(out_path):
        os.makedirs(out_path)
    gene_list = [l.strip() for l in open(path_to_gene_list, 'r').readlines()]
    for (dirpath, dirnames, filenames) in os.walk(path_to_xparr):
        for filename in filenames:
            i = filename.rfind('.xparr')
            if i != -1:
                drug = filename[0: i]
                mut_stat, mut_pos, mut_neg = read_xparr(path_to_xparr + filename)
                with open(out_path + drug + '.stat', 'w') as out_f:
                    for m, c in mut_stat:
                        pos = int(m)
                        gene_name = gene_list[pos]
                        out_f.write(gene_name + '\t' + str(c) + '\t' + str(mut_pos[m]) +
                                  '\t' + str(mut_neg[m]) + '\n')