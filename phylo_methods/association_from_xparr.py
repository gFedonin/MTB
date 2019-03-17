import operator
import os
from os.path import exists

from src.core.annotations import read_annotations
from src.core.constants import data_path, ref_len, upstream_length
from src.core.data_reading import read_h37rv
from src.phylo_methods.convert_coordinates_na_to_aa import localize_vars, read_Walker, print_pos

path_to_xparr = data_path + 'xparr/mc10_mega_MP/'
path_to_indel_list = data_path + 'xparr/mc10_mega_MP/indel_list.txt'
out_path = '../../res/mc10_mega_MP_stat/'

use_indels = True


def read_xparr(path_to_file):
    mut_pos = {}
    mut_neg = {}
    with open(path_to_file, 'r') as f:
        for line in f.readlines()[2:]:
            s = line.strip().split('\t')
            if s[5] != '':
                mut_list = s[5].split(';')
                if len(mut_list) > 0:
                    if s[4] != '':
                        for m in mut_list:
                            if not m[0].isdigit():
                                m = m[1:-1]
                            c = mut_pos.get(m)
                            if c is None:
                                mut_pos[m] = 1
                            else:
                                mut_pos[m] = c + 1
                    else:
                        for m in mut_list:
                            if not m[0].isdigit():
                                m = m[1:-1]
                            c = mut_neg.get(m)
                            if c is None:
                                mut_neg[m] = 1
                            else:
                                mut_neg[m] = c + 1
    mut_stat = []
    for m, c in mut_pos.items():
        n = mut_neg.get(m)
        if n is None:
            mut_stat.append((m, 1))
            mut_neg[m] = 0
        else:
            mut_stat.append((m, c/(c + n)))
    for m, c in mut_neg.items():
        n = mut_pos.get(m)
        if n is None:
            mut_stat.append((m, 0))
            mut_pos[m] = 0
    mut_stat.sort(key=operator.itemgetter(1), reverse=True)
    return mut_stat, mut_pos, mut_neg


if __name__ == '__main__':
    if not exists(out_path):
        os.makedirs(out_path)
    if use_indels:
        indel_list = [l.strip() for l in open(path_to_indel_list, 'r').readlines()]
    else:
        indel_list = None
    gene_to_Walker = read_Walker()
    cds_list = read_annotations(upstream_length)
    for (dirpath, dirnames, filenames) in os.walk(path_to_xparr):
        for filename in filenames:
            i = filename.rfind('.xparr')
            if i != -1:
                drug = filename[0: i]
                mut_stat, mut_pos, mut_neg = read_xparr(path_to_xparr + filename)
                snps = [int(m) for m,c in mut_stat]
                snps.sort()
                pos_to_cds, indel_pos_to_cds = localize_vars(snps, cds_list, indel_list)
                with open(out_path + drug + '.stat', 'w') as out_f:
                    for m, c in mut_stat:
                        pos = int(m)
                        if pos > ref_len:
                            del_pos = int(indel_list[pos - ref_len - 1].split('_')[1])
                            cds = indel_pos_to_cds[del_pos]
                            print_pos(out_f, cds, del_pos, gene_to_Walker, str(c) + '\t' + str(mut_pos[m]) +
                                      '\t' + str(mut_neg[m]), 'indel')
                        else:
                            cds = pos_to_cds.get(pos)
                            print_pos(out_f, cds, pos, gene_to_Walker, str(c) + '\t' + str(mut_pos[m]) +
                                      '\t' + str(mut_neg[m]), 'snp')