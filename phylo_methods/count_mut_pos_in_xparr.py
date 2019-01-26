import operator
import os
from os.path import exists

from src.core.constants import data_path, ref_len
from src.core.data_reading import read_h37rv
from src.phylo_methods.convert_coordinates_na_to_aa import localize_vars, read_Walker, print_pos

path_to_xparr = data_path + 'xparr/mc10_mega_MP/'
path_to_indel_list = data_path + 'xparr/mc10_mega_MP/indel_list.txt'
out_path = '../../res/mc10_mega_MP_stat_new/'

use_indels = True


def read_xparr(path_to_file):
    muts = set()
    mut_pos = {}
    mut_neg = {}
    with open(path_to_file, 'r') as f:
        for line in f.readlines()[2:]:
            s = line.split('\t')
            if len(s) == 6:
                if s[-1].strip() != '':
                    mut_list = s[-1].strip().split(';')
                    if s[-2] != '':
                        for m in mut_list:
                            if not m[0].isdigit():
                                m = m[1:-1]
                            m = int(m)
                            c = mut_pos.get(m)
                            muts.add(m)
                            if c is None:
                                mut_pos[m] = 1
                            else:
                                mut_pos[m] = c + 1
                    else:
                        for m in mut_list:
                            if not m[0].isdigit():
                                m = m[1:-1]
                            m = int(m)
                            c = mut_neg.get(m)
                            muts.add(m)
                            if c is None:
                                mut_neg[m] = 1
                            else:
                                mut_neg[m] = c + 1
    return list(muts), mut_pos, mut_neg


if __name__ == '__main__':
    if not exists(out_path):
        os.makedirs(out_path)
    if use_indels:
        indel_list = [l.strip() for l in open(path_to_indel_list, 'r').readlines()]
    else:
        indel_list = None
    gene_to_Walker = read_Walker()
    for (dirpath, dirnames, filenames) in os.walk(path_to_xparr):
        for filename in filenames:
            i = filename.rfind('.xparr')
            if i != -1:
                drug = filename[0: i]
                muts, mut_pos, mut_neg = read_xparr(path_to_xparr + filename)
                muts.sort()
                pos_to_cds, indel_pos_to_cds = localize_vars(muts, indel_list)
                with open(out_path + drug + '.stat', 'w') as out_f:
                    for pos in muts:
                        if pos > ref_len:
                            indel_pos = int(indel_list[pos - ref_len - 1].split('_')[1])
                            cds = indel_pos_to_cds.get(indel_pos)
                            p = mut_pos.get(pos)
                            if p is None:
                                p = 0
                            n = mut_neg.get(pos)
                            if n is None:
                                n = 0
                            print_pos(out_f, cds, indel_pos, gene_to_Walker, str(p + n) + '\t' + str(p) +
                                      '\t' + str(n), 'indel')
                        else:
                            cds = pos_to_cds.get(pos)
                            p = mut_pos.get(pos)
                            if p is None:
                                p = 0
                            n = mut_neg.get(pos)
                            if n is None:
                                n = 0
                            print_pos(out_f, cds, pos, gene_to_Walker, str(p + n) + '\t' + str(p) + '\t' + str(n), 'snp')