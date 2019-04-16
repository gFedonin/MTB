from os import makedirs
from os.path import exists

from src.core.annotations import read_annotations
from src.core.constants import ref_len, upstream_length, data_path
from src.phylo_methods.convert_coordinates_na_to_aa import path_to_indel_list, localize_vars, read_Walker, print_pos

path_to_drug_codes = data_path + 'xparr/mc10_mega_MP/drug_codes.txt'
path_to_res_dir = '../../res/tree_was/phen_conditioned/1/'
file_name = 'vars.phen.intragene.ord_pairs'

filter_out_phylogenetic_markers = False
phylogenetic_markers_list_path = data_path + 'phylogenetic_markers_converted.csv'
p_value_threshold = 0.05
use_fdr = True
fdr_threshold = 0.05
use_pvalue_and_fdr = True
upper_or_lower = '.upper'

if use_fdr:
    if use_pvalue_and_fdr:
        filter_str = 'filtered_fdr' + str(fdr_threshold) + '_pvalue' + str(p_value_threshold)
    else:
        filter_str = 'filtered_fdr' + str(fdr_threshold)
else:
    filter_str = 'filtered_pvalue' + str(p_value_threshold)
if use_fdr:
    if use_pvalue_and_fdr:
        f_str = '.pvalue.pairs.fdr' + str(fdr_threshold) + '.pvalue' + str(p_value_threshold)
    else:
        f_str = '.pvalue.pairs.fdr' + str(fdr_threshold)
else:
    f_str = '.pvalue.pairs.pvalue' + str(p_value_threshold)
if filter_out_phylogenetic_markers:
    phylo_str = '.no_phylogenetic_markers'
else:
    phylo_str = ''


def read_drug_codes():
    number_to_drug = {}
    for line in open(path_to_drug_codes).readlines():
        s = line.strip().split('\t')
        number_to_drug[s[1]] = s[0]
    return number_to_drug


def find_max_pval_by_fdr():
    fin = None
    if exists(path_to_res_dir + file_name + upper_or_lower + '.pvalue.pairs.fdr'):
        fin = open(path_to_res_dir + file_name + upper_or_lower + '.pvalue.pairs.fdr')
    elif exists(path_to_res_dir + file_name + upper_or_lower + '.pvalue.sites.fdr'):
        fin = open(path_to_res_dir + file_name + upper_or_lower + '.pvalue.sites.fdr')
    if fin is not None:
        max_pval = -1
        for line in fin.readlines()[1:]:
            s = line.split('\t')
            if float(s[2])/float(s[1]) < fdr_threshold:
                pval = float(s[0])
                if pval > max_pval:
                    max_pval = pval
        return max_pval
    else:
        print('No such file: ' + path_to_res_dir + file_name + upper_or_lower + '.pvalue.pairs.fdr')
        print('No such file: ' + path_to_res_dir + file_name + upper_or_lower + '.pvalue.sites.fdr')
        return -1


def print_filtered_by_fdr():
    if not exists(path_to_res_dir + filter_str + upper_or_lower):
        makedirs(path_to_res_dir + filter_str + upper_or_lower)
    phylogenetic_markers_pos_set = None
    if filter_out_phylogenetic_markers:
        phylogenetic_markers_pos_set = set()
        with open(phylogenetic_markers_list_path) as fin:
            for line in fin.readlines():
                phylogenetic_markers_pos_set.add(line.split('\t')[0])

    max_pval = find_max_pval_by_fdr()
    if max_pval == -1:
        print('no pval with such fdr!')
        return
    if use_pvalue_and_fdr and max_pval > p_value_threshold:
        max_pval = p_value_threshold
    path_to_coords = path_to_res_dir + file_name + '.site_pairs'
    path_to_p_values = path_to_res_dir + file_name + upper_or_lower + '.pvalue'
    out_path = path_to_res_dir + filter_str + upper_or_lower + '/' + file_name + phylo_str + upper_or_lower + '.pvalue.pairs'
    with open(out_path, 'w') as f:
        p_val_iter = open(path_to_p_values)
        coords = open(path_to_coords)
        first_line = coords.readline().strip()
        f.write(first_line + '\tp_value\n')
        for line in coords.readlines():
            p_val = p_val_iter.readline().strip()
            if float(p_val) <= max_pval:
                if phylogenetic_markers_pos_set is not None:
                    s = line.split('\t')[0]
                    if s not in phylogenetic_markers_pos_set:
                        f.write(line.strip() + '\t' + p_val + '\n')
                else:
                    f.write(line.strip() + '\t' + p_val + '\n')


def print_filtered_by_p_value():
    if not exists(path_to_res_dir + filter_str + upper_or_lower):
        makedirs(path_to_res_dir + filter_str + upper_or_lower)
    phylogenetic_markers_pos_set = None
    if filter_out_phylogenetic_markers:
        phylogenetic_markers_pos_set = set()
        with open(phylogenetic_markers_list_path) as fin:
            for line in fin.readlines():
                phylogenetic_markers_pos_set.add(line.split('\t')[0])
    path_to_coords = path_to_res_dir + file_name + '.site_pairs'
    path_to_p_values = path_to_res_dir + file_name + upper_or_lower + '.pvalue'
    if not exists(path_to_coords) or not exists(path_to_p_values):
        print('no such file: ' + path_to_p_values)
    out_path = path_to_res_dir + filter_str + upper_or_lower + '/' + file_name + \
               phylo_str + upper_or_lower + '.pvalue.pairs'
    with open(out_path, 'w') as f:
        p_val_iter = open(path_to_p_values)
        coords = open(path_to_coords)
        first_line = coords.readline().strip()
        f.write(first_line + '\tp_value\n')
        for line in coords.readlines():
            p_val = p_val_iter.readline().strip()
            if float(p_val) <= p_value_threshold:
                if phylogenetic_markers_pos_set is not None:
                    s = line.split('\t')[0]
                    if s not in phylogenetic_markers_pos_set:
                        f.write(line.strip() + '\t' + p_val + '\n')
                else:
                    f.write(line.strip() + '\t' + p_val + '\n')


def convert(target_field, drug_field):
    if not exists(path_to_res_dir + filter_str + upper_or_lower + '_converted/'):
        makedirs(path_to_res_dir + filter_str + upper_or_lower + '_converted/')
    indel_list = [l.strip() for l in open(path_to_indel_list, 'r').readlines()]
    cds_list = read_annotations(upstream_length)
    path_to_coords = path_to_res_dir + filter_str + upper_or_lower + '/' + file_name + \
                     phylo_str + upper_or_lower + '.pvalue.pairs'
    if exists(path_to_coords):
        out_path = path_to_res_dir + filter_str + upper_or_lower + '_converted/' + \
                   file_name + phylo_str + f_str + '.filtered.converted'
        snps = [int(line.strip().split('\t')[0]) for line in open(path_to_coords, 'r').readlines()[1:] if line != '\n']
        pos_to_cds, indel_pos_to_cds = localize_vars(snps, cds_list, indel_list)
        gene_to_Walker = read_Walker()
        number_to_drug = read_drug_codes()
        with open(out_path, 'w') as out_f:
            with open(path_to_coords, 'r') as f:
                drug_code_pos = -1
                target_pos = -1
                header = f.readline().strip()
                s = header.split('\t')
                for i in range(len(s)):
                    if s[i] == drug_field:
                        drug_code_pos = i
                    if s[i] == target_field:
                        target_pos = i
                out_f.write('type\tgene_name\tgene_pos\t' + header + '\tWalker\tleft_gene\tright_gene\n')
                for line in f.readlines():
                    if line != '\n':
                        s = line.strip().split('\t')
                        pos = int(s[target_pos])
                        if drug_code_pos != -1:
                            s[drug_code_pos] = number_to_drug[s[drug_code_pos]]
                        if pos > ref_len:
                            l = '\t'.join(s)
                            s = indel_list[pos - ref_len - 1].split('_')
                            del_pos = int(s[1])
                            cds = indel_pos_to_cds.get(del_pos)
                            print_pos(out_f, cds, del_pos, gene_to_Walker, l, s[0], cds_list)
                        else:
                            cds = pos_to_cds.get(pos)
                            print_pos(out_f, cds, pos, gene_to_Walker, '\t'.join(s), 'snp', cds_list)


if __name__ == '__main__':
    if use_fdr:
        find_max_pval_by_fdr()
        print_filtered_by_fdr()
        convert('target site', None)#'target site' or 'bgr_site' or 'bg site'
    else:
        print_filtered_by_p_value()
        convert('target site', None)