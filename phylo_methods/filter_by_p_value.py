from os import makedirs
from os.path import exists

from src.core.annotations import read_annotations
from src.core.constants import ref_len, upstream_length, data_path
from src.phylo_methods.convert_coordinates_na_to_aa import path_to_indel_list, localize_vars, read_Walker, print_pos

path_to_drug_codes = data_path + 'xparr/mc10_mega_MP/drug_codes.txt'
drug_list = ['PZA', 'RIF', 'AMI', 'CAP', 'CIP', 'EMB', 'ETH', 'INH', 'KAN', 'MOX', 'OFX', 'PRO', 'STR']#
phylogenetic_markers_list_path = data_path + 'phylogenetic_markers_converted.csv'
path_to_res_dir = '../../res/tree_was/mc10_mega_MP_filter1/'
filter_out_phylogenetic_markers = False

p_value_threshold = 0.05
use_fdr = True
fdr_threshold = 0.1
use_pvalue_and_fdr = True
upper_or_lower = '.lower'

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
    with open(path_to_res_dir + 'max_pvals_fdr' + str(fdr_threshold) + upper_or_lower + '.list', 'w') as f:
        for drug in drug_list:
            fin = None
            if exists(path_to_res_dir + drug + '/' + drug.lower() + upper_or_lower + '.pvalue.pairs_filter.fdr'):
                fin = open(path_to_res_dir + drug + '/' + drug.lower() + upper_or_lower + '.pvalue.pairs_filter.fdr')
            elif exists(path_to_res_dir + drug + '/' + drug.lower() + upper_or_lower + '.pvalue.sites_filter.fdr'):
                fin = open(path_to_res_dir + drug + '/' + drug.lower() + upper_or_lower + '.pvalue.sites_filter.fdr')
            if fin is not None:
                print(drug)
                max_pval = -1
                for line in fin.readlines()[1:]:
                    s = line.split('\t')
                    if float(s[2])/float(s[1]) < fdr_threshold:
                        pval = float(s[0])
                        if pval > max_pval:
                            max_pval = pval
                if max_pval >= 0:
                    f.write(drug + '\t' + str(max_pval) + '\n')


def print_filtered_by_fdr():
    if not exists(path_to_res_dir + filter_str + upper_or_lower):
        makedirs(path_to_res_dir + filter_str + upper_or_lower)
    phylogenetic_markers_pos_set = None
    if filter_out_phylogenetic_markers:
        phylogenetic_markers_pos_set = set()
        with open(phylogenetic_markers_list_path) as fin:
            for line in fin.readlines():
                phylogenetic_markers_pos_set.add(line.split('\t')[0])
    for line in open(path_to_res_dir + 'max_pvals_fdr' + str(fdr_threshold) + upper_or_lower + '.list').readlines():
        s = line.split('\t')
        drug = s[0]
        max_pval = float(s[1])
        if use_pvalue_and_fdr and max_pval > p_value_threshold:
            max_pval = p_value_threshold
        print(drug)
        path_to_coords = path_to_res_dir + drug + '/' + drug.lower() + '.site_pairs'
        path_to_p_values = path_to_res_dir + drug + '/' + drug.lower() + upper_or_lower + '.pvalue'
        out_path = path_to_res_dir + filter_str + upper_or_lower + '/' + drug.lower() + \
                   phylo_str + upper_or_lower + '.pvalue.pairs'
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
    for drug in drug_list:
        print(drug)
        path_to_coords = path_to_res_dir + drug + '/' + drug.lower() + '.site_pairs'
        path_to_p_values = path_to_res_dir + drug + '/' + drug.lower() + upper_or_lower + '.pvalue'
        if not exists(path_to_coords) or not exists(path_to_p_values):
            continue
        out_path = path_to_res_dir + filter_str + upper_or_lower + '/' + drug.lower() + \
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


def convert_for_all_drugs():
    if not exists(path_to_res_dir + filter_str + upper_or_lower + '_converted/'):
        makedirs(path_to_res_dir + filter_str + upper_or_lower + '_converted/')
    indel_list = [l.strip() for l in open(path_to_indel_list, 'r').readlines()]
    cds_list = read_annotations(upstream_length)
    number_to_drug = read_drug_codes()
    for drug in drug_list:
        path_to_coords = path_to_res_dir + filter_str + upper_or_lower + '/' + drug.lower() + \
                         phylo_str + upper_or_lower + '.pvalue.pairs'
        if exists(path_to_coords):
            out_path = path_to_res_dir + filter_str + upper_or_lower + '_converted/' + \
                       drug.lower() + phylo_str + f_str + '.filtered.converted'
            snps = [int(line.strip().split('\t')[0]) for line in open(path_to_coords, 'r').readlines()[1:] if line != '\n']
            pos_to_cds, indel_pos_to_cds = localize_vars(snps, cds_list, indel_list)
            gene_to_Walker = read_Walker()
            with open(out_path, 'w') as out_f:
                with open(path_to_coords, 'r') as f:
                    drug_code_pos = -1
                    header = f.readline().strip()
                    s = header.split('\t')
                    for i in range(len(s)):
                        if s[i] == 'target site' or s[i] == 'bgr_site':
                            drug_code_pos = i
                            break
                    out_f.write('type\tgene_name\tgene_pos\t' + header + '\tWalker\tleft_gene\tright_gene\n')
                    for line in f.readlines():
                        if line != '\n':
                            s = line.strip().split('\t')
                            pos = int(s[0])
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


def merge_lists_for_all_drugs():
    out_path = path_to_res_dir + 'pairs' + upper_or_lower + '.' + filter_str + \
               phylo_str + '.converted.csv'
    number_to_drug = read_drug_codes()
    drug_code_pos = -1
    with open(out_path, 'w') as f:
        first = True
        for drug in drug_list:
            path = path_to_res_dir + filter_str + upper_or_lower + '_converted/' + \
                       drug.lower() + phylo_str + f_str + '.filtered.converted'
            if exists(path):
                fin = open(path)
                if first:
                    first = False
                    header = fin.readline()
                    s = header.split('\t')
                    for drug_code_pos in range(len(s)):
                        if s[drug_code_pos] == 'target site':
                            break
                    f.write(header)
                else:
                    fin.readline()
                for line in fin.readlines():
                    s = line.split('\t')
                    s[drug_code_pos] = number_to_drug[s[drug_code_pos]]
                    f.write('\t'.join(s))


if __name__ == '__main__':
    if use_fdr:
        find_max_pval_by_fdr()
        print_filtered_by_fdr()
        convert_for_all_drugs()
        merge_lists_for_all_drugs()
    else:
        print_filtered_by_p_value()
        convert_for_all_drugs()
        merge_lists_for_all_drugs()