import random
from bisect import bisect_left, bisect_right
import os
from math import floor

from os.path import exists

# from ete3 import Tree
from Bio import SeqIO
from sklearn.externals.joblib import Parallel, delayed

from src.core.annotations import read_annotations
from src.core.constants import codon_table, codon_table_compl, data_path, upstream_length

drug = 'Rifampicin'
mut = 'Gene\trpoB\t450'
forward_tag = '[&!color=#ff0033]'
backward_tag = '[&!color=#00ff00]'

out_path = data_path + 'colored_trees/'

path_to_ids = data_path + 'dr_covered_with_pheno_and_snp.txt'
path_to_snps = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_mc10/'
path_to_alignment = data_path + 'ancestors_mc10_mega_merged.fasta'
path_to_pheno = data_path + 'pheno_mc5_mega/'
path_to_pheno_and_trees = data_path + 'reconstructed_mc10_mega_MP/'
path_to_ref = data_path + 'h37rv.fasta'

path_to_fig_tree = data_path + 'fig_tree.txt'

color_only_changes = True
prune = False
random_fraction = 0.1


def read_h37rv():
    fasta_sequences = SeqIO.parse(open(path_to_ref), 'fasta')
    sequence = ''
    for fasta in fasta_sequences:
        sequence = str(fasta.seq)
    return sequence.upper()


def get_aminoacids_sense(child_seq, parent_seq, ref_seq, nucleotide_pos, snps, i):
    snp_num = len(snps)
    pos = snps[i]
    alt = child_seq[i]
    if color_only_changes:
        par = parent_seq[i]
    else:
        par = ref_seq[pos - 1]
    if nucleotide_pos == 0:
        if i != snp_num - 1:
            pos1 = snps[i + 1]
            alt1 = child_seq[i + 1]
            if color_only_changes:
                par1 = parent_seq[i + 1]
            else:
                par1 = ref_seq[pos1 - 1]
            if pos1 == pos + 1:
                if i != snp_num - 2 and snps[i + 2] == pos + 2:
                    if color_only_changes:
                        aa0 = codon_table[par + par1 + parent_seq[i + 2]]
                    else:
                        aa0 = codon_table[par + par1 + ref_seq[pos + 1]]
                    aa1 = codon_table[alt + alt1 + child_seq[i + 2]]
                    return aa0, aa1, i + 2
                else:
                    aa0 = codon_table[par + par1 + ref_seq[pos + 1]]
                    aa1 = codon_table[alt + alt1 + ref_seq[pos + 1]]
                    return aa0, aa1, i + 1
            elif pos1 == pos + 2:
                aa0 = codon_table[par + ref_seq[pos] + par1]
                aa1 = codon_table[alt + ref_seq[pos] + alt1]
                return aa0, aa1, i + 1
        aa0 = codon_table[par + ref_seq[pos] + ref_seq[pos + 1]]
        aa1 = codon_table[alt + ref_seq[pos] + ref_seq[pos + 1]]
        return aa0, aa1, i
    elif nucleotide_pos == 1:
        if i != snp_num - 1 and snps[i + 1] == pos + 1:
            if color_only_changes:
                aa0 = codon_table[ref_seq[pos - 2] + par + parent_seq[i + 1]]
            else:
                aa0 = codon_table[ref_seq[pos - 2] + par + ref_seq[pos]]
            aa1 = codon_table[ref_seq[pos - 2] + alt + child_seq[i + 1]]
            return aa0, aa1, i + 1
        else:
            aa0 = codon_table[ref_seq[pos - 2] + par + ref_seq[pos]]
            aa1 = codon_table[ref_seq[pos - 2] + alt + ref_seq[pos]]
            return aa0, aa1, i
    else:
        aa0 = codon_table[ref_seq[pos - 3:pos - 1] + par]
        aa1 = codon_table[ref_seq[pos - 3:pos - 1] + alt]
        return aa0, aa1, i


def get_aminoacids_antisense(child_seq, parent_seq, ref_seq, nucleotide_pos, snps, i):
    snp_num = len(snps)
    pos = snps[i]
    alt = child_seq[i]
    if color_only_changes:
        par = parent_seq[i]
    else:
        par = ref_seq[pos - 1]
    if nucleotide_pos == 2:
        if i != snp_num - 1:
            pos1 = snps[i + 1]
            alt1 = child_seq[i + 1]
            if color_only_changes:
                par1 = parent_seq[i + 1]
            else:
                par1 = ref_seq[pos1 - 1]
            if pos1 == pos + 1:
                if i != snp_num - 2 and snps[i + 2] == pos + 2:
                    if color_only_changes:
                        aa0 = codon_table_compl[par + par1 + parent_seq[i + 2]]
                    else:
                        aa0 = codon_table_compl[par + par1 + ref_seq[pos + 1]]
                    aa1 = codon_table_compl[alt + alt1 + child_seq[i + 2]]
                    return aa0, aa1, i + 2
                else:
                    aa0 = codon_table_compl[par + par1 + ref_seq[pos + 1]]
                    aa1 = codon_table_compl[alt + alt1 + ref_seq[pos + 1]]
                    return aa0, aa1, i + 1
            elif pos1 == pos + 2:
                aa0 = codon_table_compl[par + ref_seq[pos] + par1]
                aa1 = codon_table_compl[alt + ref_seq[pos] + alt1]
                return aa0, aa1, i + 1
        aa0 = codon_table_compl[par + ref_seq[pos:pos + 2]]
        aa1 = codon_table_compl[alt + ref_seq[pos:pos + 2]]
        return aa0, aa1, i
    elif nucleotide_pos == 1:
        if i != snp_num - 1 and snps[i + 1] == pos + 1:
            if color_only_changes:
                aa0 = codon_table_compl[ref_seq[pos - 2] + par + parent_seq[i + 1]]
            else:
                aa0 = codon_table_compl[ref_seq[pos - 2] + par + ref_seq[pos]]
            aa1 = codon_table_compl[ref_seq[pos - 2] + alt + child_seq[i + 1]]
            return aa0, aa1, i + 1
        aa0 = codon_table_compl[ref_seq[pos - 2] + par + ref_seq[pos]]
        aa1 = codon_table_compl[ref_seq[pos - 2] + alt + ref_seq[pos]]
        return aa0, aa1, i
    else:
        aa0 = codon_table_compl[ref_seq[pos - 3:pos - 1] + par]
        aa1 = codon_table_compl[ref_seq[pos - 3:pos - 1] + alt]
        return aa0, aa1, i


def localize_all_snps(all_snps, cds_list):
    res = {}
    s = mut.split('\t')
    for cds in cds_list:
        if cds.type.name == s[0] and cds.name == s[1]:
            i = bisect_left(all_snps, cds.start)
            j = bisect_right(all_snps, cds.end, lo=i)
            for k in range(i, j):
                res[all_snps[k]] = cds
            break
    return res


def format_snp(fasta_seq, parent_seq, snp_pos_list, ref_seq, snp_to_cds):
    s = mut.split('\t')
    mut_pos = int(s[2])
    # print('looking for ' + str(mut_pos))
    j = 0
    while j < len(snp_pos_list):
        if color_only_changes and fasta_seq[j] == parent_seq[j]:
            j += 1
            continue
        pos = snp_pos_list[j]
        cds = snp_to_cds.get(pos)

        # if cds is not None:
        if cds.strand == 1:
            protein_pos = (pos - cds.start) // 3 + 1
            nucleotide_pos = (pos - cds.start) % 3
            if protein_pos == mut_pos:
                aa0, aa1, j = get_aminoacids_sense(fasta_seq, parent_seq, ref_seq, nucleotide_pos, snp_pos_list, j)
                if aa0 != aa1:
                    return 1
                else:
                    return 0

        else:  # if strand is '-'
            protein_pos = (cds.end - pos) // 3 + 1
            nucleotide_pos = (cds.end - pos) % 3
            print(protein_pos)
            if protein_pos == mut_pos:
                aa0, aa1, j = get_aminoacids_antisense(fasta_seq, parent_seq, ref_seq, nucleotide_pos, snp_pos_list, j)
                if aa0 != aa1:
                    return 1
                else:
                    return 0
        j += 1
    return 0


def read_snps(sample_id):
    snps = {}
    with open(path_to_snps + sample_id + '.variants', 'r') as f1:
        for line in f1.readlines():
            s = line.strip().split('\t')
            snps[int(s[0])] = s[1]
    return sample_id, snps


def filter_snp_list(snp_pos_list, snp_to_cds):
    index_list = []
    pos_list = []
    for i in range(len(snp_pos_list)):
        if snp_to_cds.get(snp_pos_list[i]) is not None:
            index_list.append(i)
            pos_list.append(snp_pos_list[i])
    return index_list, pos_list


def filter_aln(sample_id, seq, index_list):
    res = []
    for i in index_list:
        res.append(seq[i])
    return sample_id, ''.join(res)


def read_parents(drug):
    parents = []
    with open(path_to_pheno_and_trees + drug + '/parents.csv', 'r') as f:
        root = f.readline().strip()
        for line in f.readlines():
            s = line.strip().split('\t')
            parents.append((s[0], s[1], s[2]))
    return root, parents


def format_mut_lists(sample_to_seq, pos_list, ref_seq, snp_to_cds, parents):
    forward = set()
    backward = set()
    for node_id, parent_id, dist in parents:
        print(node_id)
        c = format_snp(sample_to_seq[node_id], sample_to_seq[parent_id], pos_list, ref_seq, snp_to_cds)
        if c == 1:
            forward.add(node_id)
            # print('+')
        elif c == -1:
            backward.add(node_id)
        # else:
        #     print('-')
    return forward, backward


def read_pheno(drug):

    sample_to_pheno = {}
    with open(path_to_pheno + drug + '.pheno', 'r') as f:
        for line in f.readlines():
            s = line.strip().split('\t')
            if s[1] == '1':
                sample_to_pheno[s[0]] = 'R'
            else:
                sample_to_pheno[s[0]] = 'S'
    with open(path_to_pheno_and_trees + drug + '/anc.pheno', 'r') as f:
        for line in f.readlines():
            s = line.strip().split('\t')
            if s[1] == '1':
                sample_to_pheno[s[0]] = 'R'
            else:
                sample_to_pheno[s[0]] = 'S'
    return sample_to_pheno


def main():

    if not exists(out_path):
        os.makedirs(out_path)

    h37rv = read_h37rv()

    sample_ids = [sample_id.strip() for sample_id in open(path_to_ids, 'r').readlines()]
    tasks = Parallel(n_jobs=-1)(delayed(read_snps)(sample_id) for sample_id in sample_ids)
    all_snp_pos = set()
    for sample_id, snps in tasks:
        for snp_pos in snps.keys():
            all_snp_pos.add(snp_pos)
    snp_pos_list = list(all_snp_pos)
    snp_pos_list.sort()
    print('total snp pos: ' + str(len(snp_pos_list)))

    cds = read_annotations(upstream_length)
    snp_to_cds = localize_all_snps(snp_pos_list, cds)

    index_list, pos_list = filter_snp_list(snp_pos_list, snp_to_cds)
    print(str(len(pos_list)) + ' snps after filtering')

    sample_to_seq = {}
    sample_sequences = SeqIO.parse(open(path_to_alignment), 'fasta')
    tasks = Parallel(n_jobs=-1)(delayed(filter_aln)(seq.name, seq.seq, index_list) for seq in sample_sequences)
    for sample_id, seq in tasks:
        sample_to_seq[sample_id] = seq
    print('alignment read')

    root, parents = read_parents(drug)

    sample_to_pheno = read_pheno(drug)

    forward, backward = format_mut_lists(sample_to_seq, pos_list, h37rv, snp_to_cds, parents)
    print('forward ' + str(len(forward)) + ' backward ' + str(len(backward)))

    # tree = Tree()
    # tree.name = root
    # node_name_to_node = {root: tree}
    # # if root in forward:
    # #     tree.name += forward_tag
    # # if root in backward:
    # #     tree.name += backward_tag
    # for node_id, parent_id, dist in parents:
    #     parent = node_name_to_node[parent_id]
    #     node = Tree()
    #     node_name_to_node[node_id] = node
    #     node.name = node_id
    #     # if node_id in forward:
    #     #     node.name += forward_tag
    #     # if node_id in backward:
    #     #     node.name += backward_tag
    #     parent.add_child(node, node_id, dist)
    # if prune:
    #     selected_ids = random.sample(sample_to_pheno.keys(), floor(random_fraction*len(sample_to_pheno)))
    #     tree.prune(selected_ids)
    # newick = tree.write(format=1)
    # for sample_id, pheno in sample_to_pheno.items():
    #     if sample_id in forward:
    #         newick = newick.replace(sample_id, sample_id + '_' + pheno + forward_tag)
    #     elif sample_id in backward:
    #         newick = newick.replace(sample_id, sample_id + '_' + pheno + backward_tag)
    #     else:
    #         newick = newick.replace(sample_id, sample_id + '_' + pheno)
    # with open(out_path + drug + '.nex', 'w') as f:
    #     f.write('#NEXUS\n')
    #     f.write('begin taxa;\n')
    #     f.write('\tdimensions ntax=' + str(len(sample_to_pheno)) + ';\n')
    #     f.write('\ttaxlabels\n')
    #     for sample_id, pheno in sample_to_pheno.items():
    #         f.write('\t' + sample_id + '_' + pheno + '\n')
    #     f.write(';\n')
    #     f.write('end;\n\n')
    #     f.write('begin trees;\n')
    #     f.write('\ttree tree_1 = [&R] ' + newick)
    #     f.write('\nend;\n\n')
    #     for line in open(path_to_fig_tree, 'r').readlines():
    #         f.write(line)
    # # tree.write(format=1, outfile=out_path + drug + '.nw')


if __name__ == '__main__':
    main()
