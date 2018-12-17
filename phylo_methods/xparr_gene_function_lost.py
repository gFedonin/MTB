from bisect import bisect_left, bisect_right
from ete3 import Tree
from Bio import SeqIO
from os.path import exists
from sklearn.externals.joblib import Parallel, delayed
import os

from src.core.annotations import CDSType, read_annotations
from src.core.constants import codon_table, codon_table_compl, upstream_length
from src.core.data_reading import read_dict, read_h37rv

path_to_ids = '../../data/dr_covered_with_pheno_and_snp.txt'
path_to_snps = '../../data/snps/raw_with_DR_with_pheno_and_snp_mc5/'
path_to_alignment = '../../data/ancestors_mc5_merged.fasta'
path_to_pheno = '../../data/pheno_mc5_Walker/'
path_to_pheno_and_trees = '../../data/reconstructed_mc5_Walker_ARD/'
path_to_dictionaries = '../../data/dictionaries/'

out_path = '../../data/xparr_mc5_Walker_ARD_all_genes/'
overwrite = True



def get_aminoacids_sense(child_seq, parent_seq, ref_seq, nucleotide_pos, snps, i):
    snp_num = len(snps)
    pos = snps[i]
    alt = child_seq[i]
    par = parent_seq[i]
    if nucleotide_pos == 0:
        if i != snp_num - 1:
            pos1 = snps[i + 1]
            alt1 = child_seq[i + 1]
            par1 = parent_seq[i + 1]
            if pos1 == pos + 1:
                if i != snp_num - 2 and snps[i + 2] == pos + 2:
                    aa0 = codon_table[par + par1 + parent_seq[i + 2]]
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
            aa0 = codon_table[ref_seq[pos - 2] + par + parent_seq[i + 1]]
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
    par = parent_seq[i]
    if nucleotide_pos == 2:
        if i != snp_num - 1:
            pos1 = snps[i + 1]
            alt1 = child_seq[i + 1]
            par1 = parent_seq[i + 1]
            if pos1 == pos + 1:
                if i != snp_num - 2 and snps[i + 2] == pos + 2:
                    aa0 = codon_table_compl[par + par1 + parent_seq[i + 2]]
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
            aa0 = codon_table_compl[ref_seq[pos - 2] + par + parent_seq[i + 1]]
            aa1 = codon_table_compl[ref_seq[pos - 2] + alt + child_seq[i + 1]]
            return aa0, aa1, i + 1
        aa0 = codon_table_compl[ref_seq[pos - 2] + par + ref_seq[pos]]
        aa1 = codon_table_compl[ref_seq[pos - 2] + alt + ref_seq[pos]]
        return aa0, aa1, i
    else:
        aa0 = codon_table_compl[ref_seq[pos - 3:pos - 1] + par]
        aa1 = codon_table_compl[ref_seq[pos - 3:pos - 1] + alt]
        return aa0, aa1, i


def localize_all_snps(all_snps, cds_list, keep_genes_set=None):
    res = {}
    for cds in cds_list:
        if cds.type == CDSType.Gene:
            if keep_genes_set is not None and cds.name not in keep_genes_set:
                continue
            i = bisect_left(all_snps, cds.start)
            j = bisect_right(all_snps, cds.end, lo=i)
            for k in range(i, j):
                res[all_snps[k]] = cds
    return res


def format_snp(sample_id, fasta_seq, parent_seq, snp_pos_list, ref_seq, snp_to_cds):
    nonsyn = []
    syn = []
    j = 0
    while j < len(snp_pos_list):
        if fasta_seq[j] == parent_seq[j]:
            j += 1
            continue
        pos = snp_pos_list[j]
        cds = snp_to_cds.get(pos)
        if cds.strand == 1:
            # if cds.type == CDSType.Gene:
            #     protein_pos = (pos - cds.start) // 3 + 1
                nucleotide_pos = (pos - cds.start) % 3

                # NEIGHBOURS ANALYSIS
                aa0, aa1, j = get_aminoacids_sense(fasta_seq, parent_seq, ref_seq, nucleotide_pos, snp_pos_list, j)

                if aa0 != aa1:
                    nonsyn.append(aa0 + str(pos - nucleotide_pos) + aa1)
                else:
                    syn.append(str(pos - nucleotide_pos))
        else:  # if strand is '-'
            # if cds.type == CDSType.Gene:
            #     protein_pos = (cds.end - pos) // 3 + 1
                nucleotide_pos = (cds.end - pos) % 3

                # NEIGHBOURS ANALYSIS
                aa0, aa1, j = get_aminoacids_antisense(fasta_seq, parent_seq, ref_seq, nucleotide_pos, snp_pos_list, j)

                if aa0 != aa1:
                    nonsyn.append(aa0 + str(pos - 2 + nucleotide_pos) + aa1)
                else:
                    syn.append(str(pos - 2 + nucleotide_pos))
        j += 1
    return sample_id, ';'.join(syn), ';'.join(nonsyn)


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


def format_mut_lists(drug, sample_to_seq, pos_list, ref_seq, snp_to_cds, drug_to_gene_set, parents):
    formatted_snps = Parallel(n_jobs=-1)(
        delayed(format_snp)(node_id, sample_to_seq[node_id], sample_to_seq[parent_id], pos_list, ref_seq,
                            snp_to_cds, drug_to_gene_set[drug])
        for node_id, parent_id, dist in parents
    )
    print('done with snp format for ' + drug)

    sample_to_mut = {}
    for sample_id, syn, nonsyn in formatted_snps:
        sample_to_mut[sample_id] = (syn, nonsyn)
    return sample_to_mut


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

    drug_to_number = []
    i = 0
    for (dirpath, dirnames, filenames) in os.walk(path_to_pheno_and_trees):
        for drug in dirnames:
            drug_to_number.append(drug + '\t' + str(i))
            i += 1
            if exists(out_path + drug + '.xparr') and not overwrite:
                continue

            root, parents = read_parents(drug)

            sample_to_pheno = read_pheno(drug)

            sample_to_mut = format_mut_lists(drug, sample_to_seq, pos_list, h37rv, snp_to_cds, drug_to_gene_set, parents)

            with open(out_path + drug + '.xparr', 'w') as f:
                f.write("child\tparent\tlength\n")
                f.write(root + '\n')
                for node_id, parent_id, dist in parents:
                    branch_line = [node_id, parent_id, dist]
                    syn, nonsyn = sample_to_mut[node_id]
                    branch_line.append(syn)
                    pheno = sample_to_pheno[node_id]
                    parent_pheno = sample_to_pheno[parent_id]
                    if parent_pheno != pheno:
                        branch_line.append(parent_pheno + str(i) + pheno)
                    else:
                        branch_line.append('')
                    branch_line.append(nonsyn)
                    f.write('\t'.join(branch_line))
                    f.write('\n')
    with open(out_path + 'drug_codes.txt', 'w') as f:
        f.write('\n'.join(drug_to_number))


if __name__ == '__main__':
    main()
