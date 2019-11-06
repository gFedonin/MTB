from Bio import SeqIO
from os.path import exists
from sklearn.externals.joblib import Parallel, delayed
import os

from core.annotations import read_annotations, localize_all_variants, CDSType
from core.constants import codon_table, codon_table_compl, upstream_length, data_path
from core.data_reading import read_dict, read_h37rv

path_to_ids = data_path + 'dr_covered_with_pheno_and_snp.txt'
path_to_snps = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_mc10/'
path_to_alignment = data_path + 'ancestors_mc10_mega_merged.fasta'
path_to_snps_list = data_path + 'snp_aln_with_DR_with_pheno_and_snp_mc10_old_rep.txt'
path_to_pheno = data_path + 'pheno_mc5_mega/'
path_to_pheno_and_trees = data_path + 'reconstructed_mc10_mega_MP/'
path_to_dictionaries = data_path + 'dictionaries/'

out_path = data_path + 'xparr/mc10_mega_snp_dr_only/'


overwrite = True
use_DR_genes_only = True
thread_num = 32


def extract_dr_genes(path_to_dictionaries):
    name, drug_to_mut_list = read_dict(path_to_dictionaries, 'Walker_dictionary')
    drug_to_gene_set = {}
    all_genes = set()
    for drug, mut_list in drug_to_mut_list.items():
        genes = set()
        for mut in mut_list:
            s = mut.split('\t')
            genes.add(s[1])
            all_genes.add(s[1])
        drug_to_gene_set[drug] = genes
    return drug_to_gene_set, all_genes


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


def format_snp(sample_seq, parent_seq, snp_pos_list, ref_seq, snp_to_cds):
    nonsyn = []
    syn = []
    j = 0
    while j < len(snp_pos_list):
        if sample_seq[j] == parent_seq[j]:
            j += 1
            continue
        pos = snp_pos_list[j]
        cds = snp_to_cds.get(pos)

        if cds is not None:
            if cds.type == CDSType.Gene:
                if cds.strand == 1:
                #     protein_pos = (pos - cds.start) // 3 + 1
                    nucleotide_pos = (pos - cds.start) % 3

                    # NEIGHBOURS ANALYSIS
                    aa0, aa1, j = get_aminoacids_sense(sample_seq, parent_seq, ref_seq, nucleotide_pos, snp_pos_list, j)

                    if aa0 != aa1:
                        nonsyn.append(aa0 + str(pos - nucleotide_pos) + aa1)
                    else:
                        syn.append(str(pos - nucleotide_pos))
                else:  # if strand is '-'
                # if cds.type == CDSType.Gene:
                #     protein_pos = (cds.end - pos) // 3 + 1
                    nucleotide_pos = (cds.end - pos) % 3

                    # NEIGHBOURS ANALYSIS
                    aa0, aa1, j = get_aminoacids_antisense(sample_seq, parent_seq, ref_seq, nucleotide_pos, snp_pos_list, j)

                    if aa0 != aa1:
                        nonsyn.append(aa0 + str(pos - 2 + nucleotide_pos) + aa1)
                    else:
                        syn.append(str(pos - 2 + nucleotide_pos))
            else:
                nonsyn.append(parent_seq[j] + str(pos) + sample_seq[j])
        else:
            nonsyn.append(parent_seq[j] + str(pos) + sample_seq[j])
        j += 1
    return ';'.join(syn), ';'.join(nonsyn)


def read_snps(path_to_snps, sample_id):
    snps = {}
    with open(path_to_snps + sample_id + '.variants', 'r') as f1:
        for line in f1.readlines():
            s = line.strip().split('\t')
            if s[2] == 'snp':
                snps[int(s[0])] = s[1]
    return sample_id, snps


def filter_snp_list(snp_pos_list, all_snp_pos, snp_to_cds=None):
    index_list = []
    pos_list = []
    for i in range(len(snp_pos_list)):
        pos = snp_pos_list[i]
        if pos in all_snp_pos:
            if snp_to_cds is not None:
                if snp_to_cds.get(pos) is not None:
                    index_list.append(i)
                    pos_list.append(pos)
            else:
                index_list.append(i)
                pos_list.append(pos)
    return index_list, pos_list


def read_parents(path_to_pheno_and_trees, drug):
    parents = []
    with open(path_to_pheno_and_trees + drug + '/parents.csv', 'r') as f:
        root = f.readline().strip()
        for line in f.readlines():
            s = line.strip().split('\t')
            parents.append((s[0], s[1], s[2]))
    return root, parents


def read_pheno(path_to_pheno, path_to_pheno_and_trees, drug):

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


def format_mut_for_samples(sample_pair_list, sample_to_snp_seq, snp_pos_list, ref_seq, snp_to_cds):
    res = []
    for node_id, parent_id in sample_pair_list:
        syn, nonsyn = format_snp(sample_to_snp_seq[node_id], sample_to_snp_seq[parent_id], snp_pos_list, ref_seq,
                                 snp_to_cds)
        res.append((node_id, syn, nonsyn))
    return res


def format_mut_lists(drug, sample_to_snp_seq, snp_pos_list, ref_seq, snp_to_cds, parents):
    sample_pair_lists = []
    sample_to_snp_seq_list = []
    for i in range(thread_num):
        sample_pair_lists.append([])
        sample_to_snp_seq_list.append({})
    for i in range(len(parents)):
        node_id, parent_id, dist = parents[i]
        thread_index = i % thread_num
        sample_pair_lists[thread_index].append((node_id, parent_id))
        sample_to_snp_seq_list[thread_index][node_id] = sample_to_snp_seq[node_id]
        sample_to_snp_seq_list[thread_index][parent_id] = sample_to_snp_seq[parent_id]

    formatted_snps = Parallel(n_jobs=thread_num)(
        delayed(format_mut_for_samples)(sample_pair_lists[i], sample_to_snp_seq_list[i], snp_pos_list, ref_seq, snp_to_cds)
        for i in range(len(sample_pair_lists))
    )
    print('done with snp format for ' + drug)

    sample_to_mut = {}
    for l in formatted_snps:
        for sample_id, syn, nonsyn in l:
            sample_to_mut[sample_id] = (syn, nonsyn)
    return sample_to_mut


def read_all_snps(path_to_snps, sample_ids):
    tasks = Parallel(n_jobs=thread_num)(delayed(read_snps)(path_to_snps, sample_id) for sample_id in sample_ids)
    all_snp_pos_set = set()
    for sample_id, snps in tasks:
        for snp_pos in snps.keys():
            all_snp_pos_set.add(snp_pos)
    snp_pos_list = list(all_snp_pos_set)
    print('total snp pos: ' + str(len(snp_pos_list)))
    return snp_pos_list, all_snp_pos_set


def filter_aln_list(aln_list, index_list):
    res = []
    for seq in aln_list:
        filtered_seq = []
        for i in index_list:
            filtered_seq.append(seq[i])
        res.append((seq.name, ''.join(filtered_seq)))
    return res


def filter_all_aln(sample_sequences, index_list):
    sample_to_seq = {}
    aln_lists = []
    for i in range(thread_num):
        aln_lists.append([])
    i = 0
    for seq in sample_sequences:
        aln_lists[i % thread_num].append(seq)
        i += 1
    tasks = Parallel(n_jobs=thread_num)(delayed(filter_aln_list)(aln_list, index_list) for aln_list in aln_lists)
    for l in tasks:
        for sample_id, seq in l:
            sample_to_seq[sample_id] = seq
    return sample_to_seq


def main():

    if not exists(out_path):
        os.makedirs(out_path)

    ref_seq = read_h37rv()

    sample_ids = [sample_id.strip() for sample_id in open(path_to_ids, 'r').readlines()]
    snp_pos_list_from_samples, all_snp_pos_set = read_all_snps(path_to_snps, sample_ids)

    cds_list = read_annotations(upstream_length)

    if use_DR_genes_only:
        drug_to_gene_set, all_dr_genes = extract_dr_genes(path_to_dictionaries)
        snp_to_cds = localize_all_variants(snp_pos_list_from_samples, cds_list, all_dr_genes)
    else:
        snp_to_cds = localize_all_variants(snp_pos_list_from_samples, cds_list)

    snp_pos_list_from_alignment = [int(l.strip()) for l in open(path_to_snps_list, 'r').readlines()]
    index_list, pos_list = filter_snp_list(snp_pos_list_from_alignment, all_snp_pos_set, snp_to_cds)
    print(str(len(index_list)) + ' snps after filtering')


    sample_sequences = SeqIO.parse(open(path_to_alignment, 'r'), 'fasta')
    sample_to_seq = filter_all_aln(sample_sequences, index_list)

    drug_to_number = []
    i = 0
    for (dirpath, dirnames, filenames) in os.walk(path_to_pheno_and_trees):
        for drug in dirnames:
            drug_to_number.append(drug + '\t' + str(i))
            i += 1
            if exists(out_path + drug + '.xparr') and not overwrite:
                continue

            root, parents = read_parents(path_to_pheno_and_trees, drug)

            sample_to_pheno = read_pheno(path_to_pheno, path_to_pheno_and_trees, drug)

            sample_to_mut = format_mut_lists(drug, sample_to_seq, pos_list, ref_seq, snp_to_cds, parents)

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
                        branch_line.append(parent_pheno + str(i - 1) + pheno)
                    else:
                        branch_line.append('')
                    branch_line.append(nonsyn)
                    f.write('\t'.join(branch_line))
                    f.write('\n')
    with open(out_path + 'drug_codes.txt', 'w') as f:
        f.write('\n'.join(drug_to_number))


if __name__ == '__main__':
    main()
