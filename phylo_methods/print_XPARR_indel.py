from Bio import SeqIO
from os.path import exists
from sklearn.externals.joblib import Parallel, delayed
import os

from src.core.annotations import CDSType, read_annotations, localize_all_snps
from src.core.constants import upstream_length, data_path
from src.core.data_reading import read_dict, read_h37rv
from src.phylo_methods.print_XPARR import get_aminoacids_sense, get_aminoacids_antisense, read_pheno, read_parents, \
    filter_snp_list, filter_all_aln

path_to_ids = data_path + 'dr_covered_with_pheno_and_snp.txt'
path_to_snps = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_mc10/'
path_to_alignment = data_path + 'ancestors_mc10_mega_merged.fasta'
path_to_pheno = data_path + 'pheno_mc5_mega/'
path_to_pheno_and_trees = data_path + 'reconstructed_mc10_mega_MP/'
path_to_dictionaries = data_path + 'dictionaries/'
path_to_snps_list = data_path + 'snp_aln_with_DR_with_pheno_and_snp_mc10_old_rep.txt'

out_path = data_path + 'xparr/mc10_mega_MP/'
overwrite = True

use_DR_genes_only = False


def extract_dr_genes():
    name, drug_to_mut_list = read_dict(path_to_dictionaries, 'Walker_dictionary')
    drug_to_gene_set = {}
    all_genes = set()
    for drug, mut_list in drug_to_mut_list.items():
        genes = set()
        for mut in mut_list:
            s = mut.split('\t')
            if s[0] == 'Gene':
                genes.add(s[1])
                all_genes.add(s[1])
        drug_to_gene_set[drug] = genes
    return drug_to_gene_set, all_genes


def format_variant(sample_id, sample_snp_seq, parent_snp_seq, snp_pos_list, ref_seq, snp_to_cds, sample_indel_seq,
                   parent_indel_seq):
    nonsyn = []
    syn = []
    j = 0
    while j < len(snp_pos_list):
        if sample_snp_seq[j] == parent_snp_seq[j]:
            j += 1
            continue
        pos = snp_pos_list[j]
        cds = snp_to_cds.get(pos)
        if cds is not None:
            if cds.type == CDSType.Gene:
                if cds.strand == 1:
                    nucleotide_pos = (pos - cds.start) % 3

                    aa0, aa1, j = get_aminoacids_sense(sample_snp_seq, parent_snp_seq, ref_seq, nucleotide_pos, snp_pos_list, j)

                    if aa0 != aa1:
                        nonsyn.append(aa0 + str(pos - nucleotide_pos) + aa1)
                    else:
                        syn.append(str(pos - nucleotide_pos))
                else:  # if strand is '-'
                    nucleotide_pos = (cds.end - pos) % 3

                    aa0, aa1, j = get_aminoacids_antisense(sample_snp_seq, parent_snp_seq, ref_seq, nucleotide_pos, snp_pos_list, j)

                    if aa0 != aa1:
                        nonsyn.append(aa0 + str(pos - 2 + nucleotide_pos) + aa1)
                    else:
                        syn.append(str(pos - 2 + nucleotide_pos))
            else:
                nonsyn.append(parent_snp_seq[j] + str(pos) + sample_snp_seq[j])
        else:
            nonsyn.append(parent_snp_seq[j] + str(pos) + sample_snp_seq[j])
        j += 1
    ref_len = len(ref_seq)
    for j in range(0, len(sample_indel_seq)):
        if sample_indel_seq[j] == parent_indel_seq[j]:
            continue
        else:
            nonsyn.append(str(ref_len + j + 1))
    return sample_id, ';'.join(syn), ';'.join(nonsyn)


def read_variants(sample_id):
    snps = []
    inserts = []
    deletions = []
    with open(path_to_snps + sample_id + '.variants', 'r') as f1:
        del_start = 0
        del_len = 0
        for line in f1.readlines():
            s = line.strip().split('\t')
            pos = int(s[0])
            if s[2] == 'del':
                if del_len == 0:
                    del_len = 1
                    del_start = pos
                else:
                    if del_start + del_len == pos:
                        del_len += 1
                    else:
                        deletions.append((del_start, del_len))
                        del_start = pos
                        del_len = 1
            elif s[2] == 'ins':
                if del_len > 0:
                    deletions.append((del_start, del_len))
                    del_start = 0
                    del_len = 0
                inserts.append((pos, s[1]))
            else:
                if del_len > 0:
                    deletions.append((del_start, del_len))
                    del_start = 0
                    del_len = 0
                snps.append((pos, s[1]))
        if del_len > 0:
            deletions.append((del_start, del_len))
    return sample_id, snps, inserts, deletions


def read_all_variants(sample_ids):
    tasks = Parallel(n_jobs=-1)(delayed(read_variants)(sample_id) for sample_id in sample_ids)
    all_snp_pos = set()
    all_indels = set()
    for sample_id, snps, inserts, deletions in tasks:
        for snp_pos, alt in snps:
            all_snp_pos.add(snp_pos)
        for pos, alt in inserts:
            all_indels.add('ins_' + str(pos) + '_' + alt)
        for pos, length in deletions:
            all_indels.add('del_' + str(pos) + '_' + str(length))
    snp_pos_list = list(all_snp_pos)
    snp_pos_list.sort()
    # indel_list = list(all_indels)
    # indel_list.sort()
    print('total snp pos in samples: ' + str(len(snp_pos_list)))
    print('total indels in samples: ' + str(len(all_indels)))
    return snp_pos_list, all_snp_pos, all_indels, set(all_indels)


def filter_indel_list(indel_list, filter_set, indel_to_cds=None):
    index_list = []
    filtered_indel_list = []
    for i in range(len(indel_list)):
        indel = indel_list[i]
        if indel in filter_set:
            if indel_to_cds != None:
                pos = int(indel.split('_')[1])
                if indel_to_cds.get(pos) is not None:
                    index_list.append(i)
                    filtered_indel_list.append(indel)
            else:
                index_list.append(i)
                filtered_indel_list.append(indel)
    return index_list, filtered_indel_list


def filter_alignments(snp_index_list, indel_index_list):
    sample_sequences = SeqIO.parse(open(path_to_alignment, 'r'), 'fasta')
    sample_to_snp_seq = filter_all_aln(sample_sequences, snp_index_list)
    print('snp alignment read')

    sample_sequences = SeqIO.parse(open(path_to_pheno_and_trees + 'anc_indels.fasta', 'r'), 'fasta')
    sample_to_indel_seq = filter_all_aln(sample_sequences, indel_index_list)
    print('indel alignment read')

    return sample_to_snp_seq, sample_to_indel_seq


def filter_all_variants(snp_pos_list_from_alignment, snp_pos_set_from_samples, indel_list_from_alignment,
                        all_indels_set_from_samples, use_DR_genes_only=False):

    cds_list = read_annotations(upstream_length)
    if use_DR_genes_only:
        drug_to_gene_set, all_dr_genes = extract_dr_genes()
        snp_to_cds = localize_all_snps(snp_pos_list_from_alignment, cds_list, all_dr_genes)
        indel_pos_set = set()
        for m in indel_list_from_alignment:
            indel_pos_set.add(int(m.split('_')[1]))
        indel_pos_list = list(indel_pos_set)
        indel_to_cds = localize_all_snps(indel_pos_list, cds_list, all_dr_genes)
        snp_index_list, pos_list = filter_snp_list(snp_pos_list_from_alignment, snp_to_cds)
        indel_index_list, filtered_indel_list = filter_indel_list(indel_list_from_alignment, all_indels_set_from_samples, indel_to_cds)
    else:
        snp_to_cds = localize_all_snps(snp_pos_list_from_alignment, cds_list)
        snp_index_list, pos_list = filter_snp_list(snp_pos_list_from_alignment, snp_pos_set_from_samples)
        indel_index_list, filtered_indel_list = filter_indel_list(indel_list_from_alignment, all_indels_set_from_samples)
    print(str(len(pos_list)) + ' snps after filtering')
    print(str(len(filtered_indel_list)) + ' indels after filtering')
    return cds_list, snp_to_cds, snp_index_list, pos_list, indel_index_list, filtered_indel_list


def format_mut_lists(drug, sample_to_snp_seq, snp_pos_list, ref_seq, snp_to_cds, parents, sample_to_indel_seq):
    formatted_snps = Parallel(n_jobs=-1)(
        delayed(format_variant)(node_id, sample_to_snp_seq[node_id], sample_to_snp_seq[parent_id], snp_pos_list, ref_seq,
                                snp_to_cds, sample_to_indel_seq[node_id], sample_to_indel_seq[parent_id])
        for node_id, parent_id, dist in parents
    )
    print('done with variant format for ' + drug)

    sample_to_mut = {}
    for sample_id, syn, nonsyn in formatted_snps:
        sample_to_mut[sample_id] = (syn, nonsyn)
    return sample_to_mut


def main():

    if not exists(out_path):
        os.makedirs(out_path)

    ref_seq = read_h37rv()

    sample_ids = [sample_id.strip() for sample_id in open(path_to_ids, 'r').readlines()]

    snp_pos_list_from_samples, all_snp_pos_set, all_indels_from_samples, all_indels_set_from_samples = \
        read_all_variants(sample_ids)

    snp_pos_list_from_alignment = [int(l.strip()) for l in open(path_to_snps_list, 'r').readlines()]
    indel_list_from_alignment = [l.strip() for l in open(path_to_pheno_and_trees + 'indel_list.txt')]

    cds_list, snp_to_cds, snp_index_list, pos_list, indel_index_list, filtered_indel_list = \
        filter_all_variants(snp_pos_list_from_alignment, all_snp_pos_set, indel_list_from_alignment,
                        all_indels_set_from_samples, use_DR_genes_only)

    with open(out_path + 'indel_list.txt', 'w') as f:
        f.write('\n'.join(filtered_indel_list))
        f.write('\n')

    sample_to_snp_seq, sample_to_indel_seq = filter_alignments(snp_index_list, indel_index_list)

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

            sample_to_mut = format_mut_lists(drug, sample_to_snp_seq, pos_list, ref_seq, snp_to_cds, parents,
                                             sample_to_indel_seq)

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
