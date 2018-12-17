from ete3 import Tree
from os.path import exists
from sklearn.externals.joblib import Parallel, delayed
import os

from src.core.annotations import CDSType, read_annotations, localize_all_snps
from src.core.constants import upstream_length, data_path
from src.core.data_reading import read_h37rv
from src.phylo_methods.print_XPARR import get_aminoacids_sense, get_aminoacids_antisense, read_pheno, filter_snp_list
from src.phylo_methods.print_XPARR_indel import read_all_variants, filter_alignments, read_parents, \
    filter_indel_list

path_to_ids = data_path + 'dr_covered_with_pheno_and_snp.txt'
path_to_snps = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_mc10/'
path_to_alignment = data_path + 'ancestors_mc10_mega_merged.fasta'
path_to_pheno = data_path + 'pheno_mc5_mega/'
path_to_pheno_and_trees = data_path + 'reconstructed_mc10_mega_MP/'
path_to_dictionaries = data_path + 'dictionaries/'
path_to_snps_list = data_path + 'snp_aln_with_DR_with_pheno_and_snp_mc10_old_rep.txt'

out_path = data_path + 'xparr/mc10_mega_MP_genes/'
overwrite = True


def format_variant(sample_id, sample_snp_seq, parent_snp_seq, snp_pos_list, ref_seq, snp_to_cds, sample_indel_seq,
                   parent_indel_seq, indel_list, indel_to_cds, gene_to_id):
    nonsyn = set()
    syn = set()
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
                        nonsyn.add(gene_to_id[cds.name])
                    else:
                        syn.add(gene_to_id[cds.name])
                else:  # if strand is '-'
                    nucleotide_pos = (cds.end - pos) % 3

                    aa0, aa1, j = get_aminoacids_antisense(sample_snp_seq, parent_snp_seq, ref_seq, nucleotide_pos, snp_pos_list, j)

                    if aa0 != aa1:
                        nonsyn.add(gene_to_id[cds.name])
                    else:
                        syn.add(gene_to_id[cds.name])
            else:
                nonsyn.add(gene_to_id[cds.name])
        j += 1
    for j in range(0, len(sample_indel_seq)):
        if sample_indel_seq[j] == parent_indel_seq[j]:
            continue
        else:
            s = indel_list[j].split('_')
            cds = indel_to_cds.get(int(s[1]))
            if cds is not None:
                nonsyn.add(gene_to_id[cds.name])
    syn = list(syn)
    syn.sort()
    nonsyn = list(nonsyn)
    nonsyn.sort()
    return sample_id, ';'.join(syn), ';'.join(nonsyn)


def format_mut_lists(drug, sample_to_snp_seq, snp_pos_list, ref_seq, snp_to_cds, parents, sample_to_indel_seq, indel_list, indel_to_cds, gene_to_id):
    formatted_snps = Parallel(n_jobs=-1)(
        delayed(format_variant)(node_id, sample_to_snp_seq[node_id], sample_to_snp_seq[parent_id], snp_pos_list, ref_seq,
                                snp_to_cds, sample_to_indel_seq[node_id], sample_to_indel_seq[parent_id], indel_list, indel_to_cds, gene_to_id)
        for node_id, parent_id, dist in parents
    )
    print('done with variant format for ' + drug)

    sample_to_mut = {}
    for sample_id, syn, nonsyn in formatted_snps:
        sample_to_mut[sample_id] = (syn, nonsyn)
    return sample_to_mut


def filter_all_variants(snp_pos_list_from_alignment, snp_pos_set_from_samples, indel_list_from_alignment,
                        all_indels_set_from_samples):
    cds_list = read_annotations(upstream_length)
    snp_to_cds = localize_all_snps(snp_pos_list_from_alignment, cds_list)
    indel_pos_set = set()
    for m in indel_list_from_alignment:
        indel_pos_set.add(int(m.split('_')[1]))
    indel_pos_list = list(indel_pos_set)
    indel_to_cds = localize_all_snps(indel_pos_list, cds_list)
    snp_index_list, pos_list = filter_snp_list(snp_pos_list_from_alignment, snp_pos_set_from_samples)
    indel_index_list, filtered_indel_list = filter_indel_list(indel_list_from_alignment, all_indels_set_from_samples,
                                                              indel_to_cds)
    print(str(len(pos_list)) + ' snps after filtering')
    print(str(len(filtered_indel_list)) + ' indels after filtering')
    return cds_list, snp_to_cds, snp_index_list, pos_list, indel_index_list, filtered_indel_list, indel_to_cds


def main():

    if not exists(out_path):
        os.makedirs(out_path)

    ref_seq = read_h37rv()

    sample_ids = [sample_id.strip() for sample_id in open(path_to_ids, 'r').readlines()]

    snp_pos_list_from_samples, all_snp_pos_set, all_indels_from_samples, all_indels_set_from_samples = \
        read_all_variants(sample_ids)

    snp_pos_list_from_alignment = [int(l.strip()) for l in open(path_to_snps_list, 'r').readlines()]
    indel_list_from_alignment = [l.strip() for l in open(path_to_pheno_and_trees + 'indel_list.txt')]

    cds_list, snp_to_cds, snp_index_list, pos_list, indel_index_list, filtered_indel_list, indel_to_cds = \
        filter_all_variants(snp_pos_list_from_alignment, all_snp_pos_set, indel_list_from_alignment,
                        all_indels_set_from_samples)

    gene_to_id = {}
    gene_list = []
    for cds in cds_list:
        if cds.type != CDSType.upstream:
            gene_list.append(cds.name)
    for i in range(len(gene_list)):
        gene_to_id[gene_list[i]] = str(i)

    with open(out_path + 'gene_list.txt', 'w') as f:
        f.write('\n'.join(gene_list))
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
                                             sample_to_indel_seq, filtered_indel_list, indel_to_cds, gene_to_id)

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
