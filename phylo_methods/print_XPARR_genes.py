from ete3 import Tree
from Bio import SeqIO
from os.path import exists
from sklearn.externals.joblib import Parallel, delayed
import os

from src.core.annotations import read_annotations, localize_all_variants, CDSType
from src.core.constants import upstream_length, data_path
from src.core.data_reading import read_h37rv
from src.phylo_methods.print_XPARR import read_all_snps, filter_snp_list, get_aminoacids_sense, \
    get_aminoacids_antisense, read_parents, read_pheno, filter_all_aln

path_to_ids = data_path + 'dr_covered_with_pheno_and_snp.txt'
path_to_snps = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_mc10/'
path_to_alignment = data_path + 'ancestors_mc10_mega_merged.fasta'
path_to_snps_list = data_path + 'snp_aln_with_DR_with_pheno_and_snp_mc10_old_rep.txt'
path_to_pheno = data_path + 'pheno_mc5_mega/'
path_to_pheno_and_trees = data_path + 'reconstructed_mc10_mega_MP/'

out_path = data_path + 'xparr/mc10_mega_snp_genes/'

overwrite = True
thread_num = 32


def format_snp(sample_seq, parent_seq, snp_pos_list, ref_seq, snp_to_cds, gene_to_id):
    nonsyn = set()
    syn = set()
    j = 0
    while j < len(snp_pos_list):
        if sample_seq[j] == parent_seq[j]:
            j += 1
            continue
        pos = snp_pos_list[j]
        cds = snp_to_cds.get(pos)
        if cds is None:
            j += 1
            continue

        if cds is not None:
            if cds.type == CDSType.Gene:
                if cds.strand == 1:
                        nucleotide_pos = (pos - cds.start) % 3

                        aa0, aa1, j = get_aminoacids_sense(sample_seq, parent_seq, ref_seq, nucleotide_pos, snp_pos_list, j)

                        if aa0 != aa1:
                            nonsyn.add(gene_to_id[cds.name])
                        else:
                            syn.add(gene_to_id[cds.name])
                else:  # if strand is '-'
                        nucleotide_pos = (cds.end - pos) % 3

                        aa0, aa1, j = get_aminoacids_antisense(sample_seq, parent_seq, ref_seq, nucleotide_pos, snp_pos_list, j)

                        if aa0 != aa1:
                            nonsyn.add(gene_to_id[cds.name])
                        else:
                            syn.add(gene_to_id[cds.name])
            else:
                nonsyn.add(gene_to_id[cds.name])
        j += 1
    syn = list(syn)
    syn.sort()
    nonsyn = list(nonsyn)
    nonsyn.sort()
    return ';'.join(syn), ';'.join(nonsyn)


def format_mut_for_samples(sample_pair_list, sample_to_snp_seq, snp_pos_list, ref_seq, snp_to_cds, gene_to_id):
    res = []
    for node_id, parent_id in sample_pair_list:
        syn, nonsyn = format_snp(sample_to_snp_seq[node_id], sample_to_snp_seq[parent_id], snp_pos_list, ref_seq,
                                 snp_to_cds, gene_to_id)
        res.append((node_id, syn, nonsyn))
    return res


def format_mut_lists(drug, sample_to_snp_seq, snp_pos_list, ref_seq, snp_to_cds, parents, gene_to_id):
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
        delayed(format_mut_for_samples)(sample_pair_lists[i], sample_to_snp_seq_list[i], snp_pos_list, ref_seq,
                                        snp_to_cds, gene_to_id)
        for i in range(len(sample_pair_lists))
    )
    print('done with snp format for ' + drug)

    sample_to_mut = {}
    for l in formatted_snps:
        for sample_id, syn, nonsyn in l:
            sample_to_mut[sample_id] = (syn, nonsyn)
    return sample_to_mut


def main():

    if not exists(out_path):
        os.makedirs(out_path)

    ref_seq = read_h37rv()

    sample_ids = [sample_id.strip() for sample_id in open(path_to_ids, 'r').readlines()]
    snp_pos_list, all_snp_pos = read_all_snps(path_to_snps, sample_ids)

    cds_list = read_annotations(upstream_length)

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

    snp_to_cds = localize_all_variants(snp_pos_list, cds_list)

    index_list, pos_list = filter_snp_list(path_to_snps_list, all_snp_pos, snp_to_cds)

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

            sample_to_mut = format_mut_lists(drug, sample_to_seq, pos_list, ref_seq, snp_to_cds, parents, gene_to_id)

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
