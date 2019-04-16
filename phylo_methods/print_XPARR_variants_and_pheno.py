from sklearn.externals.joblib import Parallel, delayed

from src.core.annotations import CDSType
from src.core.constants import data_path
from src.core.data_reading import read_h37rv
from src.phylo_methods.print_XPARR import get_aminoacids_sense, get_aminoacids_antisense
from src.phylo_methods.print_XPARR_indel import filter_all_variants, filter_alignments, read_all_variants

path_to_ids = data_path + 'dr_covered_with_pheno_and_snp.txt'
path_to_snps = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_mc10/'
path_to_alignment = data_path + 'ancestors_mc10_mega_merged.fasta'
path_to_dictionaries = data_path + 'dictionaries/'
path_to_snps_list = data_path + 'snp_aln_with_DR_with_pheno_and_snp_mc10_old_rep.txt'

first_line = False

if first_line:
    path_to_pheno = data_path + 'pheno_mc5_mega_first_line/'
    path_to_pheno_and_trees = data_path + 'reconstructed_mc10_mega_MP_first_line/'
else:
    path_to_pheno = data_path + 'pheno_mc5_mega_mix/'
    path_to_pheno_and_trees = data_path + 'reconstructed_mc10_mega_MP_mix/'

if first_line:
    drug_names = ('Isoniazid', 'Rifampicin', 'Ethambutol', 'Pyrazinamide', 'Streptomycin')
else:
    drug_names = ('Isoniazid', 'Rifampicin', 'Ethambutol', 'Pyrazinamide', 'Streptomycin', 'Moxifloxacin', 'Ofloxacin',
                  'Amikacin', 'Capreomycin', 'Kanamycin')

path_to_target_snp_pos_list = '../../res/merged_snp_pos_RR.list'
path_to_target_indel_pos_list = '../../res/merged_indel_pos_RR.list'

out_path = data_path + 'xparr/mc10_mega_MP_RR_vars_vs_vars_filtered_frd0.1_pvalue0.05_pheno_array_10drugs_fixed.xparr'
path_to_drug_codes = data_path + 'xparr/mc10_mega_MP/drug_codes.txt'
overwrite = True

print_RR = True
print_SS = False
use_DR_genes_only = False


def read_parents(path_to_pheno_and_trees, drug):
    parents = []
    with open(path_to_pheno_and_trees + drug + '/parents.csv', 'r') as f:
        root = f.readline().strip()
        for line in f.readlines():
            s = line.strip().split('\t')
            parents.append((s[0], s[1], s[2]))
    return root, parents


def format_variant(sample_id, sample_snp_seq, parent_snp_seq, snp_pos_list, ref_seq, snp_to_cds, sample_indel_seq,
                   parent_indel_seq, filtered_snp_pos_set, filtered_indel_pos_set):
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
                        if (pos - nucleotide_pos) in filtered_snp_pos_set:
                            nonsyn.append(aa0 + str(pos - nucleotide_pos) + aa1)
                    else:
                        syn.append(str(pos - nucleotide_pos))
                else:  # if strand is '-'
                    nucleotide_pos = (cds.end - pos) % 3

                    aa0, aa1, j = get_aminoacids_antisense(sample_snp_seq, parent_snp_seq, ref_seq, nucleotide_pos, snp_pos_list, j)

                    if aa0 != aa1:
                        if (pos - 2 + nucleotide_pos) in filtered_snp_pos_set:
                            nonsyn.append(aa0 + str(pos - 2 + nucleotide_pos) + aa1)
                    else:
                        syn.append(str(pos - 2 + nucleotide_pos))
            else:
                if pos in filtered_snp_pos_set:
                    nonsyn.append(parent_snp_seq[j] + str(pos) + sample_snp_seq[j])
        else:
            if pos in filtered_snp_pos_set:
                nonsyn.append(parent_snp_seq[j] + str(pos) + sample_snp_seq[j])
        j += 1
    ref_len = len(ref_seq)
    for j in range(0, len(sample_indel_seq)):
        if sample_indel_seq[j] == parent_indel_seq[j]:
            continue
        else:
            i = ref_len + j + 1
            if i in filtered_indel_pos_set:
                nonsyn.append(str(i))
    return sample_id, ';'.join(syn), ';'.join(nonsyn)


def format_mut_lists(drug, sample_to_snp_seq, snp_pos_list, ref_seq, snp_to_cds, parents, sample_to_indel_seq,
                     filtered_snp_pos_set, filtered_indel_pos_set):
    formatted_snps = Parallel(n_jobs=-1)(
        delayed(format_variant)(node_id, sample_to_snp_seq[node_id], sample_to_snp_seq[parent_id], snp_pos_list, ref_seq,
                                snp_to_cds, sample_to_indel_seq[node_id], sample_to_indel_seq[parent_id],
                                filtered_snp_pos_set, filtered_indel_pos_set)
        for node_id, parent_id, dist in parents
    )
    print('done with variant format for ' + drug)

    sample_to_mut = {}
    for sample_id, syn, nonsyn in formatted_snps:
        sample_to_mut[sample_id] = (syn, nonsyn)
    return sample_to_mut


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


def main():

    ref_seq = read_h37rv()

    sample_ids = [sample_id.strip() for sample_id in open(path_to_ids, 'r').readlines()]

    snp_pos_list_from_samples, all_snp_pos_set, all_indels_from_samples, all_indels_set_from_samples = \
        read_all_variants(sample_ids)

    snp_pos_list_from_alignment = [int(l.strip()) for l in open(path_to_snps_list, 'r').readlines()]
    indel_list_from_alignment = [l.strip() for l in open(path_to_pheno_and_trees + 'indel_list.txt')]

    cds_list, snp_to_cds, snp_index_list, pos_list, indel_index_list, filtered_indel_list = \
        filter_all_variants(snp_pos_list_from_alignment, all_snp_pos_set, indel_list_from_alignment,
                        all_indels_set_from_samples, use_DR_genes_only)

    filtered_snp_pos_set = set(int(line.strip()) for line in open(path_to_target_snp_pos_list).readlines())
    filtered_indel_pos_set = set(int(line.strip()) for line in open(path_to_target_indel_pos_list).readlines())

    sample_to_snp_seq, sample_to_indel_seq = filter_alignments(snp_index_list, indel_index_list)

    drug_to_number = {}
    for line in open(path_to_drug_codes).readlines():
        s = line.strip().split('\t')
        drug_to_number[s[0]] = s[1]
    drug_to_pheno = {}
    for drug in drug_names:
        sample_to_pheno = read_pheno(path_to_pheno, path_to_pheno_and_trees, drug)
        drug_to_pheno[drug] = sample_to_pheno

    root, parents = read_parents(path_to_pheno_and_trees, drug_names[0])
    sample_to_mut = format_mut_lists('all', sample_to_snp_seq, pos_list, ref_seq, snp_to_cds, parents,
                                     sample_to_indel_seq, filtered_snp_pos_set, filtered_indel_pos_set)

    with open(out_path, 'w') as f:
        f.write("child\tparent\tlength\tsyn\tnonsyn\tnonsyn\tphenotypes:")
        f.write(','.join([drug_to_number[drug] for drug in drug_names]) + '\n')
        f.write(root + '\n')
        for node_id, parent_id, dist in parents:
            branch_line = [node_id, parent_id, dist]
            syn, nonsyn = sample_to_mut[node_id]
            branch_line.append(syn)
            branch_line.append(nonsyn)
            branch_line.append(nonsyn)
            pheno_array = []
            for drug in drug_names:
                sample_to_pheno = drug_to_pheno[drug]
                pheno = sample_to_pheno[node_id]
                if pheno == 'R':
                    pheno_array.append('1')
                else:
                    pheno_array.append('0')
            branch_line.append(''.join(pheno_array))
            f.write('\t'.join(branch_line))
            f.write('\n')


if __name__ == '__main__':
    main()
