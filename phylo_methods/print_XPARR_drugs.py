from Bio import SeqIO
from os.path import exists
from sklearn.externals.joblib import Parallel, delayed

from src.core.annotations import CDSType, read_annotations, localize_all_snps
from src.core.constants import data_path, upstream_length
from src.core.data_reading import read_h37rv
from src.phylo_methods.print_XPARR import get_aminoacids_sense, get_aminoacids_antisense, read_all_snps, \
    filter_snp_list, filter_all_aln

path_to_pheno = data_path + 'pheno_mc5_mega_mix/'
path_to_pheno_and_trees = data_path + 'reconstructed_mc10_mega_MP_mix/'
path_to_ids = data_path + 'dr_covered_with_pheno_and_snp.txt'
path_to_snps_list = data_path + 'snp_aln_with_DR_with_pheno_and_snp_mc10_old_rep.txt'
path_to_snps = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_mc10/'
path_to_alignment = data_path + 'ancestors_mc10_mega_merged.fasta'

out_path = data_path + 'xparr/mc10_mega_drugs_vs_drugs_10drugs.xparr'

overwrite = True
thread_num = 32

# drug_names = ('Isoniazid', 'Rifampicin', 'Ethambutol', 'Pyrazinamide', 'Streptomycin')
drug_names = ('Isoniazid', 'Rifampicin', 'Ethambutol', 'Pyrazinamide', 'Streptomycin', 'Moxifloxacin', 'Ofloxacin',
              'Amikacin', 'Capreomycin', 'Kanamycin')


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


def format_variant(sample_id, sample_snp_seq, parent_snp_seq, snp_pos_list, ref_seq, snp_to_cds):
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

                    if aa0 == aa1:
                        syn.append(str(pos - nucleotide_pos))
                else:  # if strand is '-'
                    nucleotide_pos = (cds.end - pos) % 3

                    aa0, aa1, j = get_aminoacids_antisense(sample_snp_seq, parent_snp_seq, ref_seq, nucleotide_pos, snp_pos_list, j)

                    if aa0 == aa1:
                        syn.append(str(pos - 2 + nucleotide_pos))
        j += 1
    return sample_id, ';'.join(syn)


def format_mut_lists(sample_to_snp_seq, snp_pos_list, ref_seq, snp_to_cds, parents):
    formatted_snps = Parallel(n_jobs=-1)(
        delayed(format_variant)(node_id, sample_to_snp_seq[node_id], sample_to_snp_seq[parent_id], snp_pos_list, ref_seq,
                                snp_to_cds)
        for node_id, parent_id, dist in parents
    )

    sample_to_mut = {}
    for sample_id, syn in formatted_snps:
        sample_to_mut[sample_id] = syn
    return sample_to_mut


def main():

    ref_seq = read_h37rv()

    sample_ids = [sample_id.strip() for sample_id in open(path_to_ids, 'r').readlines()]
    snp_pos_list_from_samples, all_snp_pos_set = read_all_snps(path_to_snps, sample_ids)

    cds_list = read_annotations(upstream_length)

    snp_to_cds = localize_all_snps(snp_pos_list_from_samples, cds_list)

    snp_pos_list_from_alignment = [int(l.strip()) for l in open(path_to_snps_list, 'r').readlines()]
    index_list, pos_list = filter_snp_list(snp_pos_list_from_alignment, all_snp_pos_set, snp_to_cds)
    print(str(len(index_list)) + ' snps after filtering')

    sample_sequences = SeqIO.parse(open(path_to_alignment, 'r'), 'fasta')
    sample_to_seq = filter_all_aln(sample_sequences, index_list)

    drug_to_number = []
    i = 0

    drug_to_pheno = {}
    drug_to_parents = {}
    for drug in drug_names:
        drug_to_number.append(drug + '\t' + str(i))
        i += 1
        sample_to_pheno = read_pheno(path_to_pheno, path_to_pheno_and_trees, drug)
        drug_to_pheno[drug] = sample_to_pheno
        root, parents = read_parents(path_to_pheno_and_trees, drug)
        drug_to_parents[drug] = (root, parents)

    sample_to_mut = format_mut_lists(sample_to_seq, pos_list, ref_seq, snp_to_cds, drug_to_parents[drug_names[0]][1])

    if exists(out_path) and not overwrite:
        return

    with open(out_path, 'w') as f:
        f.write("child\tparent\tlength\n")
        root, parents = drug_to_parents[drug_names[0]]
        f.write(root + '\n')
        for node_id, parent_id, dist in parents:
            branch_line = [node_id, parent_id, dist]
            syn = sample_to_mut[node_id]
            branch_line.append(syn)
            nonsyn = []
            for i in range(len(drug_names)):
                sample_to_pheno = drug_to_pheno[drug_names[i]]
                pheno = sample_to_pheno[node_id]
                parent_pheno = sample_to_pheno[parent_id]
                if parent_pheno != pheno:
                    nonsyn.append(parent_pheno + str(i) + pheno)
            branch_line.append(';'.join(nonsyn))
            branch_line.append(';'.join(nonsyn))
            f.write('\t'.join(branch_line))
            f.write('\n')
    with open(out_path + '.drug_codes', 'w') as f:
        f.write('\n'.join(drug_to_number))


if __name__ == '__main__':
    main()