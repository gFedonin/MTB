from Bio import SeqIO
from os.path import exists
from sklearn.externals.joblib import Parallel, delayed
from ete3 import Tree
from src.core.annotations import CDSType, read_annotations, localize_all_variants
from src.core.constants import data_path, upstream_length
from src.core.data_reading import read_h37rv
from src.phylo_methods.print_XPARR import get_aminoacids_sense, get_aminoacids_antisense, read_all_snps, \
    filter_snp_list, filter_all_aln

first_line = False

# if first_line:
#     path_to_pheno = data_path + 'pheno_mc5_mega_first_line/'
#     path_to_pheno_and_trees = data_path + 'reconstructed_mc10_mega_MP_first_line/'
# else:
#     path_to_pheno = data_path + 'pheno_mc5_mega_mix/'
#     path_to_pheno_and_trees = data_path + 'reconstructed_mc10_mega_MP_mix/'
path_to_pheno = data_path + 'pheno_mc5_mega/'
path_to_pheno_and_trees = data_path + 'reconstructed_mc10_mega_MP/'

# path_to_ids = data_path + 'dr_covered_with_pheno_and_snp.txt'
path_to_ids = data_path + '10drugs.sample_list'
path_to_snps_list = data_path + 'snp_aln_with_DR_with_pheno_and_snp_mc10_old_rep.txt'
path_to_snps = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_mc10/'
path_to_alignment = data_path + 'ancestors_mc10_mega_merged.fasta'
path_to_drug_codes = data_path + 'xparr/mc10_mega_MP/drug_codes.txt'

print_RR = False

if first_line:
    if print_RR:
        out_path = data_path + 'xparr/mc10_mega_drugs_vs_drugs_RR_first_line_fixed.xparr'
    else:
        out_path = data_path + 'xparr/mc10_mega_drugs_vs_drugs_first_line_fixed.xparr'
else:
    if print_RR:
        out_path = data_path + 'xparr/mc10_mega_drugs_vs_drugs_RR_10drugs_fixed.xparr'
    else:
        out_path = data_path + 'xparr/mc10_mega_drugs_vs_drugs_10drugs_fixed.xparr'

overwrite = True
thread_num = 32


if first_line:
    drug_names = ('Isoniazid', 'Rifampicin', 'Ethambutol', 'Pyrazinamide', 'Streptomycin')
else:
    drug_names = ('Isoniazid', 'Rifampicin', 'Ethambutol', 'Pyrazinamide', 'Streptomycin', 'Moxifloxacin', 'Ofloxacin',
                  'Amikacin', 'Capreomycin', 'Kanamycin')


def read_parents(path_to_pheno_and_trees, drug, sample_ids):
    parents = []
    t = Tree()
    name_to_node = {}
    with open(path_to_pheno_and_trees + drug + '/parents.csv', 'r') as f:
        t.name = f.readline().strip()
        name_to_node[t.name] = t
        for line in f.readlines():
            s = line.strip().split('\t')
            parent = name_to_node[s[1]]
            child = parent.add_child(name=s[0], dist=float(s[2]))
            name_to_node[s[0]] = child
    t.prune(sample_ids)
    for node in t.iter_descendants("levelorder"):
        parents.append((node.name, node.up.name, str(node.dist)))
    return t.name, parents


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

    snp_to_cds = localize_all_variants(snp_pos_list_from_samples, cds_list)

    snp_pos_list_from_alignment = [int(l.strip()) for l in open(path_to_snps_list, 'r').readlines()]
    index_list, pos_list = filter_snp_list(snp_pos_list_from_alignment, all_snp_pos_set, snp_to_cds)
    print(str(len(index_list)) + ' snps after filtering')

    sample_sequences = SeqIO.parse(open(path_to_alignment, 'r'), 'fasta')
    sample_to_seq = filter_all_aln(sample_sequences, index_list)

    drug_to_number = {}
    if path_to_drug_codes is None:
        with open(out_path + '.drug_codes', 'w') as f:
            i = 0
            for drug in drug_names:
                drug_to_number[drug] = str(i)
                f.write(drug + '\t' + str(i))
                i += 1
    else:
        for line in open(path_to_drug_codes).readlines():
            s = line.strip().split('\t')
            drug_to_number[s[0]] = s[1]

    drug_to_pheno = {}
    root, parents = read_parents(path_to_pheno_and_trees, 'Streptomycin', sample_ids)
    for drug in drug_names:
        sample_to_pheno = read_pheno(path_to_pheno, path_to_pheno_and_trees, drug)
        drug_to_pheno[drug] = sample_to_pheno

    sample_to_mut = format_mut_lists(sample_to_seq, pos_list, ref_seq, snp_to_cds, parents)

    if exists(out_path) and not overwrite:
        return

    with open(out_path, 'w') as f:
        f.write("child\tparent\tlength\n")
        f.write(root + '\n')
        for node_id, parent_id, dist in parents:
            branch_line = [node_id, parent_id, dist]
            syn = sample_to_mut[node_id]
            branch_line.append(syn)
            nonsyn = []
            for drug in drug_names:
                sample_to_pheno = drug_to_pheno[drug]
                pheno = sample_to_pheno[node_id]
                parent_pheno = sample_to_pheno[parent_id]
                if parent_pheno != pheno or (print_RR and pheno == 'R'):
                    nonsyn.append(parent_pheno + drug_to_number[drug] + pheno)
            branch_line.append(';'.join(nonsyn))
            branch_line.append(';'.join(nonsyn))
            f.write('\t'.join(branch_line))
            f.write('\n')


if __name__ == '__main__':
    main()