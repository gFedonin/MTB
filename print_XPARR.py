from bisect import bisect_left, bisect_right
from ete3 import Tree
from Bio import SeqIO
from sklearn.externals.joblib import Parallel, delayed
import os

from annotations import read_annotations
from constants import codon_table, codon_table_compl

path_to_ids = './data/Full_subset_filtered_snp_pheno2.txt'
path_to_snps = './data/snps/raw_with_DR/'
path_to_tries = './data/tries/'
path_to_pheno = './data/pheno/'
out_path = './data/xparr/'
path_to_annotations = './data/AL123456_rev.gff'
path_to_ref = './data/h37rv.fasta'

upstream_length = 100


def read_h37rv():
    fasta_sequences = SeqIO.parse(open(path_to_ref), 'fasta')
    sequence = ''
    for fasta in fasta_sequences:
        sequence = str(fasta.seq)
    return sequence.upper()


def read_snps(sample_id):
    snps = []
    with open(path_to_snps + sample_id + '.snp', 'r') as f1:
        lines = f1.readlines()
        for line in lines:
            s = line[:-1].split('\t')
            snps.append((int(s[0]), s[1]))
    return sample_id, snps


def get_aminoacids_sense(ref_seq, nucleotide_pos, snps, i):
    snp_num = len(snps)
    pos, alt = snps[i]
    if nucleotide_pos == 0:
        aa0 = codon_table[ref_seq[pos - 1:pos + 2]]
        if i != snp_num - 1:
            pos1, alt1 = snps[i + 1]
            if pos1 == pos + 1:
                if i != snp_num - 2 and snps[i + 2][0] == pos + 2:
                    aa1 = codon_table[alt + alt1 + snps[i + 2][1]]
                    return aa0, aa1, i + 2
                else:
                    aa1 = codon_table[alt + alt1 + ref_seq[pos + 1]]
                    return aa0, aa1, i + 1
            elif pos1 == pos + 2:
                aa1 = codon_table[alt + ref_seq[pos] + alt1]
                return aa0, aa1, i + 1
        aa1 = codon_table[alt + ref_seq[pos] + ref_seq[pos + 1]]
        return aa0, aa1, i
    elif nucleotide_pos == 1:
        aa0 = codon_table[ref_seq[pos - 2:pos + 1]]
        if i != snp_num - 1 and snps[i + 1][0] == pos + 1:
            aa1 = codon_table[ref_seq[pos - 2] + alt + snps[i + 1][1]]
            return aa0, aa1, i + 1
        else:
            aa1 = codon_table[ref_seq[pos - 2] + alt + ref_seq[pos]]
            return aa0, aa1, i
    else:
        aa0 = codon_table[ref_seq[pos - 3:pos]]
        aa1 = codon_table[ref_seq[pos - 3:pos - 1] + alt]
        return aa0, aa1, i


def get_aminoacids_antisense(ref_seq, nucleotide_pos, snps, i):
    snp_num = len(snps)
    pos, alt = snps[i]
    if nucleotide_pos == 2:
        aa0 = codon_table_compl[ref_seq[pos - 1:pos + 2]]
        if i != snp_num - 1:
            pos1, alt1 = snps[i + 1]
            if pos1 == pos + 1:
                if i != snp_num - 2 and snps[i + 2][0] == pos + 2:
                    aa1 = codon_table_compl[alt + alt1 + snps[i + 2][1]]
                    return aa0, aa1, i + 2
                else:
                    aa1 = codon_table_compl[alt + alt1 + ref_seq[pos + 1]]
                    return aa0, aa1, i + 1
            elif pos1 == pos + 2:
                aa1 = codon_table_compl[alt + ref_seq[pos] + alt1]
                return aa0, aa1, i + 1
        aa1 = codon_table_compl[alt + ref_seq[pos:pos + 2]]
        return aa0, aa1, i
    elif nucleotide_pos == 1:
        aa0 = codon_table_compl[ref_seq[pos - 2:pos + 1]]
        if i != snp_num - 1 and snps[i + 1][0] == pos + 1:
            aa1 = codon_table_compl[ref_seq[pos - 2] + alt + snps[i + 1][1]]
            return aa0, aa1, i + 1
        aa1 = codon_table_compl[ref_seq[pos - 2] + alt + ref_seq[pos]]
        return aa0, aa1, i
    else:
        aa0 = codon_table_compl[ref_seq[pos - 3:pos]]
        aa1 = codon_table_compl[ref_seq[pos - 3:pos - 1] + alt]
        return aa0, aa1, i


def localize_all_snps(all_snps, cds_list):
    res = {}
    for cds in cds_list:
        i = bisect_left(all_snps, cds.start)
        j = bisect_right(all_snps, cds.end, lo=i)
        for k in range(i, j):
            res[all_snps[k]] = cds
    return res


def format_snp(sample_id, snps, ref_seq, snp_to_cds):
    nonsyn = []
    syn = []
    j = 0
    while j < len(snps):
        pos, alt = snps[j]
        try:
            cds = snp_to_cds[pos]
        except:
            cds = None

        if cds is not None:
            if cds.strand == 1:
                if cds.type == CDSType.Gene:
                    protein_pos = (pos - cds.start) // 3 + 1
                    nucleotide_pos = (pos - cds.start) % 3

                    # NEIGHBOURS ANALYSIS
                    old_aminoacid, new_aminoacid, j = get_aminoacids_sense(ref_seq, nucleotide_pos, snps, j)

                    if old_aminoacid != new_aminoacid:
                        nonsyn.append(old_aminoacid + cds.name + '_' + str(protein_pos) + new_aminoacid)
                    else:
                        syn.append(cds.name + '_' + str(protein_pos))
            else:  # if strand is '-'
                if cds.type == CDSType.Gene:
                    protein_pos = (cds.end - pos) // 3 + 1
                    nucleotide_pos = (cds.end - pos) % 3

                    # NEIGHBOURS ANALYSIS
                    old_aminoacid, new_aminoacid, j = get_aminoacids_antisense(ref_seq, nucleotide_pos, snps, j)

                    if old_aminoacid != new_aminoacid:
                        nonsyn.append(old_aminoacid + cds.name + '_' + str(protein_pos) + new_aminoacid)
                    else:
                        syn.append(cds.name + '_' + str(protein_pos))
        j += 1
    return sample_id, ';'.join(syn), ';'.join(nonsyn)


def main():
    sample_to_snps = {}
    all_snp_pos = set()
    h37rv = read_h37rv()

    sample_ids = [sample_id[:-1] for sample_id in open(path_to_ids, 'r').readlines()]

    cds = read_annotations(path_to_annotations, upstream_length)

    tasks = Parallel(n_jobs=-1)(
        delayed(read_snps)(sample_id) for sample_id in sample_ids)
    for sample_id, snps in tasks:
        sample_to_snps[sample_id] = snps
        for snp_pos, alt in snps:
            all_snp_pos.add(snp_pos)
    all_snps = list(all_snp_pos)
    all_snps.sort()
    print('done with snp')

    snp_to_cds = localize_all_snps(all_snps, cds)

    formatted_snps = Parallel(n_jobs=-1)(
        delayed(format_snp)(sample_id, snps, h37rv, snp_to_cds)
        for sample_id, snps in sample_to_snps.items()
    )
    print('done with snp format')

    sample_to_mut = {}
    for sample_id, syn, nonsyn in formatted_snps:
        sample_to_mut[sample_id] = (syn, nonsyn)

    drugs = []
    for (dirpath, dirnames, filenames) in os.walk(path_to_tries):
        for filename in filenames:
            drugs.append(filename.split('.')[0])
    for drug in drugs:
        t = Tree(path_to_tries + drug + '.nw')
        sample_to_pheno = {}
        with open(path_to_pheno + drug + '.pheno', 'r') as f:
            for line in f.readlines():
                s = line.split('\t')
                if s[1] == '1':
                    sample_to_pheno[s[0]] = 'R'
                else:
                    sample_to_pheno[s[0]] = 'S'
        with open(out_path, 'w') as f:
            f.write("child\tparent\tlength\n")
            root = t.get_tree_root()
            f.write(root.name)
            queue = [root]
            while (len(queue) != 0):
                node = queue.pop(0)
                if node.is_leaf():
                    branch_line = [node.name, node.up.name, str(node.dist)]
                    syn, nonsyn = sample_to_mut[node.name]
                    pheno = sample_to_pheno[node.name]
                    f.write('\t'.join(branch_line))
                    f.write('\n')
                else:
                    f.write(node.name + '\t' + node.up.name + '\t' + str(node.dist) + '\n')


if __name__ == '__main__':
    main()
