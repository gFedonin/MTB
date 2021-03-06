from bisect import bisect_left

from os.path import exists, isfile

from Bio import SeqIO
from Bio.Seq import Seq


import os

from src.core.annotations import read_annotations, CDSType
from src.core.constants import codon_table, codon_table_compl, complement, upstream_length

path_to_old_dicts = '../../data/chukreev_dictionaries/'
out_path = '../../data/dictionaries_test/'
path_to_annotations = '../../data/AL123456_rev.gff'
path_to_ref = '../../data/h37rv.fasta'

nucl_list = ['A', 'T', 'G', 'C']

check_ref_seq = True


def get_aa_list():
    aa_set = set()
    for aa in codon_table.values():
        aa_set.add(aa)
    return list(aa_set)


def localize_snp(snp_pos, cds_list, cds_starts):
    i = bisect_left(cds_starts, snp_pos)
    if snp_pos == cds_list[i].start:
        return cds_list[i]
    if snp_pos < cds_list[i - 1].end:
        return cds_list[i - 1]
    else:
        return None


def get_aminoacids_sense(ref_seq, nucleotide_pos, pos, alt):
    if nucleotide_pos == 0:
        aa0 = codon_table[ref_seq[pos - 1:pos + 2]]
        aa1 = codon_table[alt + ref_seq[pos] + ref_seq[pos + 1]]
        return aa0, aa1
    elif nucleotide_pos == 1:
        aa0 = codon_table[ref_seq[pos - 2:pos + 1]]
        aa1 = codon_table[ref_seq[pos - 2] + alt + ref_seq[pos]]
        return aa0, aa1
    else:
        aa0 = codon_table[ref_seq[pos - 3:pos]]
        aa1 = codon_table[ref_seq[pos - 3:pos - 1] + alt]
        return aa0, aa1


def get_aminoacids_antisense(ref_seq, nucleotide_pos, pos, alt):
    if nucleotide_pos == 2:
        aa0 = codon_table_compl[ref_seq[pos - 1:pos + 2]]
        aa1 = codon_table_compl[alt + ref_seq[pos:pos + 2]]
        return aa0, aa1
    elif nucleotide_pos == 1:
        aa0 = codon_table_compl[ref_seq[pos - 2:pos + 1]]
        aa1 = codon_table_compl[ref_seq[pos - 2] + alt + ref_seq[pos]]
        return aa0, aa1
    else:
        aa0 = codon_table_compl[ref_seq[pos - 3:pos]]
        aa1 = codon_table_compl[ref_seq[pos - 3:pos - 1] + alt]
        return aa0, aa1


def translate_snp(cds, pos, ref_seq, alt):
    if cds.strand == 1:
        protein_pos = (pos - cds.start) // 3 + 1
        nucleotide_pos = (pos - cds.start) % 3
        old_aminoacid, new_aminoacid = get_aminoacids_sense(ref_seq, nucleotide_pos, pos, alt)
    else:
        protein_pos = (cds.end - pos) // 3 + 1
        nucleotide_pos = (cds.end - pos) % 3
        old_aminoacid, new_aminoacid = get_aminoacids_antisense(ref_seq, nucleotide_pos, pos, alt)
    return protein_pos, old_aminoacid, new_aminoacid


def convert_dict(name, cds_list, cds_starts, name_to_cds, ref_seq, ref_seq_compl, aa_list):
    ref_seq_len = len(ref_seq)
    list = []
    print('working with ' + name)
    with open(path_to_old_dicts + name, 'r') as f:
        for line in f.readlines():
            s = line.strip().split('\t')
            if s[0] == 'Gene':
                cds = name_to_cds[s[1]]
                if cds.type == CDSType.Gene:
                    s[1] = cds.name
                    if s[3] == '*' and s[4] == '*':
                        # Pankhurst ;)
                        protein_pos = int(s[2])
                        if cds.strand == 1:
                            pos = 3*(protein_pos - 1) + cds.start
                            old_aminoacid = codon_table[ref_seq[pos - 1:pos + 2]]
                        else:
                            pos = cds.end - 3 * (protein_pos - 1)
                            old_aminoacid = codon_table_compl[ref_seq[pos - 3:pos]]
                        s[3] = old_aminoacid
                        for aa in aa_list:
                            if aa != old_aminoacid:
                                s[4] = aa
                                list.append('\t'.join(s) + '\n')
                    else:
                        if check_ref_seq:
                            protein_pos = int(s[2])
                            if cds.strand == 1:
                                pos = 3 * (protein_pos - 1) + cds.start
                                old_aminoacid = codon_table[ref_seq[pos - 1:pos + 2]]
                            else:
                                pos = cds.end - 3 * (protein_pos - 1)
                                old_aminoacid = codon_table_compl[ref_seq[pos - 3:pos]]
                            if s[3] != old_aminoacid:
                                print('wrong ref_seq aa: ' + line)
                        list.append('\t'.join(s) + '\n')
                else:
                    s[0] = cds.type.name
                    s[1] = cds.name
                    if s[3] == '*' and s[4] == '*':
                        # Pankhurst ;)
                        if cds.strand == 1:
                            old = ref_seq[cds.start + int(s[2]) - 2]
                        else:
                            old = complement[ref_seq[cds.end - int(s[2])]]
                        s[3] = old
                        for na in nucl_list:
                            if na != old:
                                s[4] = na
                                list.append('\t'.join(s) + '\n')
                    else:
                        if check_ref_seq:
                            if cds.strand == 1:
                                old = ref_seq[cds.start + int(s[2]) - 2]
                            else:
                                old = complement[ref_seq[cds.end - int(s[2])]]
                            if s[3] != old:
                                print('wrong ref_seq na: ' + line)
                        list.append('\t'.join(s)+'\n')
            else:
                pos = int(s[2])
                cds = localize_snp(pos, cds_list, cds_starts)
                if cds is not None:
                    if cds.type == CDSType.Gene:
                        if s[3] == '*' and s[4] == '*':
                            # Pankhurst ;)
                            s[3] = ref_seq[pos - 1]
                            for na in nucl_list:
                                if na != s[3]:
                                    protein_pos, old_aminoacid, new_aminoacid = translate_snp(cds, pos, ref_seq, na)
                                    list.append(cds.type.name + '\t' + cds.name + '\t' + str(protein_pos) + '\t' +
                                                old_aminoacid + '\t' + new_aminoacid + '\t' + s[5] + '\n')
                        else:
                            protein_pos, old_aminoacid, new_aminoacid = translate_snp(cds, pos, ref_seq, s[4])
                            if old_aminoacid == new_aminoacid:
                                print('synonymous! ' + line)
                            list.append(cds.type.name + '\t' + cds.name + '\t' + str(protein_pos) + '\t' + old_aminoacid +
                                        '\t' + new_aminoacid + '\t' + s[5] + '\n')
                    elif cds.type == CDSType.upstream:
                        if cds.strand == 1:
                            upstream_pos = pos - cds.end - 1
                        else:
                            upstream_pos = cds.start - pos - 1
                        if s[3] == '*' and s[4] == '*':
                            # Pankhurst ;)
                            if cds.strand == 1:
                                old = ref_seq[pos - 1]
                            else:
                                old = ref_seq_compl[ref_seq_len - pos]
                            for na in nucl_list:
                                if na != old:
                                    list.append(cds.type.name + '\t' + cds.name + '\t' + str(upstream_pos) + '\t' +
                                                old + '\t' + na + '\t' + s[5] + '\n')
                        else:
                            if check_ref_seq:
                                if cds.strand == 1:
                                    old = ref_seq[pos - 1]
                                else:
                                    old = ref_seq_compl[ref_seq_len - pos]
                                    s[3] = complement[s[3]]
                                    s[4] = complement[s[4]]
                                if s[3] != old:
                                    print('wrong ref_seq na: ' + line)
                            list.append(cds.type.name + '\t' + cds.name + '\t' + str(upstream_pos) + '\t' + s[3] + '\t' +
                                        s[4] + '\t' + s[5] + '\n')
                    else:
                        if cds.strand == 1:
                            pos = pos - cds.start + 1
                        else:
                            pos = cds.end - pos + 1
                        if s[3] == '*' and s[4] == '*':
                            # Pankhurst ;)
                            if cds.strand == 1:
                                old = ref_seq[int(s[2]) - 1]
                            else:
                                old = ref_seq_compl[ref_seq_len - int(s[2])]
                            for na in nucl_list:
                                if na != old:
                                    list.append(cds.type.name + '\t' + cds.name + '\t' + str(pos) + '\t' +
                                                old + '\t' + na + '\t' + s[5] + '\n')
                        else:
                            if check_ref_seq:
                                if cds.strand == 1:
                                    old = ref_seq[int(s[2]) - 1]
                                else:
                                    old = ref_seq_compl[ref_seq_len - int(s[2])]
                                    s[3] = complement[s[3]]
                                    s[4] = complement[s[4]]
                                if s[3] != old:
                                    print('wrong ref_seq na: ' + line)
                            list.append(cds.type.name + '\t' + cds.name + '\t' + str(pos) + '\t' + s[3] + '\t' +
                                        s[4] + '\t' + s[5] + '\n')
                else:
                    pos = ref_seq_len - pos
                    if s[3] == '*' and s[4] == '*':
                        # Pankhurst ;)
                        for na in nucl_list:
                            if na != ref_seq_compl[pos]:
                                list.append('non_cds\t \t' + str(pos) + '\t' + ref_seq_compl[pos] + '\t' + na + '\t' +
                                            s[5] + '\n')
                    else:
                        s[3] = complement[s[3]]
                        if check_ref_seq:
                            if s[3] != ref_seq_compl[pos]:
                                print('wrong ref_seq na: ' + line)
                        list.append('non_cds\t \t' + str(pos) + '\t' + s[3] + '\t' + complement[s[4]] + '\t' + s[5] + '\n')
    with open(out_path + name, 'w') as f:
        for line in list:
            f.write(line)


def read_h37rv():
    fasta_sequences = SeqIO.parse(open(path_to_ref), 'fasta')
    sequence = ''
    for fasta in fasta_sequences:
        sequence = str(fasta.seq)
    return sequence.upper()


def main():
    if not exists(out_path):
        os.makedirs(out_path)
    h37rv = read_h37rv()
    h37rv_compl = str(Seq(h37rv).reverse_complement())
    aa_list = get_aa_list()
    cds_list = read_annotations(path_to_annotations, upstream_length)
    cds_starts = []
    name_to_cds = {}
    for cds in cds_list:
        cds_starts.append(cds.start)
        if cds.type != CDSType.upstream:
            name_to_cds[cds.name] = cds
            if cds.synonym != '':
                name_to_cds[cds.synonym] = cds
    filenames = [f for f in os.listdir(path_to_old_dicts) if isfile(path_to_old_dicts + f)]
    for filename in filenames:
        if '.txt' in filename:
            convert_dict(filename, cds_list, cds_starts, name_to_cds, h37rv, h37rv_compl, aa_list)


if __name__ == '__main__':
    main()
