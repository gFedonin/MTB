from Bio import SeqIO
from Bio.Seq import Seq

from src.core.annotations import read_annotations, localize_all_snps, CDSType
from src.core.constants import codon_table, codon_table_compl, complement, upstream_length

path_to_annotations = './data/AL123456_rev.gff'
path_to_ref = './data/h37rv.fasta'
path_to_snps = './res/Rifampicin_snps.csv'
out_path = './res/Rifampicin_snps_translated.csv'
path_to_snp_list = './data/snp_pos_with_DR.csv'

def read_h37rv():
    fasta_sequences = SeqIO.parse(open(path_to_ref), 'fasta')
    sequence = ''
    for fasta in fasta_sequences:
        sequence = str(fasta.seq)
    return sequence.upper()


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


def format_snp(snps, ref_seq, ref_seq_compl, snp_to_cds):

    res = []
    j = 0
    ref_seq_len = len(ref_seq)
    with open(out_path, 'w') as f:
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
                        old_aminoacid, new_aminoacid = get_aminoacids_sense(ref_seq, nucleotide_pos, pos, alt)

                        if old_aminoacid != new_aminoacid:
                            res.append('\t'.join(
                                (cds.type.name, cds.name, str(protein_pos), old_aminoacid, new_aminoacid)))
                    elif cds.type == CDSType.upstream:
                        res.append('\t'.join(
                            (cds.type.name, cds.name, str(pos - cds.end - 1), ref_seq[pos - 1], alt)))
                    else:
                        res.append('\t'.join(
                            (cds.type.name, cds.name, str(pos - cds.start + 1), ref_seq[pos - 1], alt)))
                else:  # if strand is '-'
                    if cds.type == CDSType.Gene:
                        protein_pos = (cds.end - pos) // 3 + 1
                        nucleotide_pos = (cds.end - pos) % 3

                        # NEIGHBOURS ANALYSIS
                        old_aminoacid, new_aminoacid = get_aminoacids_antisense(ref_seq, nucleotide_pos, pos, alt)

                        if old_aminoacid != new_aminoacid:
                            res.append('\t'.join(
                                (cds.type.name, cds.name, str(protein_pos), old_aminoacid, new_aminoacid)))
                    elif cds.type == CDSType.upstream:
                        res.append('\t'.join(
                            (cds.type.name, cds.name, str(cds.start - pos - 1), ref_seq[pos - 1], alt)))
                    else:
                        res.append('\t'.join(
                            (cds.type.name, cds.name, str(cds.end - (pos - 1)), ref_seq_compl[ref_seq_len - pos], complement[alt])))
            else:
                coord = ref_seq_len - pos
                res.append('\t'.join(('non_cds', ' ', str(coord + 1), ref_seq_compl[coord], complement[alt])))
            j += 1
        if len(res) > 0:
            f.write('\n'.join(res))
            f.write('\n')


def main():

    h37rv = read_h37rv()
    h37rv_compl = str(Seq(h37rv).reverse_complement())

    cds = read_annotations(path_to_annotations, upstream_length)

    genome_pos_list = [int(line.strip()) for line in open(path_to_snp_list, 'r').readlines()]

    all_snps = []
    all_snp_pos = []
    with open(path_to_snps, 'r') as f:
        for line in f.readlines():
            s = line.strip().split('\t')
            pos = genome_pos_list[int(s[0])]
            # pos = int(s[0]) - 1
            all_snp_pos.append(pos)
            all_snps.append((pos, s[1].upper()))

    snp_to_cds = localize_all_snps(all_snp_pos, cds)

    format_snp(all_snps, h37rv, h37rv_compl, snp_to_cds)


if __name__ == '__main__':
    main()
