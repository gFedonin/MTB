from Bio import SeqIO

from src.core.annotations import read_annotations, CDSType
from src.core.constants import codon_table, complement, codon_table_compl

path_to_annotations = './data/AL123456_rev.gff'
path_to_ref = './data/h37rv.fasta'

upstream_length = 100

def read_h37rv():
    fasta_sequences = SeqIO.parse(open(path_to_ref), 'fasta')
    sequence = ''
    for fasta in fasta_sequences:
        sequence = str(fasta.seq)
    return sequence.upper()


def fix_gene_annotations(cds_list, ref_seq):
    for cds in cds_list:
        if cds.type != CDSType.Gene:
            continue
        cds_len = cds.end - cds.start + 1
        if cds_len % 3 == 0:
            continue
        gene_seq = ref_seq[cds.start - 1: cds.end]
        if cds.strand == 1:
            i = cds_len // 3
            last_codon = gene_seq[3*(i - 1): 3*i]
            next_codon = ref_seq[cds.start - 1 + 3*i: cds.start - 1 + 3*i + 3]
            if codon_table[last_codon] == '*':
                cds.end -= cds_len % 3
                print('fixed cds ' + cds.name)
            elif codon_table[next_codon] == '*':
                cds.end += 3 - cds_len % 3
                print('fixed cds ' + cds.name)
            else:
                print('can\'t fix the cds ' + cds.name)
        else:
            gene_seq = gene_seq.translate(complement)[::-1]
            i = cds_len // 3
            last_codon = gene_seq[3*(i - 1): 3*i]
            next_codon = ref_seq[cds.start - 1 - 3 + cds_len % 3: cds.start - 1 - 3 + cds_len % 3 + 3]
            if codon_table[last_codon] == '*':
                cds.start += cds_len % 3
                print('fixed cds ' + cds.name)
            elif codon_table_compl[next_codon] == '*':
                cds.start -= 3 - cds_len % 3
                print('fixed cds ' + cds.name)
            else:
                print('can\'t fix the cds ' + cds.name)


h37rv = read_h37rv()
cds_list = read_annotations(path_to_annotations, upstream_length)
fix_gene_annotations(cds_list, h37rv)
