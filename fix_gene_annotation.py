from Bio import SeqIO

from annotations import read_annotations, fix_gene_annotations

path_to_annotations = './data/AL123456_rev.gff'
path_to_ref = './data/h37rv.fasta'

upstream_length = 100

def read_h37rv():
    fasta_sequences = SeqIO.parse(open(path_to_ref), 'fasta')
    sequence = ''
    for fasta in fasta_sequences:
        sequence = str(fasta.seq)
    return sequence.upper()


h37rv = read_h37rv()
cds_list = read_annotations(path_to_annotations, upstream_length)
fix_gene_annotations(cds_list, h37rv)
