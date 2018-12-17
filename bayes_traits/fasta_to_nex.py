import sys
from Bio import SeqIO
from Bio.Alphabet import generic_dna

if __name__ == '__main__':
    sequences = SeqIO.parse(open(sys.argv[1]), 'fasta', generic_dna)
    with open(sys.argv[2], "w") as f:
        SeqIO.write(sequences, f, "nexus")
