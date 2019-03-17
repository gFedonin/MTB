from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.pairwise2 import format_alignment


for a in pairwise2.align.globalds("KEVLA", "EVL", matlist.blosum62, -7, -1):
    print(format_alignment(*a))