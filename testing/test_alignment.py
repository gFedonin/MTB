from Bio import SeqIO
from sklearn.externals.joblib import Parallel, delayed

path_to_alignment = './data/snp_aln_no_DR_filtered2.fasta'
out_path = './data/snp_aln_no_DR_filtered2.dist'

def check_seq(h37rv, seq):
    c = 0
    for i in range(len(h37rv)):
        if h37rv[i] != seq.seq[i]:
            c += 1
    return seq.name + '\t' + str(c) + '\n'


def main():
    fasta_sequences = SeqIO.parse(open(path_to_alignment), 'fasta')
    h37rv = ''
    seq_list = []
    for seq in fasta_sequences:
        if seq.name == 'H37Rv':
            h37rv = str(seq.seq)
        else:
            seq_list.append(seq)
    tasks = Parallel(n_jobs=-1)(delayed(check_seq)(h37rv, seq) for seq in seq_list)
    with open(out_path, 'w') as f:
        for line in tasks:
            f.write(line)


if __name__ == '__main__':
    main()