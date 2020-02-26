from Bio import SeqIO
from sklearn.externals.joblib import Parallel, delayed

from core.constants import data_path

path_to_alignment = data_path + 'combined_mq40_keep_complex_std_names_filtered_no_DR.fasta'
out_path = '../../res/combined_mq40_keep_complex_std_names_filtered_no_DR.dist'


def check_seq(h37rv, seq):
    c = 0
    for i in range(len(h37rv)):
        if h37rv[i] != seq.seq[i]:
            c += 1
    return seq.name + '\t' + str(c) + '\n'


def check_pairs():
    fasta_sequences = [r for r in SeqIO.parse(open(path_to_alignment), 'fasta')]
    for i in range(len(fasta_sequences)):
        seq_i = fasta_sequences[i].seq
        for j in range(i + 1, len(fasta_sequences)):
            seq_j = fasta_sequences[j].seq
            if seq_i == seq_j:
                print(fasta_sequences[i].id + ' ' + fasta_sequences[j].id)


def compute_dist_to_ref():
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
    # compute_dist_to_ref()
    check_pairs()