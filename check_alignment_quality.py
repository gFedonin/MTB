import math
#from Bio import SeqIO
from sklearn.externals.joblib import Parallel, delayed

path_to_aln = './data/ancestors_sparse_merged.fasta'


def entropy(i, seqs):
    stat = 0
    counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    for seq in seqs:
        counts[seq[i]] += 1
    for c in counts.values():
        if c > 0:
            stat += -c / len(seqs) * math.log(c / len(seqs))
    return stat


def main():
    # seqs = SeqIO.parse(open(path_to_aln), 'fasta')
    stat_real = 0
    stat_mega = 0
    real = []
    mega = []
    # for seq in seqs:
    #     if 'Node' in seq.name:
    #         mega.append(seq.seq)
    #     else:
    #         real.append(seq.seq)
    with open(path_to_aln, 'r') as f:
        is_mega = True
        for line in f.readlines():
            if line[0] == '>':
                is_mega = 'Node' in line
            else:
                if is_mega:
                    mega.append(line)
                else:
                    real.append(line)

    aln_len = len(real[0])
    tasks = Parallel(n_jobs=-1)(delayed(entropy)(i, real) for i in range(aln_len))
    for task in tasks:
        stat_real += task
    stat_real /= aln_len
    print('real : %f' % stat_real)
    tasks = Parallel(n_jobs=-1)(delayed(entropy)(i, mega) for i in range(aln_len))
    for task in tasks:
        stat_mega += task
    stat_mega /= aln_len
    print('mega : %f' % stat_mega)


if __name__ == '__main__':
    main()