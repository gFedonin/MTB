from Bio import SeqIO
import numpy as np
from bisect import bisect_left, bisect_right

from sklearn.externals.joblib import Parallel, delayed

path_to_ids = './data/Full_subset_filtered_snp_pheno2.txt'
path_to_snps = '/data/snps/raw_no_DR/'
out_path_aln = './data/snp_aln_no_DR_filtered2.fasta'
path_to_mummer = './data/mummer1.aligns'
path_to_ref = './data/h37rv.fasta'


def read_h37rv():
    fasta_sequences = SeqIO.parse(open(path_to_ref), 'fasta')
    sequence = ''
    for fasta in fasta_sequences:
        sequence = str(fasta.seq)
    return sequence.upper()


def parse_mummer():
    aln = []
    aln1 = []
    aln2 = []
    h37rv_from = 0
    h37rv_to = 0
    with open(path_to_mummer, 'r') as f:
        for line in f:
            if line.startswith('-- BEGIN'):
                if len(aln1) > 0:
                    aln.append((h37rv_from, h37rv_to, ''.join(aln1), ''.join(aln2)))
                    aln1.clear()
                    aln2.clear()
                s = line.split(' ')
                h37rv_from = int(s[5])
                h37rv_to = int(s[7])
            else:
                if line[0] not in ('\n', ' ', '-', '=', '/'):
                    s = line.split()
                    if len(aln1) == len(aln2):
                        aln1.append(s[1])
                    else:
                        aln2.append(s[1])
    aln.append((h37rv_from, h37rv_to, ''.join(aln1), ''.join(aln2)))
    aln.sort(key=lambda tup: tup[0])
    curr_segment = aln[0]
    aln_filtered = []
    for segment in aln:
        if segment[0] <= curr_segment[1]:
            # intersection
            if curr_segment[1] - curr_segment[0] < segment[1] - segment[0]:
                curr_segment = segment
        else:
            aln_filtered.append(curr_segment)
            curr_segment = segment
    aln_filtered.append(curr_segment)
    # with open('./test.txt', 'w') as f:
    #     for segment in aln_filtered:
    #         f.write(segment[2] + '\n')
    #         f.write(segment[3] + '\n')
    return aln_filtered


def preprocess_aln(aln):
    aln_new = []
    for segment in aln:
        cannetti_seq = []
        aln1 = segment[2]
        aln2 = segment[3]
        if len(aln1) != len(aln2):
            print('wrong lens: ' + str(segment[0]) + ' ' + str(segment[1]))
        start = 0
        is_gap = False
        for i in range(len(aln1)):
            if aln1[i] == '.':
                if not is_gap:
                    cannetti_seq.append(aln2[start:i])
                    is_gap = True
            else:
                if is_gap:
                    start = i
                    is_gap = False
        if start > 0:
            cannetti_seq.append(aln2[start:len(aln2)])
            aln_new.append((segment[0], segment[1], ''.join(cannetti_seq)))
        else:
            aln_new.append((segment[0], segment[1], aln2))
    # for segment in aln_new:
    #     if segment[1] - segment[0] != len(segment[2]) - 1:
    #         print('post proc err: ' + str(segment[0]) + ' ' + str(segment[1]))
    return aln_new


def h37rv_to_canetti(aln, pos):
    for segment in aln:
        if segment[0] <= pos <= segment[1]:
            if segment[2][pos - segment[0]] == '.':
                return '-'
            return segment[2][pos - segment[0]]
    return ''


def read_snps(sample_id):
    snps = {}
    with open(path_to_snps + sample_id + '.snp', 'r') as f1:
        lines = f1.readlines()
        for line in lines:
            s = line[:-1].split('\t')
            snps[int(s[0])] = s[1]
    return sample_id, snps


def main():
    sample_to_snps = {}
    all_snp_pos = set()

    sample_ids = [sample_id[:-1] for sample_id in open(path_to_ids, 'r').readlines()]

    h37rv = read_h37rv()

    aln = parse_mummer()
    aln = preprocess_aln(aln)
    print('mummer parsed')

    tasks = Parallel(n_jobs=-1)(delayed(read_snps)(sample_id) for sample_id in sample_ids)
    for sample_id, snps in tasks:
        sample_to_snps[sample_id] = snps
        for snp_pos in snps.keys():
            all_snp_pos.add(snp_pos)
    all_snps = list(all_snp_pos)
    all_snps.sort()

    with open(out_path_aln, 'w') as f:
        f.write('>H37Rv\n')
        for snp_pos in all_snps:
            f.write(h37rv[snp_pos - 1])
        f.write('\n')
        f.write('>canetti\n')
        for snp_pos in all_snps:
            n = h37rv_to_canetti(aln, snp_pos)
            if n != '':
                f.write(n.upper())
            else:
                f.write(h37rv[snp_pos - 1])
        f.write('\n')
        for sample_id, snps in sample_to_snps.items():
            f.write('>' + sample_id + '\n')
            for snp_pos in all_snps:
                snp_letter = snps.get(snp_pos)
                if snp_letter is not None:
                    f.write(snp_letter)
                else:
                    f.write(h37rv[snp_pos - 1])
            f.write('\n')


if __name__ == '__main__':
    main()
