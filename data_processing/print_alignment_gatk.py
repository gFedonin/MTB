from os import listdir

from Bio import SeqIO
from sklearn.externals.joblib import Parallel, delayed

from src.core.constants import data_path
from src.core.data_reading import read_h37rv

path_to_snps = data_path + 'snps/gatk_before_cortex/raw_variants/'
out_path_aln = data_path + 'gatk_before_cortex.fasta'
out_path_snps = data_path + 'gatk_before_cortex_snps.txt'
path_to_mummer = data_path + 'mummer1.aligns'

thread_num = 32


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
                    aln.append(
                        (h37rv_from, h37rv_to, ''.join(aln1), ''.join(aln2)))
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


def read_snps(fname):
    sample_id = fname[0: fname.rfind('.variants')]
    variants = {}
    with open(path_to_snps + fname) as f1:
        for line in f1.readlines():
            if line[0] == '#':
                continue
            s = line.strip().split('\t')
            if int(s[-1]) == 1:
                variants[int(s[0])] = s[1]
    return sample_id, variants


def main():
    sample_to_snps = {}
    all_snp_pos = set()

    fnames = [fname for fname in listdir(path_to_snps)]

    h37rv = read_h37rv()

    aln = parse_mummer()
    aln = preprocess_aln(aln)
    print('mummer parsed')

    tasks = Parallel(n_jobs=thread_num, batch_size=len(
        fnames)//thread_num + 1)(delayed(read_snps)(fname) for fname in fnames)
    for sample_id, snps in tasks:
        sample_to_snps[sample_id] = snps
        for pos, alt in snps.items():
            all_snp_pos.add(pos)
    all_snps = list(all_snp_pos)
    all_snps.sort()

    with open(out_path_snps, 'w') as f:
        for pos in all_snps:
            f.write(str(pos))
            f.write('\n')

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
