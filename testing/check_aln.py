import math

path_to_aln = './data/ancestors_merged.fasta'


def main():
    print('sequential')
    stat_real = 0
    stat_mega = 0
    real_counts = []
    mega_counts = []
    seq_num_mega = 0
    seq_num_real = 0
    with open(path_to_aln, 'r') as f:
        is_mega = True
        for line in f.readlines():
            if line[0] == '>':
                is_mega = 'Node' in line
            else:
                if len(real_counts) == 0:
                    aln_len = len(line) - 1
                    for i in range(aln_len):
                        real_counts.append([0, 0, 0, 0])
                        mega_counts.append([0, 0, 0, 0])
                if is_mega:
                    for i in range(aln_len):
                        n = line[i]
                        if n == 'A':
                            mega_counts[i][0] += 1
                        elif n == 'T':
                            mega_counts[i][1] += 1
                        elif n == 'G':
                            mega_counts[i][2] += 1
                        elif n == 'C':
                            mega_counts[i][3] += 1
                    seq_num_mega += 1
                else:
                    for i in range(aln_len):
                        n = line[i]
                        if n == 'A':
                            real_counts[i][0] += 1
                        elif n == 'T':
                            real_counts[i][1] += 1
                        elif n == 'G':
                            real_counts[i][2] += 1
                        elif n == 'C':
                            real_counts[i][3] += 1
                    seq_num_real += 1

    for i in range(aln_len):
        for c in real_counts[i]:
            if c > 0:
                stat_real += -c / seq_num_real * math.log(c / seq_num_real)
    stat_real /= aln_len
    print('real : %f' % stat_real)
    for i in range(aln_len):
        for c in mega_counts[i]:
            if c > 0:
                stat_mega += -c / seq_num_mega * math.log(c / seq_num_mega)
    stat_mega /= aln_len
    print('mega : %f' % stat_mega)


if __name__ == '__main__':
    main()