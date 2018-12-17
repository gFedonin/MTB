path_to_splits = '../../data/temp/'
out_path = '../../data/ancestors_sparse_mc10_mega_merged.fasta'

sequential = False
split_num = 144


def main():

    seq_splits = []
    with open(path_to_splits + 'split_' + str(0) + '.fasta', 'r') as f:
        name = ''
        for line in f.readlines():
            if line[0] == '>':
                name = line
            else:
                fragments = []
                seq_splits.append((name, fragments))
                fragments.append(line.strip())
    for i in range(1, split_num):
        with open(path_to_splits + 'split_' + str(i) + '.fasta', 'r') as f:
            j = 0
            for line in f.readlines():
                if line[0] != '>':
                    seq_splits[j][1].append(line.strip())
                    j += 1

    with open(out_path, 'w') as f:
        if sequential:
            for name, fragments in seq_splits:
                f.write(name)
                f.write(''.join(fragments))
                f.write('\n')
        else:
            for name, fragments in seq_splits:
                f.write(name)
                list = []
                aln_len = len(fragments[0])
                for i in range(aln_len):
                    for s in fragments:
                        if i < len(s):
                            list.append(s[i])
                list.append('\n')
                f.write(''.join(list))


if __name__ == '__main__':
    main()
