from os.path import exists

from os import makedirs

path_to_aln = './data/snp_aln_with_DR_with_pheno_and_snp_mc5.fasta'
out_path = './data/split_aln_mc5/'
path_to_subset = ''

sequential = True
split_num = 144


def main():
    if not exists(out_path):
        makedirs(out_path)
    seq_splits = []
    with open(path_to_aln, 'r') as f:
        name = ''
        for line in f.readlines():
            if line[0] == '>':
                name = line
            else:
                fragments = []
                seq_splits.append((name, fragments))
                if sequential:
                    block_size = len(line)//split_num + 1
                    for i in range(0, len(line) - 1, block_size):
                        if i + block_size < len(line):
                            fragments.append(line[i:i + block_size])
                        else:
                            fragments.append(line[i:len(line) - 1])
                else:
                    lists = []
                    for i in range(split_num):
                        lists.append([])
                    for i in range(len(line) - 1):
                        lists[i % split_num].append(line[i])
                    for i in range(split_num):
                        fragments.append(''.join(lists[i]))
    for i in range(split_num):
        with open(out_path + 'split_' + str(i) + '.fasta', 'w') as f:
            for name, fragments in seq_splits:
                f.write(name)
                f.write(fragments[i])
                f.write('\n')


if __name__ == '__main__':
    main()
