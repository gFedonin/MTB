from src.core.constants import ref_len, complement

in_path_triplets = '../../data/phylogenetic markers.csv'
out_path = '../../data/phylogenetic_markers_converted.csv'
in_path_single = '../../data/phylogenetic markers Beijing CAO.csv'


def convert_triplets():
    snps = []
    with open(in_path_triplets) as fin:
        for line in fin.readlines()[1:]:
            s = line.strip().split('\t')
            pos = int(s[1])
            triplet = s[2].split('/')[1]
            for i in range(3):
                if triplet[i].isupper():
                    pos += i
                    alt = triplet[i]
                    break
            snps.append(str(ref_len - pos + 1) + '\t' + complement[alt])
    with open(out_path, 'w') as fout:
        fout.write('\n'.join(snps))
        fout.write('\n')


def converts_single():
    snps = []
    with open(in_path_single) as fin:
        for line in fin.readlines()[1:]:
            s = line.strip().split('\t')
            pos = int(s[1])
            alt = s[2][2]
            snps.append(str(ref_len - pos + 1) + '\t' + complement[alt])
    with open(out_path, 'a') as fout:
        fout.write('\n'.join(snps))
        fout.write('\n')


if __name__ == '__main__':
    convert_triplets()
    converts_single()
