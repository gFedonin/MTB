from os import makedirs

from os.path import exists

from src.core.annotations import read_annotations, CDSType
from src.core.constants import upstream_length

path_to_Walker_dict = './data/dictionaries_raw/Walker_dictionary.txt'
path_to_annotations = './data/AL123456_rev.gff'
out_path = './data/dictionaries_indels/Walker_dictionary.txt'


def main():
    cds_list = read_annotations(path_to_annotations, upstream_length)
    name_to_cds = {}
    for cds in cds_list:
        if cds.type == CDSType.upstream:
            continue
        name_to_cds[cds.name] = cds
        if cds.synonym != '':
            name_to_cds[cds.synonym] = cds
    with open(out_path, 'w') as f:
        for line in open(path_to_Walker_dict, 'r').readlines():
            l = line.strip().split('\t')
            s = l[0].split('_')
            name = s[0]
            cds = name_to_cds[name]
            if len(s) == 2:
                #snp
                pos = int(s[1][1:-1])
                if pos < 0:
                    f.write('\t'.join((CDSType.upstream.name, cds.name, s[1][1:-1], s[1][0], s[1][-1], 'snp', l[-1])))
                    f.write('\n')
                else:
                    f.write('\t'.join((cds.type.name, cds.name, s[1][1:-1], s[1][0], s[1][-1], 'snp', l[-1])))
                    f.write('\n')
            else:
                #indel
                if 'del' in s[2]:
                    f.write('\t'.join((cds.type.name, cds.name, s[1], s[2][3:], '-', 'del', l[-1])))
                    f.write('\n')
                else:
                    f.write('\t'.join((cds.type.name, cds.name, s[1], '-', s[2][3:], 'ins', l[-1])))
                    f.write('\n')


if __name__ == '__main__':
    main()
