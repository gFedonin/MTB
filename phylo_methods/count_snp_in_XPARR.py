from os import listdir

from src.core.constants import data_path, ref_len

path_to_muts = '../../res/tree_was/filtered_no_phylogenetc_markers_no_pgrs/'
path_to_XPARR = data_path + 'xparr/mc10_mega_MP_vars_vs_vars.xparr'


def read_selected_pos(path):
    snps = []
    indels = []
    for line in open(path).readlines()[1:]:
        s = line.strip().split('\t')[0]
        pos = int(s)
        if pos > ref_len:
            indels.append(s)
        else:
            snps.append(s)
    return snps, indels


def count_mut_in_file(path_to_XPARR, snp_counts, indel_counts):
    for line in open(path_to_XPARR).readlines()[2:]:
        for mut in line.strip().split('\t')[-1].split(';'):
            try:
                pos = int(mut)
                if mut in indel_counts:
                    indel_counts[mut] += 1
            except:
                if mut[1:-1] in snp_counts:
                    snp_counts[mut[1:-1]] += 1


if __name__ == '__main__':
    for f in listdir(path_to_muts):
        print(f[:f.find('.')].upper())
        snps, indels = read_selected_pos(path_to_muts + f)
        snp_counts = {snp: 0 for snp in snps}
        indel_counts = {indel: 0 for indel in indels}
        if path_to_XPARR[-1] == '/':
            for path in path_to_XPARR:
                count_mut_in_file(path, snp_counts, indel_counts)
        else:
            count_mut_in_file(path_to_XPARR, snp_counts, indel_counts)
        print('snps:')
        for snp, count in snp_counts.items():
            print(snp + ' ' + str(count))
        print('indels:')
        for indel, count in indel_counts.items():
            print(indel + ' ' + str(count))
        print()

