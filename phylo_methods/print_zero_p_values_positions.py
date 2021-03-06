from src.core.annotations import read_annotations, localize_all_variants
from src.core.constants import upstream_length, data_path, ref_len

phylogenetic_markers_list_path = '../../data/phylogenetic_markers_converted.csv'
path_to_indel_list = data_path + 'xparr/mc10_mega_MP/indel_list.txt'
drug_list = ['PZA', 'RIF', 'AMI', 'CAP', 'CIP', 'EMB', 'ETH', 'INH', 'KAN', 'MOX', 'OFX', 'PRO', 'STR']#

filter_out_phylogenetic_markers = True
filter_genes = 'PGRS'


def print_zero_p_values():
    phylogenetic_markers_pos_set = None
    if filter_out_phylogenetic_markers:
        phylogenetic_markers_pos_set = set()
        with open(phylogenetic_markers_list_path) as fin:
            for line in fin.readlines():
                phylogenetic_markers_pos_set.add(line.split('\t')[0])
    for drug in drug_list:
        print(drug)
        path_to_coords = '../../res/tree_was/' + drug + '/' + drug.lower() + '.site_pairs'
        path_to_p_values = '../../res/tree_was/' + drug + '/' + drug.lower() + '.upper.pvalue'
        out_path = '../../res/tree_was/filtered_no_phylogenetc_markers/' + drug.lower() + '.upper.pvalue.pairs'
        with open(out_path, 'w') as f:
            p_val_iter = open(path_to_p_values)
            coords = open(path_to_coords)
            f.write(coords.readline())
            for line in coords.readlines():
                p_val = p_val_iter.readline().strip()
                if float(p_val) == 0:
                    if phylogenetic_markers_pos_set is not None:
                        s = line.split('\t')[0]
                        if s not in phylogenetic_markers_pos_set:
                            f.write(line)
                        else:
                            print(line)
                    else:
                        f.write(line)


def print_no_filtered_genes():
    indel_list = [line.strip() for line in open(path_to_indel_list).readlines()]
    cds_list = read_annotations(upstream_length)
    for drug in drug_list:
        snp_list = []
        indel_pos_list = []
        pos_to_indel = {}
        with open('../../res/tree_was/filtered_no_phylogenetc_markers/' + drug.lower() + '.upper.pvalue.pairs') as fin:
            for line in fin.readlines()[1:]:
                pos = int(line.split('\t')[0])
                if pos <= ref_len:
                    snp_list.append(pos)
                else:
                    s = indel_list[pos - ref_len - 1].split('_')
                    indel_pos_list.append(int(s[1]))
                    pos_to_indel[int(s[1])] = str(pos)

        snp_pos_to_cds = localize_all_variants(snp_list, cds_list)
        indel_pos_to_cds = localize_all_variants(indel_pos_list, cds_list)
        c = 0
        with open('../../res/tree_was/filtered_no_phylogenetc_markers_no_pgrs/' + drug.lower() + '.upper.pvalue.pairs', 'w') as fout:
            for pos, cds in snp_pos_to_cds.items():
                if filter_genes not in cds.name:
                    fout.write(str(pos))
                    fout.write('\n')
                else:
                    c += 1
            for pos, cds in indel_pos_to_cds.items():
                if filter_genes not in cds.name:
                    fout.write(pos_to_indel[pos])
                    fout.write('\n')
                else:
                    c += 1
        print('found ' + str(c) + ' ' + filter_genes + ' in ' + drug)


def merge_pos_lists():
    snp_pos_set = set()
    indel_pos_set = set()
    indel_list = [line.strip() for line in open(path_to_indel_list).readlines()]
    for drug in drug_list:
        with open('../../res/tree_was/filtered_no_phylogenetc_markers_no_pgrs/' + drug.lower() + '.upper.pvalue.pairs') as fin:
            for line in fin.readlines()[1:]:
                pos = int(line.split('\t')[0])
                if pos <= ref_len:
                    snp_pos_set.add(pos)
                else:
                    s = indel_list[pos - ref_len - 1].split('_')
                    indel_pos_set.add(int(s[1]))
    filtered_snp_list = list(snp_pos_set)
    filtered_snp_list.sort()
    filtered_snp_list = list(map(int, filtered_snp_list))
    with open('../../res/tree_was/merged_snp.pos', 'w') as fout:
        fout.write('\n'.join(filtered_snp_list))
        fout.write('\n')
    filtered_indel_list = list(indel_pos_set)
    filtered_indel_list.sort()
    filtered_indel_list = list(map(int, filtered_indel_list))
    with open('../../res/tree_was/merged_indel.pos', 'w') as fout:
        fout.write('\n'.join(filtered_indel_list))
        fout.write('\n')


if __name__ == '__main__':
    # print_no_filtered_genes()
    merge_pos_lists()