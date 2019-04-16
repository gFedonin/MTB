from src.core.annotations import read_annotations, localize_all_variants, CDSType
from src.core.constants import data_path, upstream_length, ref_len

filter_path = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_no_win_qual_mqm_std3_mqm30_no_highcov_str10/'
path_to_raw_pos_list = filter_path + 'filtered_raw_variants_pos.csv'
path_to_indels = data_path + 'xparr/mc10_mega_MP/indel_list.txt'

out_path = filter_path + 'filtered_xparr_variant_pos.list'


def convert_raw_pos_to_xparr(pos, cds):
    if cds is None:
        return pos
    if cds.type == CDSType.Gene:
        if cds.strand == 1:
            nucleotide_pos = (pos - cds.start) % 3
            return pos - nucleotide_pos
        else:  # if strand is '-'
            nucleotide_pos = (cds.end - pos) % 3
            return pos - 2 + nucleotide_pos
    else:
        return pos


if __name__ == '__main__':
    raw_pos_list = [int(l.strip()) for l in open(path_to_raw_pos_list).readlines()]
    raw_pos_set = set(raw_pos_list)
    indel_list = [l.strip() for l in open(path_to_indels).readlines()]
    cds_list = read_annotations(upstream_length)
    pos_to_cds = localize_all_variants(raw_pos_list, cds_list)
    converted = set(convert_raw_pos_to_xparr(pos, pos_to_cds.get(pos)) for pos in raw_pos_list)
    i = 0
    for indel in indel_list:
        s = indel.split('_')
        pos = int(s[1])
        if pos in raw_pos_set:
            converted.add(ref_len + i + 1)
        i += 1
    converted_list = list(converted)
    converted_list.sort()
    with open(out_path, 'w') as f:
        f.write('\n'.join(map(str, converted_list)))
        f.write('\n')