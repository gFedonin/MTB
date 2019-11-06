from sklearn.externals.joblib import Parallel, delayed

from src.core.annotations import read_annotations, localize_all_variants, CDSType
from src.core.constants import data_path, upstream_length, ref_len

filter_path = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_no_win_qual_mqm_std3_mqm30_no_highcov_bam_filtered/'
path_to_raw_pos_list = filter_path + 'filtered_raw_variants_pos.csv'
path_to_indels = data_path + 'xparr/mc10_mega_MP/indel_list.txt'
path_to_ids = filter_path + 'samples_filtered.list'#'test.list'

out_path = filter_path + 'filtered_xparr_variant_pos_bam_filtered.list'
out_path_indels = filter_path + 'filtered_indels.list'

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


def read_indels_from_sample(sample_id):
    indels = []
    in_del = False
    del_pos = -1
    del_len = 0
    for line in open(filter_path + sample_id + '.variants').readlines():
        s = line.strip().split('\t')
        if s[-1] == 'ins':
            if in_del:
                in_del = False
            indels.append('ins_' + s[0] + '_' + s[1])
        elif s[-1] == 'del':
            if in_del:
               del_len += 1
            else:
                in_del = True
                del_len = 1
                del_pos = s[0]
        else:
            if in_del:
                indels.append('del_' + del_pos + '_' + str(del_len))
                in_del = False
    if in_del:
        indels.append('del_' + del_pos + '_' + str(del_len))
    return indels


def read_indels_from_all_samples():
    tasks = Parallel(n_jobs=-1)(delayed(read_indels_from_sample)(l.strip()) for l in open(path_to_ids).readlines())
    indels = set()
    for task in tasks:
        indels.update(task)
    return indels


if __name__ == '__main__':
    raw_pos_list = [int(l.strip()) for l in open(path_to_raw_pos_list).readlines()]
    raw_pos_set = set(raw_pos_list)
    indel_list = [l.strip() for l in open(path_to_indels).readlines()]
    cds_list = read_annotations(upstream_length)
    pos_to_cds = localize_all_variants(raw_pos_list, cds_list)
    converted = set(convert_raw_pos_to_xparr(pos, pos_to_cds.get(pos)) for pos in raw_pos_list)
    converted_list = list(converted)
    converted_list.sort()
    filtered_indels = read_indels_from_all_samples()
    with open(out_path_indels, 'w') as fout:
        i = 0
        for indel in indel_list:
            if indel in filtered_indels:
                converted_list.append(ref_len + i + 1)
                fout.write(indel + '\n')
            i += 1
    with open(out_path, 'w') as f:
        f.write('\n'.join(map(str, converted_list)))
        f.write('\n')