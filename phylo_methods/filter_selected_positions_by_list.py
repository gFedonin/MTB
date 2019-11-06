from src.core.constants import data_path

path_to_res_dir = '../../res/tree_was/mc10_mega_MP_RR_bam_filtered/'
p_value_threshold = 0.05
use_fdr = True
fdr_threshold = 0.1
use_pvalue_and_fdr = True
upper_or_lower = '.lower'

if use_fdr:
    if use_pvalue_and_fdr:
        filter_str = 'filtered_fdr' + str(fdr_threshold) + '_pvalue' + str(p_value_threshold)
    else:
        filter_str = 'filtered_fdr' + str(fdr_threshold)
else:
    filter_str = 'filtered_pvalue' + str(p_value_threshold)

path_to_selected_positions_list = path_to_res_dir + 'pairs' + upper_or_lower + '.' + filter_str + \
                                  '.converted.csv'
path_to_filter_positions_list = data_path + \
                                'snps/raw_with_DR_with_indel_with_pheno_and_snp_no_win_qual_mqm_std3_mqm30_no_highcov_bam_filtered/filtered_annotated_variants.csv'

out_path = path_to_res_dir + 'pairs' + upper_or_lower + '.' + filter_str + '.converted.filtered_by_mqm_bam_filtered.csv'
out_path_trash = path_to_res_dir + 'pairs' + upper_or_lower + '.' + filter_str + '.converted.trash.filtered_by_mqm_bam_filtered.csv'


if __name__ == '__main__':
    pos_to_filter = set()
    for line in open(path_to_filter_positions_list).readlines()[1:]:
        s = line.split('\t')
        pos_to_filter.add(s[1] + '\t' + s[2])
    with open(out_path, 'w') as f:
        with open(out_path_trash, 'w') as f1:
            fin = open(path_to_selected_positions_list)
            f.write(fin.readline())
            for line in fin.readlines():
                s = line.split('\t')
                if s[1] + '\t' + s[2] in pos_to_filter:
                    f.write(line)
                else:
                    f1.write(line)