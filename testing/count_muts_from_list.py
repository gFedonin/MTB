from core.constants import data_path
from core.data_reading import read_all_variants

path_to_variants = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_no_win_qual_mqm_std3_mqm30_filter_samples_first/'
# path_to_list = path_to_variants + 'filtered_raw_variants_pos_converted.csv'
path_to_list = data_path + 'mc10_mega_MP_RR_filter/pairs.lower.filtered_fdr0.05.converted.filtered_by_mqm.csv'
path_to_data = data_path + 'snps/annotated_with_DR_with_indel_with_pheno_and_snp_no_win_qual_mqm_std3_mqm30_long_del_pg_filter_samples_first/'
path_to_ids = path_to_variants + 'samples_filtered.list'

out_path = data_path + 'mc10_mega_MP_RR_filter/counts.lower.filtered_fdr0.05.converted.filtered_by_mqm.csv'


if __name__ == '__main__':
    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    sample_to_variants = read_all_variants(path_to_data, sample_ids)
    variant_counts = {'\t'.join(l.split('\t')[0:3]): 0 for l in open(path_to_list).readlines()[1:]}
    for sample_id, variant_list in sample_to_variants.items():
        for var in variant_list:
            s = var.strip().split('\t')
            if s[0] == 'non_cds':
                mut = s[-1] + '\t' + s[0] + '\t' + s[2]
            else:
                mut = s[-1] + '\t' + s[1] + '\t' + s[2]
            c = variant_counts.get(mut)
            if c is not None:
                variant_counts[mut] = c + 1
    with open(out_path, 'w') as fout:
        for var, c in variant_counts.items():
            fout.write('%s\t%d\n' % (var, c))