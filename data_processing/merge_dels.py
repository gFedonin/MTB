from os import listdir, makedirs
from os.path import exists

from sklearn.externals.joblib import Parallel, delayed

from core.constants import data_path
from core.data_reading import read_all_variants

# path_to_vars = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_no_win_qual_mqm_std3_mqm30_no_highcov_bam_filtered/'
# out_path = data_path + 'snps/raw_no_win_qual_mqm_std3_mqm30_no_highcov_bam_filtered_ld/'
# path_to_vars = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_no_win_qual_mqm_std3_mqm30_no_highcov/'
# out_path = data_path + 'snps/raw_no_win_qual_mqm_std3_mqm30_no_highcov_ld/'
# path_to_vars = data_path + 'snps/raw_freebayes_mqm_std3_mqm30_no_highcov_no_repeats_no_cov/'
# out_path = data_path + 'snps/raw_no_win_qual_mqm_std3_mqm30_no_highcov_no_repeats_no_cov_ld/'
# path_to_vars = data_path + 'snps/gatk_before_cortex/raw_variants_fixed_no_rep_gatk/'
# out_path = data_path + 'snps/gatk_before_cortex/raw_variants_fixed_no_rep_gatk_ld/'
# path_to_vars = data_path + 'snps/freebayes_after_cortex/raw_no_win_qual_mqm_std3_mqm30_no_highcov/'
# out_path = data_path + 'snps/freebayes_after_cortex/raw_no_win_qual_mqm_std3_mqm30_no_highcov_ld/'
# path_to_vars = data_path + 'snps/gatk_before_cortex/raw_variants_no_gvcf_mc10/'
# out_path = data_path + 'snps/gatk_before_cortex/raw_variants_no_gvcf_mc10_ld/'
# path_to_vars = data_path + 'snps/gatk_minimap/raw_variants/'
# out_path = data_path + 'snps/gatk_minimap/raw_variants_ld/'
# path_to_vars = data_path + 'snps/cortex/raw_variants/'
# out_path = data_path + 'snps/cortex/raw_variants_ld/'
path_to_vars = data_path + 'snps/skesa_bwa_mapq_raw/'
out_path = data_path + 'snps/skesa_bwa_mapq_raw_ld/'


def print_with_merged_dels(sample_id, variants):
    with open(out_path + sample_id + '.variants', 'w') as f:
        is_del = False
        start = -1
        l = 0
        for var in variants:
            s = var.split('\t')
            if s[-1] == 'del':
                pos = int(s[0])
                if is_del:
                    if pos == start + l:
                        l += 1
                    else:
                        f.write(str(start) + '\t' + str(l) + '\tdel\n')
                        start = pos
                        l = 1
                else:
                    is_del = True
                    start = pos
                    l = 1
            else:
                if is_del:
                    f.write(str(start) + '\t' + str(l) + '\tdel\n')
                    start = -1
                    l = 0
                    is_del = False
                f.write(var + '\n')
    return 1


def merge_all_dels():
    if not exists(out_path):
        makedirs(out_path)
    sample_ids = [fname[:-len('.variants')] for fname in listdir(path_to_vars) if fname.endswith('.variants')]
    sample_to_variants = read_all_variants(path_to_vars, sample_ids)
    tasks = Parallel(n_jobs=-1)(delayed(print_with_merged_dels)(sample, variants)
                                for sample, variants in sample_to_variants.items())
    c = 0
    for task in tasks:
        c += task
    print('%d samples processed' % c)


if __name__ == "__main__":
    merge_all_dels()
