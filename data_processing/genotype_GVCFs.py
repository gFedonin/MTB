from os import system, makedirs
from os.path import exists

from core.constants import data_path, ref_len

path_to_GATK = '/export/home/fedonin/gatk-4.1.4.0/gatk'
path_to_h37rv = data_path + 'h37rv.fasta'
# gatk_out_path = data_path + 'snps/combined_nocoll/'
# gatk_out_path = data_path + 'snps/gatk_before_cortex/'
# gatk_out_path = data_path + 'snps/gatk_combined/'
# gatk_out_path = data_path + 'coll18/snps/gatk/'
# gatk_out_path = data_path + 'walker18/snps/gatk/'
gatk_out_path = data_path + 'missing/snps/gatk/'
# gatk_out_path = data_path + 'snps/gatk_bowtie_rev_l15n1/'
# gvcf_path = gatk_out_path + 'gvcf_force_active_no_opt/'
gvcf_path = gatk_out_path + 'gvcf_bwa/'
# log_path_db = gatk_out_path + 'genomics_DB_import_force_active_no_opt.log'
# log_path_db = gatk_out_path + 'genomics_DB_import_bwa_combined3.log'
# log_path_db = gatk_out_path + 'genomics_DB_import_bwa_trimmed34.log'
log_path_db = gatk_out_path + 'genomics_DB_import_bwa_missing.log'
# log_path_genotype = gatk_out_path + 'genotype_force_active_no_opt.log'
# log_path_genotype = gatk_out_path + 'genotype_bwa_combined3.log'
# log_path_genotype = gatk_out_path + 'genotype_bwa_trimmed34.log'
log_path_genotype = gatk_out_path + 'genotype_bwa_missing.log'
temp_path = gatk_out_path + 'tmp/'
# path_to_ids = data_path + 'all_with_pheno.txt'
# path_to_ids = data_path + 'coll18/coll18_supp.samples'
# path_to_ids = data_path + 'walker18/trimming34.list'
path_to_ids = data_path + 'missing/missing_trimmed.list'
# path_to_ids = data_path + 'debug2.list'
# path_to_map_file = gatk_out_path + 'gvcf_force_active_no_opt.map'
# path_to_map_file = gatk_out_path + 'gvcf_bwa.map'
# path_to_map_file = gatk_out_path + 'combined_ids.map'
# path_to_map_file = gatk_out_path + 'combined_no_coll3.map'
# path_to_map_file = gatk_out_path + 'trimmed34.map'
path_to_map_file = gatk_out_path + 'missing.map'
# path_to_map_file = gatk_out_path + 'combined_ids_small.map'
# path_to_db = gatk_out_path + 'variant_force_active_no_opt.db'
# path_to_db = gatk_out_path + 'variant_bwa.db'
# path_to_db = gatk_out_path + 'variant_bwa_combined3.db'
# path_to_db = gatk_out_path + 'trimmed34.db'
path_to_db = gatk_out_path + 'mising.db'
# merged_vcf_path = gatk_out_path + 'merged_force_active_no_opt.vcf.gz'
# merged_vcf_path = gatk_out_path + 'merged_bwa.vcf.gz'
# merged_vcf_path = gatk_out_path + 'merged_bwa_combined3.vcf.gz'
# merged_vcf_path = gatk_out_path + 'merged_bwa_trimmed34.vcf.gz'
merged_vcf_path = gatk_out_path + 'merged_bwa_missing.vcf.gz'

thread_num = '32'

#
def gen_map_file():
    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    with open(path_to_map_file, 'w') as f:
        for sample_id in sample_ids:
            if exists(gvcf_path + sample_id + '_h37rv.g.vcf'):
                f.write(sample_id + '\t' + gvcf_path + sample_id + '_h37rv.g.vcf\n')


def run_GenomicsDBImport():
    if not exists(temp_path):
        makedirs(temp_path)
    system(path_to_GATK + ' --java-options "-Xmx4g -Xms4g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenomicsDBImport '
                          '--genomicsdb-workspace-path ' +
           path_to_db + ' --batch-size 1000 -L ch1:1-' + str(ref_len) + ' --sample-name-map ' + path_to_map_file +
           ' --tmp-dir=' + temp_path + ' --reader-threads ' + thread_num + ' > ' + log_path_db + ' 2>> ' + log_path_db)#


def run_GenotypeGVCFs():
    system(path_to_GATK + ' --java-options "-Xmx4g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenotypeGVCFs --sample-ploidy 1 -R ' + 
    path_to_h37rv + ' -V gendb://' + path_to_db + ' -O ' + merged_vcf_path + ' ' + '--tmp-dir=' + temp_path + ' > ' +
    log_path_genotype + ' 2>> ' + log_path_genotype)


if __name__ == '__main__':
    # gen_map_file()
    run_GenomicsDBImport()
    run_GenotypeGVCFs()
