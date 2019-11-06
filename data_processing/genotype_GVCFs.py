from os import system, makedirs
from os.path import exists

from src.core.constants import data_path, ref_len

path_to_GATK = '/export/home/fedonin/gatk-4.1.2.0/gatk'
path_to_h37rv = data_path + 'h37rv.fasta'
# gatk_out_path = data_path + 'snps/gatk_before_cortex/'
gatk_out_path = data_path + 'snps/gatk_bowtie_rev_l15n1/'
gvcf_path = gatk_out_path + 'gvcf_force_active_no_opt/'
log_path_db = gatk_out_path + 'genomics_DB_import_force_active_no_opt.log'
log_path_genotype = gatk_out_path + 'genotype_force_active_no_opt.log'
temp_path = gatk_out_path + 'tmp/'
# path_to_ids = data_path + 'all_with_pheno.txt'
path_to_ids = data_path + 'debug2.list'
path_to_map_file = gatk_out_path + 'gvcf_force_active_no_opt.map'
path_to_db = gatk_out_path + 'variant_force_active_no_opt.db'
merged_vcf_path = gatk_out_path + 'merged_force_active_no_opt.vcf.gz'

thread_num = '32'


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
    gen_map_file()
    run_GenomicsDBImport()
    run_GenotypeGVCFs()
