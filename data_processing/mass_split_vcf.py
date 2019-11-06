from os import system, makedirs
from os.path import exists

from sklearn.externals.joblib import delayed, Parallel

from core.constants import data_path

path_to_GATK = '/export/home/fedonin/gatk-4.1.2.0/gatk'
path_to_h37rv = data_path + 'h37rv.fasta'
# out_path_gatk = data_path + 'snps/gatk_before_cortex/'
out_path_gatk = data_path + 'snps/gatk_bwa_mem_rev/'
path_to_merged_vcf = out_path_gatk + 'merged_force_active_no_opt.vcf'
# path_to_ids = out_path_gatk + 'sample.list'
path_to_ids = data_path + 'debug2.list'
out_path_vcf = out_path_gatk + 'vcf_split_gatk_fa_no_opt/'
# out_path_vcf_filtered = out_path_gatk + 'vcf_filtered/'
# out_path_log = out_path_gatk + 'filter_log/'

minAltFrac = '0.75'
minCoverage = '10'

force_filter = True


def extract_sample(sample_id, out_path_vcf, big_vcf):
    if not exists(out_path_vcf + sample_id + '.vcf.gz'):
        system('bcftools view -q ' + minAltFrac + ' -O z -s ' + sample_id + ' ' + big_vcf + ' > ' + out_path_vcf + sample_id + '.vcf.gz')
        return 1
    return 0


def extract_sample_gatk(sample_id, out_path_vcf, big_vcf):
    if not exists(out_path_vcf + sample_id + '.vcf'):
        system(path_to_GATK + ' --java-options "-Xmx1g -Xms1g" SelectVariants --exclude-non-variants true ' +
                              '--remove-unused-alternates true -OVI false -R ' + path_to_h37rv + ' -sn ' +
               sample_id + ' -V ' + big_vcf + ' -O ' + out_path_vcf + sample_id + '.vcf')
        return 1
    return 0


# def filter_sample(sample_id):
#     if not exists(out_path_vcf_filtered + sample_id + '.vcf') or force_filter:
#         system(path_to_GATK + ' --java-options "-Xmx1g -Xms1g" VariantFiltration -R ' + path_to_h37rv + ' -V ' + out_path_vcf + sample_id + '.vcf.gz' +
#             ' --filter-name "QD" --filter-expression "QD < 2.0" --filter-name "FS" --filter-expression "FS > 60.0" \
#             --filter-name "MQ" --filter-expression "MQ < 40.0" --filter-name "MQRankSum"\
#             --filter-expression "MQRankSum < -12.5" --filter-name "ReadPosRankSum"\
#             --filter-expression "ReadPosRankSum < -8.0" --filter-name "SOR" ' +
#             '--filter-expression "SOR > 3.0" -O ' + out_path_vcf_filtered + sample_id + '.vcf > ' + out_path_log +
#             sample_id + '.log 2>> ' + out_path_log + sample_id + '.log')
#         return 1
#     else:
#         return 0


def split_gatk():
    if not exists(out_path_vcf):
        makedirs(out_path_vcf)
    # if not exists(out_path_vcf_filtered):
    #     makedirs(out_path_vcf_filtered)
    # if not exists(out_path_log):
    #     makedirs(out_path_log)
    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    tasks = Parallel(n_jobs=10)(delayed(extract_sample_gatk)(sample_id, out_path_vcf,  path_to_merged_vcf)
                                for sample_id in sample_ids)
    c = 0
    for task in tasks:
        c += task
    # for sample_id in sample_ids:
    #     c += extract_sample_gatk(sample_id, out_path_vcf,  out_path_gatk + 'merged.vcf.gz')
    print('%d samples extracted' % c)
    # tasks = Parallel(n_jobs=-1)(delayed(filter_sample)(sample_id)
    #                             for sample_id in sample_ids)
    # c = 0
    # for task in tasks:
    #     c += 1
    # print('%d samples filtered' % c)


def split_freebayes():
    path_to_big_vcf = data_path + 'snps/freebayes_after_cortex_merged.vcf'
    out_path_bayes = data_path + 'snps/freebayes_after_cortex_vcf/'
    if not exists(out_path_bayes):
        makedirs(out_path_bayes)
    sample_ids = [l.strip() for l in open(data_path + 'snps/freebayes_before_sample.list').readlines()]
    tasks = Parallel(n_jobs=-1)(delayed(extract_sample)(sample_id, out_path_bayes, path_to_big_vcf)
                                for sample_id in sample_ids)
    c = 0
    for task in tasks:
        c += 1
    print('%d samples extracted' % c)


if __name__ == '__main__':
    split_gatk()
    # split_freebayes()
