import os
from os import makedirs
from os.path import exists
from subprocess import check_call

from sklearn.externals.joblib import Parallel, delayed
from core.constants import data_path

path_to_GATK = '/export/home/fedonin/gatk-4.1.4.0/gatk'
path_to_piccard = '/export/home/fedonin/picard.jar'
path_to_h37rv = data_path + 'h37rv.fasta'
# 'gatk_problem.list'
# path_to_ids = data_path + 'all_with_pheno.txt'
# path_to_ids = data_path + 'debug1.list'
# path_to_ids = data_path + 'list5.txt'
path_to_ids = data_path + 'coll18/coll5.txt'
# path_to_ids = data_path + 'coll18/coll18_supp.samples'
# path_to_ids = '/export/home/fedonin/MTB/src/testing/errors_before.list'
# path_to_prepared_bam = '/export/data/kchukreev/data/bam_files/gatk/prepared_bams/'
# path_to_prepared_bam = data_path + 'bwa_mem_rev_rm_dup/'
path_to_prepared_bam = data_path + 'coll18/bwa_mem_rev_rm_dup/'
# path_to_bam = data_path + 'coll18/coll18_trimmed_bam/'
path_to_bam = data_path + 'coll18/bwa_mem_rev/'
# path_to_bam = data_path + 'bwa_mem_rev/'
# path_to_prepared_bam = data_path + 'coll18/coll18_prepared_bam/'
# path_to_duplication_metrix = data_path + 'coll18/coll18_metrix/'
path_to_duplication_metrix = data_path + 'coll18/bwa_mem_rev__metrix/'
# path_to_duplication_metrix = data_path + 'bwa_mem_rev_metrix/'
path_to_prep_logs = data_path + 'coll18/bwa_mem_rev_prep_logs/'
# path_to_prep_logs = data_path + 'bwa_mem_rev_logs/'
# path_to_bam = '/export/data/kchukreev/data/bam_files/gatk/realigned_bams_with_known/'
# out_path = data_path + 'snps/gatk_before_cortex/'
out_path = data_path + 'snps/gatk_before_cortex_debug4/'
# out_path = data_path + 'snps/gatk_bowtie_rev_l15n1/'
# out_path = data_path + 'coll18/snps/gatk/'
# out_path = data_path + 'snps/gatk_after_cortex/'

suffix = '_h37rv'
# suffix = ''

gvcf = True


if gvcf:
    out_path_vcf = out_path + 'gvcf_force_active_no_opt/'
    out_path_bam = out_path + 'gvcf_bam_force_active_no_opt/'
    out_path_log = out_path + 'gvcf_log_force_active_no_opt/'
    out_path_activity = out_path + 'gvcf_activity_force_active_no_opt/'
else:
    out_path_vcf = out_path + 'vcf/'
    out_path_bam = out_path + 'bam/'
    out_path_log = out_path + 'log/'


thread_num = 144


def snp_call_gvcf(sample_id):
    if exists(path_to_prepared_bam + sample_id + suffix + '.bam'):
        if not exists(out_path_vcf + sample_id + '_h37rv.g.vcf'):
            cln = path_to_GATK + ' --java-options "-Xmx4g" HaplotypeCaller -VS SILENT -R ' + path_to_h37rv + ' -I ' + \
                           path_to_prepared_bam + sample_id + suffix + '.bam -O ' + out_path_vcf + sample_id + \
                           '_h37rv.g.vcf -ERC GVCF -bamout ' + out_path_bam + sample_id + \
                           '_h37rv_gatk.bam --smith-waterman FASTEST_AVAILABLE --force-active true ' \
                           '--disable-optimizations true -ploidy 1 ' + \
                  ' --activity-profile-out ' + out_path_activity + sample_id + \
                           '.activity_profile > ' + out_path_log + sample_id + '.log 2>> ' + out_path_log + sample_id + '.log'
            print(cln)
            os.system(cln)  # --native-pair-hmm-threads 4
            return 1
        else:
            return 0
    else:
        print(path_to_prepared_bam + sample_id + suffix + '.bam' + ' not found')
        return 0


# def snp_call_gvcf(sample_id):
#     if exists(path_to_prepared_bam + sample_id + suffix + '.bam'):
#         if not exists(out_path_vcf + sample_id + '_h37rv.g.vcf'):
#             cln = path_to_GATK + ' --java-options "-Xmx4g" HaplotypeCaller --debug-graph-transformations -VS SILENT -R ' + \
#                   path_to_h37rv + ' -I ' + \
#                            path_to_prepared_bam + sample_id + suffix + '.bam -O ' + out_path_vcf + sample_id + \
#                            '_h37rv.g.vcf -ERC GVCF -bamout ' + out_path_bam + sample_id + \
#                            '_h37rv_gatk.bam --smith-waterman FASTEST_AVAILABLE --force-active true ' \
#                            '--disable-optimizations true -ploidy 1 ' + \
#                             ' --activity-profile-out ' + out_path_activity + sample_id + \
#                            '.activity_profile --graph-output ' + out_path + sample_id + '.dot > ' + out_path_log + \
#                   sample_id + '.log 2>> ' + out_path_log + sample_id + '.log'
#             # print(cln)
#             os.system(cln)  # --native-pair-hmm-threads 4
#             return 1
#         else:
#             return 0
#     else:
#         print(sample_id + ' not found')
#         return 0


# def snp_call_gvcf(sample_id):
#     if exists(path_to_prepared_bam + sample_id + suffix + '.bam'):
#         if not exists(out_path_vcf + sample_id + '_h37rv.g.vcf'):
#             cln = path_to_GATK + ' --java-options "-Xmx4g" HaplotypeCaller --debug-graph-transformations -VS SILENT -R ' + \
#                   path_to_h37rv + ' -I ' + \
#                            path_to_prepared_bam + sample_id + suffix + '.bam -O ' + out_path_vcf + sample_id + \
#                            '_h37rv.g.vcf -ERC GVCF -bamout ' + out_path_bam + sample_id + \
#                            '_h37rv_gatk.bam --smith-waterman FASTEST_AVAILABLE --force-active true ' \
#                            '-L ch1:2248000-2249000 --disable-optimizations true -ploidy 1 ' + \
#                             ' --activity-profile-out ' + out_path_activity + sample_id + \
#                            '.activity_profile > ' + out_path_log + \
#                   sample_id + '.log 2>> ' + out_path_log + sample_id + '.log'
#             # print(cln)
#             os.system(cln)  # --native-pair-hmm-threads 4
#             return 1
#         else:
#             return 0
#     else:
#         print(sample_id + ' not found')
#         return 0


def snp_call(sample_id):
    if exists(path_to_prepared_bam + sample_id + suffix + '.bam'):
        if not exists(out_path_vcf + sample_id + '_h37rv.vcf'):
            os.system(path_to_GATK + ' --java-options "-Xmx1g" HaplotypeCaller -R ' + path_to_h37rv + ' -I ' +
                      path_to_prepared_bam + sample_id + suffix + '.bam -O ' + out_path_vcf + sample_id +
                      '_h37rv.vcf -bamout ' +
                      out_path_bam + sample_id + '_h37rv_gatk.bam --create-output-variant-index false ' +
                      '--smith-waterman FASTEST_AVAILABLE ' +
                      '-ploidy 1 ' + ' > ' + out_path_log + sample_id
                      + '.log 2>> ' + out_path_log + sample_id + '.log')  # --native-pair-hmm-threads 1
            return 1
        else:

            return 0
    else:
        print(sample_id + ' not found')
        return 0


def call_snps():
    if not exists(out_path):
        os.makedirs(out_path)
    if not exists(out_path_vcf):
        os.makedirs(out_path_vcf)
    if not exists(out_path_bam):
        os.makedirs(out_path_bam)
    if not exists(out_path_log):
        os.makedirs(out_path_log)
    if not exists(out_path_activity):
        os.makedirs(out_path_activity)
    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    print('total ids read %d' % len(sample_ids))
    if gvcf:
        tasks = Parallel(n_jobs=min(thread_num, len(sample_ids)))(delayed(snp_call_gvcf)(sample_id)
                                                                  for sample_id in sample_ids)
    else:
        tasks = Parallel(n_jobs=min(thread_num, len(sample_ids)))(delayed(snp_call)(sample_id)
                                                                  for sample_id in sample_ids)
    c = 0
    for task in tasks:
        c += task
    print('%d files processed' % c)


def prepare_bam(sample_id):
    if not exists(path_to_bam + sample_id + suffix + '.bam'):
        print('no such sample: ' + sample_id)
        return 0
    c = 'java -Xmx1g -jar ' + path_to_piccard + ' AddOrReplaceReadGroups I=' + path_to_bam + sample_id + suffix + '.bam O=' + \
        path_to_prepared_bam + sample_id + '.tm.bam' + \
        ' RGPL=illumina RGLB=LaneX RGPU=NONE SORT_ORDER=coordinate RGSM=' + sample_id + ' > ' + path_to_prep_logs + sample_id + '.log'
    os.system(c)  # , shell=True)
    c = 'java -Xmx4g -jar ' + path_to_piccard + ' MarkDuplicates INPUT=' + path_to_prepared_bam + sample_id + '.tm.bam' + ' OUTPUT=' + \
        path_to_prepared_bam + sample_id + suffix + '.bam' + \
        ' METRICS_FILE=' + path_to_duplication_metrix + sample_id + \
        '.metrix REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT' + \
        ' > ' + path_to_prep_logs + sample_id + '.log'
    os.system(c)
    os.system('rm ' + path_to_prepared_bam + sample_id + '.tm.bam ')
    os.system('samtools index ' + path_to_prepared_bam + sample_id + suffix + '.bam')
    return 1


def prepare_all_bams():
    if not exists(path_to_prepared_bam):
        os.makedirs(path_to_prepared_bam)
    if not exists(path_to_duplication_metrix):
        os.makedirs(path_to_duplication_metrix)
    if not exists(path_to_prep_logs):
        os.makedirs(path_to_prep_logs)
    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    print('total ids read %d' % len(sample_ids))
    tasks = Parallel(n_jobs=min(thread_num, len(sample_ids)))(
        delayed(prepare_bam)(sample_id) for sample_id in sample_ids)
    c = 0
    for task in tasks:
        c += task
    print('%d files processed' % c)
    # c = 0
    # for sample_id in sample_ids:
    #     c += prepare_bam(sample_id)
    # print('%d files processed' % c)


def mark_duplicates(sample_id):
    if not exists(path_to_bam + sample_id + suffix + '.bam'):
        print('no such sample: ' + sample_id)
        return 0
    c = 'java -Xmx4g -jar ' + path_to_piccard + ' MarkDuplicates INPUT=' + path_to_bam + sample_id + suffix + \
        '.bam OUTPUT=' + \
        path_to_prepared_bam + sample_id + '.bam' + \
        ' METRICS_FILE=' + path_to_duplication_metrix + sample_id + \
        '.metrix REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT' + \
        ' > ' + path_to_prep_logs + sample_id + '.log 2>> ' + path_to_prep_logs + sample_id + '.log'
    check_call(c, shell=True)
    check_call('samtools index ' + path_to_prepared_bam + sample_id + '.bam', shell=True)
    return 1


def mark_dup_all_bams():
    if not exists(path_to_prepared_bam):
        os.makedirs(path_to_prepared_bam)
    if not exists(path_to_duplication_metrix):
        os.makedirs(path_to_duplication_metrix)
    if not exists(path_to_prep_logs):
        os.makedirs(path_to_prep_logs)
    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    print('total ids read %d' % len(sample_ids))
    tasks = Parallel(n_jobs=min(thread_num, len(sample_ids)))(
        delayed(mark_duplicates)(sample_id) for sample_id in sample_ids)
    c = 0
    for task in tasks:
        c += task
    print('%d files processed' % c)
    # c = 0
    # for sample_id in sample_ids:
    #     c += prepare_bam(sample_id)
    # print('%d files processed' % c)


if __name__ == '__main__':
    # prepare_all_bams()
    # call_snps()
    mark_dup_all_bams()