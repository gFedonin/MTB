import os
from os.path import exists

# from sklearn.externals.joblib import Parallel, delayed

# from src.core.constants import data_path
data_path = '/export/data/fedonin/MTB/data/'

# path_to_freebayes = '/export/home/fedonin/freebayes/bin/freebayes'
# path_to_freebayes = '/export/home/fedonin/freebayes1.1/bin/freebayes'
path_to_freebayes = '/export/home/fedonin/freebayes-v1.3.1'
# path_to_freebayes_scripts = '/export/home/fedonin/freebayes-1.3.1/scripts/'
path_to_vcflib = '/export/home/fedonin/freebayes/vcflib/'
path_to_bamaddrg = '/export/home/fedonin/bamaddrg/'
path_to_h37rv = data_path + 'h37rv.fasta'
# path_to_ids = data_path + 'list4.txt'
path_to_ids = data_path + 'all_with_pheno.txt'
# path_to_bam = data_path + 'filtered_bams/'
# path_to_bam = '/export/data/kchukreev/data/bam_files/gatk/prepared_bams/'
path_to_bam = '/export/data/kchukreev/data/bam_files/gatk/realigned_bams_with_known/'
out_path_independent = data_path + 'snps/vcfs_str_filtered/'
out_path_simultaneous = data_path + 'snps/freebayes_simultaneous_after_cortex/'
path_to_merged_vcfs = data_path + 'snps/freebayes_before_cortex_merged_vcf/'
path_to_merged_vcf = data_path + 'snps/freebayes_before_cortex_merged.vcf'
path_to_bam_file_list = data_path + 'merged_bams.list'
path_to_regions = data_path + 'snps/regions.list'
path_to_merged_bams = data_path + 'merged_bams_after_cortex/'
minAltFrac = '0.75'
minCoverage = '10'

thread_num = 144


def snp_call(sample_id):
    os.system(path_to_freebayes + ' --bam ' + path_to_bam + sample_id + '_h37rv.bam --vcf ' + out_path_independent +
              sample_id + '_h37rv.vcf --fasta-reference ' + path_to_h37rv + ' --pvar 0.0001 --ploidy 1 ' +
              '--min-mapping-quality 0 --min-base-quality 20 --min-alternate-fraction ' + minAltFrac +
              ' --min-coverage ' + minCoverage)
    return 1


def multisample_independent():
    if not exists(out_path_independent):
        os.makedirs(out_path_independent)
    tasks = Parallel(n_jobs=thread_num)(delayed(snp_call)(l.strip()) for l in open(path_to_ids).readlines())
    c = 0
    for task in tasks:
        c += task
    print('%d files processed' % c)


def gen_bam_file_list():
    sample_ids = [path_to_bam + l.strip() + '_h37rv.bam' for l in open(path_to_ids).readlines()
                  if exists(path_to_bam + l.strip() + '_h37rv.bam')]
    with open(path_to_bam_file_list, 'w') as f:
        f.write('\n'.join(sample_ids))


# def multisample_simultaneous():
#     if not exists(out_path_simultaneous):
#         os.makedirs(out_path_simultaneous)
#     os.system(gen_bamaddrg_str() + ' | ./freebayes-parallel <(./fasta_generate_regions.py ' + path_to_h37rv +
#               '.fai 10000) ' + str(thread_num) + ' -f ' + path_to_h37rv +
#               ' --stdin --pvar 0.0001 --ploidy 1 --min-mapping-quality 0 --min-base-quality 20 --min-alternate-fraction '
#               + minAltFrac + ' --min-coverage ' + minCoverage + ' > ' + out_path_simultaneous + 'all_with_pheno.vcf')


def snp_call_merge_bam(region):
    os.system('samtools merge -R ' + region + ' -u -b ' + path_to_bam_file_list + ' - | ' +
              path_to_freebayes + ' -f ' + path_to_h37rv +
              ' --stdin --pvar 0.0001 --ploidy 1 --min-mapping-quality 0 --min-base-quality 20 --min-alternate-fraction '
              + minAltFrac + ' --min-coverage ' + minCoverage + ' > ' + out_path_simultaneous + 'reg_' +
              region.replace(':', '_') + '.vcf')
    return 1


def multisample_simultaneous():
    if not exists(out_path_simultaneous):
        os.makedirs(out_path_simultaneous)
    tasks = Parallel(n_jobs=thread_num)(delayed(snp_call_merge_bam)(l.strip()) for l in open(path_to_regions).readlines())
    c = 0
    for task in tasks:
        c += task
    print('%d files processed' % c)


#
# def multisample_simultaneous_one_thread():
#     if not exists(out_path_simultaneous):
#         os.makedirs(out_path_simultaneous)
#     os.system(path_to_bamaddrg + 'bamaddrg -b < ' + path_to_bam_file_lits + ' | ' + path_to_freebayes + ' -f ' + path_to_h37rv +
#               ' --stdin --pvar 0.0001 --ploidy 1 --min-mapping-quality 0 --min-base-quality 20 --min-alternate-fraction '
#               + minAltFrac + ' --min-coverage ' + minCoverage + ' > ' + out_path_simultaneous + 'all_with_pheno.vcf')
#


def merge_batch_bamaddrg(i, step, sample_ids):
    if i + step < len(sample_ids):
        os.system(path_to_bamaddrg + 'bamaddrg -b ' + ' -b '.join(
            sample_ids[i: i + step]) + ' > ' + path_to_merged_bams + 'merged_' +
                  str(i) + '.bam')
    else:
        os.system(
            path_to_bamaddrg + 'bamaddrg -b ' + ' -b '.join(sample_ids[i:]) + ' > ' + path_to_merged_bams + 'merged_' +
            str(i) + '.bam')
    return 1


def merge_batch(i, step, sample_ids):
    if i + step < len(sample_ids):
        os.system('samtools merge -r ' + path_to_merged_bams + 'merged_' +
                  str(i) + '.bam ' + ' '.join(sample_ids[i: i + step]))
    else:
        os.system('samtools merge -r ' + path_to_merged_bams + 'merged_' +
                  str(i) + '.bam ' + ' '.join(sample_ids[i:]))
    return 1


def merge_file_groups(step):
    if not exists(path_to_merged_bams):
        os.makedirs(path_to_merged_bams)
    sample_ids = [path_to_bam + l.strip() + '_h37rv.bam' for l in open(path_to_ids).readlines()
                  if exists(path_to_bam + l.strip() + '_h37rv.bam')]
    tasks = Parallel(n_jobs=thread_num)(delayed(merge_batch)(i, step, sample_ids) for i in range(0, len(sample_ids), step))
    c = 0
    for task in tasks:
        c += task
    print('%d batches processed' % c)


def index_merged_bam_files():
    if not exists(path_to_merged_bams):
        os.makedirs(path_to_merged_bams)
    sample_ids = [fname for fname in os.listdir(path_to_merged_bams)]
    tasks = Parallel(n_jobs=thread_num)(delayed(os.system)('samtools index ' + path_to_merged_bams + fname)
                                        for fname in sample_ids if fname.endswith('.bam'))


def merge_vcf():
    # if not exists(path_to_merged_vcfs):
    #     os.makedirs(path_to_merged_vcfs)
    # file_list = [out_path_simultaneous + fname for fname in os.listdir(out_path_simultaneous)]
    # for i in range(len(file_list)//1000 + 1):
    #     os.system('cat ' + ' '.join(file_list[i*1000: min((i + 1)*1000, len(file_list))]) + ' | ' + path_to_vcflib +
    #               '/scripts/vcffirstheader | ' + path_to_vcflib + '/bin/vcfstreamsort -w 1000 | ' + path_to_vcflib +
    #               'bin/vcfuniq > ' + path_to_merged_vcfs + str(i) + '.vcf')
    file_list = [path_to_merged_vcfs + fname for fname in os.listdir(path_to_merged_vcfs)]
    os.system('cat ' + ' '.join(file_list) + ' | ' + path_to_vcflib + '/scripts/vcffirstheader | ' + path_to_vcflib +
              '/bin/vcfstreamsort -w 1000 | ' + path_to_vcflib + 'bin/vcfuniq > ' + path_to_merged_vcf)


if __name__ == '__main__':
    # gen_bam_file_list()
    # multisample_independent()
    # multisample_simultaneous()
    # gen_free_bayes_str()
    # print(gen_bamaddrg_str())
    # multisample_simultaneous_one_thread()
    # merge_file_groups(100)
    merge_vcf()
    # index_merged_bam_files()

