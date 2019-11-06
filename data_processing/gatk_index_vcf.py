from os import system, makedirs, listdir
from os.path import exists

from sklearn.externals.joblib import delayed, Parallel

from core.constants import data_path

path_to_GATK = '/export/home/fedonin/gatk-4.1.2.0/gatk'
path_to_ids = data_path + 'all_with_pheno.txt'#'gatk_problem.list'
out_path = data_path + 'snps/gatk_before_cortex/'
# out_path_vcf = out_path + 'vcf_split/'
out_path_vcf = out_path + 'gvcf_force_active/'
# out_path_log = out_path + 'index_split_log/'
out_path_log = out_path + 'index_gvcf_force_active_log/'


def index_gvcf_file(sample_id):
    if not exists(out_path_vcf + sample_id + '_h37rv.g.vcf.idx'):
        log_path = out_path_log + sample_id + '.log'
        system(path_to_GATK + ' --java-options "-Xmx1g -Xms1g" IndexFeatureFile -F ' + out_path_vcf + sample_id + '_h37rv.g.vcf > ' +
               log_path + ' 2>> ' + log_path)
        return 1
    else:
        return 0


def index_vcfgz_file(sample_id):
    if not exists(out_path_vcf + sample_id + '.vcf.gz.tbi'):
        log_path = out_path_log + sample_id + '.log'
        if exists(out_path_vcf + sample_id + '.vcf.gz'):
            system(path_to_GATK + ' --java-options "-Xmx1g -Xms1g" IndexFeatureFile -F ' + out_path_vcf + sample_id +
                   '.vcf.gz > ' + log_path + ' 2>> ' + log_path)
            return 1
        else:
            return 0
    else:
        return 0


def index_all_gvcf_files():
    if not exists(out_path_log):
        makedirs(out_path_log)
    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    tasks = Parallel(n_jobs=-1)(delayed(index_gvcf_file)(sample_id) for sample_id in sample_ids)
    c = 0
    for task in tasks:
        c += task
    print('%d files processed' % c)


def index_all_vcf_gz_files():
    if not exists(out_path_log):
        makedirs(out_path_log)
    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    tasks = Parallel(n_jobs=-1)(delayed(index_vcfgz_file)(sample_id) for sample_id in sample_ids)
    c = 0
    for task in tasks:
        c += task
    print('%d files processed' % c)


def parse_logs():
    for fname in listdir(out_path_log):
        id = fname[:fname.rfind('.')]
        for l in open(out_path_log + fname).readlines():
            if 'A USER ERROR has occurred' in l:
                print(id)
                break


if __name__ == '__main__':
    # parse_logs()
    index_all_gvcf_files()
    # index_all_vcf_gz_files()