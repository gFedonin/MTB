from os import listdir
from os.path import exists

from src.core.constants import data_path

gatk_out_path = data_path + 'snps/gatk_after_cortex/'
out_path_vcf = gatk_out_path + 'vcf/'
out_path_bam = gatk_out_path + 'bam/'
out_path_log = gatk_out_path + 'log/'
path_to_ids = data_path + 'all_with_pheno.txt'
path_to_bam = '/export/data/kchukreev/data/bam_files/gatk/realigned_bams_with_known/'#prepared_bams


def check_logs():
    for fname in listdir(out_path_log):
        for l in open(out_path_log + fname).readlines():
            if 'error' in l.lower():
                print(fname)
                break


def check_vcf():
    for l in open(path_to_ids).readlines():
        sample_id = l.strip()
        if not exists(out_path_vcf + sample_id + '_h37rv.g.vcf'):
            print(sample_id)


def check_if_bam_exists():
    for l in open(path_to_ids).readlines():
        sample_id = l.strip()
        if not exists(path_to_bam + sample_id + '_h37rv.bam'):
            print(sample_id)


if __name__ == '__main__':
    # check_vcf()
    # check_if_bam_exists()
    check_logs()