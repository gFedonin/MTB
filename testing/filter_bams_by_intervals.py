import os
from os import makedirs
from os.path import exists

from sklearn.externals.joblib import Parallel, delayed

from core.constants import data_path

path_to_list = data_path + 'all_with_pheno.txt'
path_to_bam = '/export/data/kchukreev/data/bam_files/gatk/prepared_bams/'
path_to_bed = data_path + 'repeats.bed'
out_path = data_path + 'prepared_bams_in_repeats/'


def run_samtools(sample_id):
    if not exists(out_path + sample_id + '_h37rv_repeats.bam'):
        os.system('samtools view -b -L ' + path_to_bed + ' ' + path_to_bam + sample_id + '_h37rv.bam > ' + out_path + sample_id + 
        '_h37rv_repeats.bam')
        return 1
    return 0


def filter_all():
    sample_ids = [l.strip() for l in open(path_to_list).readlines()]
    tasks = Parallel(n_jobs=-1)(delayed(run_samtools)(sample_id)
                                for sample_id in sample_ids)
    c = 0
    for task in tasks:
        c += task
    print('%d samples processed' % c)


if __name__ == "__main__":
    if not exists(out_path):
        makedirs(out_path)
    sample_ids = [l.strip() for l in open(path_to_list).readlines()]
    c = 0
    for sample_id in sample_ids:
        if not exists(out_path + sample_id + '_h37rv_repeats.bam'):
            os.system('samtools view -b -L ' + path_to_bed + ' ' + path_to_bam + sample_id + '_h37rv.bam > ' + out_path + sample_id +
                        '_h37rv_repeats.bam')
            c += 1
    print('%d samples processed' % c)
