from os import listdir, makedirs, system
from os.path import exists
from subprocess import check_call

from core.constants import data_path

path_to_bam_files = '/export/data/kchukreev/data/bam_files/gatk/realigned_bams_with_known/'  # prepared_bams
# path_to_bam_files = data_path + 'filtered_bams/'
# out_path = data_path + 'realigned_all_ppe_pgrs/'
out_path = data_path + 'realigned_8141_8167/'
# path_to_intervals = data_path + 'ppe_pgrs.bed'
path_to_intervals = data_path + 'Spurious.bed'
path_to_list = data_path + 'list0.txt'


if __name__ == '__main__':
    if not exists(out_path):
        makedirs(out_path)
    # sample_ids = [l.strip() for l in open(path_to_list).readlines()]
    sample_ids = ['SAMN03648890', 'SAMN03648892', 'SAMN03648885']
    for sample_id in sample_ids:
        if not exists(out_path + sample_id + '.bam'):
            # print('samtools view --threads 32 -b -L ' + path_to_intervals + ' ' + path_to_bam_files + sample_id +
            #        '_h37rv.bam > ' + out_path + sample_id + '.bam')
            system('samtools view --threads 32 -b -L ' + path_to_intervals + ' ' + path_to_bam_files + sample_id +
                   '_h37rv.bam > ' + out_path + sample_id + '.bam')
            system('samtools index ' + out_path + sample_id + '.bam')
