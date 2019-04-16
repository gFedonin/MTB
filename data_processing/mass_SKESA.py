from os import system

from sklearn.externals.joblib import Parallel, delayed

from src.core.constants import data_path

path_to_skesa = '../skesa.centos6.10'
path_to_samples = '/export/data/kkuleshov/myc/sra/'
path_to_list = data_path + 'problem_samples.list'
out_path = data_path + 'skesa/'
log_path = data_path + 'skesa_log/'


def run(sample_id):
    system(path_to_skesa + ' --fastq ' + path_to_samples + sample_id + '/' + sample_id + '_R1.fastq.gz,' +
           path_to_samples + sample_id + '/' + sample_id + '_R2.fastq.gz --memory 16 > ' + out_path + sample_id +
           '.skesa.fa 2> ' + log_path + sample_id + '.skesa.log')
    return 1


if __name__ == '__main__':
    tasks = Parallel(n_jobs=-1)(delayed(run)(l.strip()) for l in open(path_to_list).readlines())