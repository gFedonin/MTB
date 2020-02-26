import os
from os.path import exists
from subprocess import call

from sklearn.externals.joblib import Parallel, delayed

from core.constants import data_path

data_set = 'patric'
path_to_list = data_path + data_set + '/patric_new.list'
# path_to_list = data_path + '/problem_for_download.list'
path_to_fastq_dump = '/export/home/fedonin/sratoolkit.2.9.6-1-centos_linux64/bin/fasterq-dump.2.9.6'
# out_path = data_path + 'coll18/'
out_path = data_path + data_set + '/' + data_set + '_raw/'
# out_path = data_path + 'problem_raw/'
temp_path = data_path + 'temp/'

thread_num = 3


def download_and_gzip(id):
    # os.system('%s --split-files -O %s --gzip %s' % (path_to_fastq_dump, out_path, id))
    if not exists(out_path + id + '_1.fastq') or not exists(out_path + id + '_2.fastq'):
        try:
            call('%s -e 1 -O %s -t %s %s' % (path_to_fastq_dump, out_path, temp_path, id), shell=True)
        except:
            print('can\'t download ' + id)
            return 0
    call('pigz ' + out_path + id + '_1.fastq', shell=True)
    if exists(out_path + id + '_1.fastq.gz') and exists(out_path + id + '_1.fastq'):
        call('rm ' + out_path + id + '_1.fastq', shell=True)
    call('pigz ' + out_path + id + '_2.fastq', shell=True)
    if exists(out_path + id + '_2.fastq.gz') and exists(out_path + id + '_2.fastq'):
        call('rm ' + out_path + id + '_2.fastq', shell=True)
    return 1


if __name__ == '__main__':
    if not exists(out_path):
        os.makedirs(out_path)
    if not exists(temp_path):
        os.makedirs(temp_path)
    samples = [l.strip() for l in open(path_to_list, 'r').readlines()]
    tasks = Parallel(n_jobs=min(thread_num, len(samples)), batch_size=len(samples) // thread_num + 1)(
        delayed(download_and_gzip)(sample) for sample in samples
                        if not exists(out_path + sample + '_1.fastq.gz') or not exists(out_path + sample + '_2.fastq.gz'))
    c = 0
    for task in tasks:
        c += task
    print('%d samples processed' % c)
    # for l in open(path_to_list, 'r').readlines():
    #     id = l.strip()
    #     if not exists(out_path + id + '_1.fastq.gz') and not exists(out_path + id + '_1.fastq'):
    #         # os.system('%s --split-files -O %s --gzip %s' % (path_to_fastq_dump, out_path, id))
    #         os.system('%s -e 32 -O %s %s' % (path_to_fastq_dump, out_path, id))