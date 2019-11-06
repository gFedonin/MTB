import os
from os.path import exists
from subprocess import check_call

from sklearn.externals.joblib import Parallel, delayed

from core.constants import data_path

# path_to_list = data_path + 'coll18_supp_not_in_our.txt'
# path_to_list = data_path + 'walker2018.list'
# path_to_list = data_path + 'Yang2017.list'
path_to_list = data_path + 'Farhat19.list'
path_to_fastq_dump = '/export/home/fedonin/sratoolkit.2.9.6-1-centos_linux64/bin/fasterq-dump.2.9.6'
# out_path = data_path + 'coll18/'
# out_path = data_path + 'walker18/'
# out_path = data_path + 'yang17/'
out_path = data_path + 'farhat19/'
temp_path = data_path + 'temp/'

thread_num = 3


def download_and_gzip(id):
    if not exists(out_path + id + '_1.fastq.gz') or not exists(out_path + id + '_2.fastq.gz'):
        # os.system('%s --split-files -O %s --gzip %s' % (path_to_fastq_dump, out_path, id))
        if not exists(out_path + id + '_1.fastq') or not exists(out_path + id + '_2.fastq'):
            try:
                check_call('%s -e 1 -O %s -t %s %s' % (path_to_fastq_dump, out_path, temp_path, id), shell=True)
            except:
                print('can\'t download ' + id)
                return 0
        check_call('gzip ' + out_path + id + '_1.fastq', shell=True)
        check_call('gzip ' + out_path + id + '_2.fastq', shell=True)
        return 1
    return 0


if __name__ == '__main__':
    if not exists(out_path):
        os.makedirs(out_path)
    if not exists(temp_path):
        os.makedirs(temp_path)
    samples = [l.strip() for l in open(path_to_list, 'r').readlines()]
    tasks = Parallel(n_jobs=min(thread_num, len(samples)), batch_size=len(samples) // thread_num + 1)(delayed(download_and_gzip)(sample)
                                                                                   for sample in samples)
    c = 0
    for task in tasks:
        c += task
    print('%d samples processed' % c)
    # for l in open(path_to_list, 'r').readlines():
    #     id = l.strip()
    #     if not exists(out_path + id + '_1.fastq.gz') and not exists(out_path + id + '_1.fastq'):
    #         # os.system('%s --split-files -O %s --gzip %s' % (path_to_fastq_dump, out_path, id))
    #         os.system('%s -e 32 -O %s %s' % (path_to_fastq_dump, out_path, id))