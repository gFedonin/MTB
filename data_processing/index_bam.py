import os
from subprocess import check_call

from sklearn.externals.joblib import Parallel, delayed

from core.constants import data_path

path_to_bam = data_path + 'filtered_bams/'
path_to_sorted_bam = data_path + 'filtered_sorted_bams/'
path_to_ids = data_path + 'list0.txt'

thread_num = 32
sort = True

def index_file(sample_id):
    if sort:
        check_call('samtools sort ' + path_to_bam + sample_id +
                '_h37rv.bam > ' + path_to_sorted_bam + '_h37rv.bam')
        check_call('samtools index ' + path_to_sorted_bam + sample_id + '_h37rv.bam')
    else:
        check_call('samtools index ' + path_to_bam + sample_id + '_h37rv.bam')
    return 1


if __name__ == '__main__':
    tasks = Parallel(n_jobs=thread_num)(delayed(index_file)(l.strip()) for l in open(path_to_ids).readlines())
    c = 0
    for task in tasks:
        c += task
    print('%d files indexed' % c)
