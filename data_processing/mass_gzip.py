import os
from os.path import exists
from sklearn.externals.joblib import Parallel, delayed
from src.core.constants import data_path

path_to_list = data_path + 'coll18_unzipped.list'
path_to_fastq = data_path + 'coll18/'


def zip_file(id, index):
    os.system('gzip %s%s_%s.fastq' % (path_to_fastq, id, index))


if __name__ == '__main__':
    tasks = Parallel(n_jobs=-1)(delayed(zip_file)(l.strip(), '1') for l in open(path_to_list, 'r').readlines()
                                if not exists(path_to_fastq + l.strip() + '_1.fastq.gz'))
    c = 0
    for task in tasks:
        c += 1
    tasks = Parallel(n_jobs=-1)(delayed(zip_file)(l.strip(), '2') for l in open(path_to_list, 'r').readlines()
                                if not exists(path_to_fastq + l.strip() + '_2.fastq.gz'))
    for task in tasks:
        c += 1
    print(str(c) + ' files processed')
