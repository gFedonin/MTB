import os
from os.path import exists

from src.core.constants import data_path

path_to_list = data_path + 'coll18_supp_not_in_our.txt'
path_to_fastq_dump = '/export/home/fedonin/sratoolkit.2.9.1-1-centos_linux64/bin/fasterq-dump.2.9.1'
out_path = data_path + 'coll18/'


if __name__ == '__main__':
    for l in open(path_to_list, 'r').readlines():
        id = l.strip()
        if not exists(out_path + id + '_1.fastq.gz') and not exists(out_path + id + '_1.fastq'):
            # os.system('%s --split-files -O %s --gzip %s' % (path_to_fastq_dump, out_path, id))
            os.system('%s -e 32 -O %s %s' % (path_to_fastq_dump, out_path, id))