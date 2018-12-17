import os
from os.path import exists

from sklearn.externals.joblib import Parallel, delayed

from src.core.constants import data_path

path_to_list = data_path + 'coll18_unmapped.list'
path_to_fastq = data_path + 'coll18/'
path_to_smalt = '../../../smalt/bin/smalt'
path_to_index = data_path + 'h37rv'

out_path = data_path + 'coll18_bam/'


def map_file(id):
    os.system('%s map -c 0.8 -x -i 800 -n 10 %s %s%s_1.fastq.gz %s%s_2.fastq.gz | samtools view -Sb - > %s%s.bam' %
              (path_to_smalt, path_to_index, path_to_fastq, id, path_to_fastq, id, out_path, id))


if __name__ == '__main__':
    tasks = Parallel(n_jobs=15)(delayed(map_file)(l.strip()) for l in open(path_to_list, 'r').readlines()
                                if not exists(out_path + l.strip() + '.bam'))
    c = 0
    for task in tasks:
        c += 1
    print(str(c) + ' files processed')
    # for l in open(path_to_list, 'r').readlines():
    #     id = l.strip()
    #     os.system('%s map -c 0.8 -x -i 800 -n 144 %s %s%s_1.fastq.gz %s%s_2.fastq.gz | samtools view -Sb - > %s%s.bam' %
    #               (path_to_smalt, path_to_index, path_to_fastq, id, path_to_fastq, id, out_path, id))