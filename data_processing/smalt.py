import os
from os.path import exists

from sklearn.externals.joblib import Parallel, delayed

from core.constants import data_path

path_to_list = data_path + 'coll4.txt'
path_to_fastq = data_path + 'coll18_trimmed/'
path_to_smalt = '/export/home/fedonin/smalt/bin/smalt'
path_to_index = data_path + 'h37rv'

out_path = data_path + 'coll18_bam/'


def map_file(id):
    os.system(path_to_smalt + ' map -c 0.8 -x -i 800 -n 10 ' + path_to_index + ' ' + path_to_fastq + id + '/' + id +
              '_p1.fastq.gz ' + path_to_fastq + id + '/' + id +'_p2.fastq.gz | samtools view -Sb - > ' + out_path +
              id + '.bam')


if __name__ == '__main__':
    if not exists(out_path):
        os.makedirs(out_path)
    tasks = Parallel(n_jobs=4)(delayed(map_file)(l.strip()) for l in open(path_to_list, 'r').readlines()
                                if not exists(out_path + l.strip() + '.bam'))
    c = 0
    for task in tasks:
        c += 1
    print(str(c) + ' files processed')
    # for l in open(path_to_list, 'r').readlines():
    #     id = l.strip()
    #     os.system('%s map -c 0.8 -x -i 800 -n 144 %s %s%s_1.fastq.gz %s%s_2.fastq.gz | samtools view -Sb - > %s%s.bam' %
    #               (path_to_smalt, path_to_index, path_to_fastq, id, path_to_fastq, id, out_path, id))