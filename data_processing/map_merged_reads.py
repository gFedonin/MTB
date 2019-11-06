import os
from os import makedirs
from os.path import exists
from subprocess import check_call

from sklearn.externals.joblib import Parallel, delayed

from core.constants import data_path


path_to_smalt = '/export/home/fedonin/smalt/bin/'

# path_to_list = data_path + 'coll4.txt'
path_to_list = data_path + 'debug1.list'
# path_to_fastq = data_path + 'coll18_trimmed/'
# path_to_reads = data_path + '/tadpole/'
path_to_reads = data_path + '/pilon/'

path_to_index = data_path + 'h37rv_rev'

# out_path = data_path + 'coll18_bam/'
# out_path = data_path + 'bowtie_rev_l15n1/'
# out_path = data_path + 'bwa_mem_rev/'
# out_path = data_path + 'smalt_k10s1c0.9/'
# out_path = data_path + 'tadpole_minimap/'
out_path = data_path + 'pilon_minimap/'

# raw_suffix1 = '_p1.fastq.gz'
# raw_suffix2 = '_p2.fastq.gz'
# suffix = '_tadpole'
suffix = ''

thread_num = '32'


# def smalt(id):
#     check_call(path_to_smalt + 'smalt map -c 0.8 -x -i 800 -n ' + thread_num + ' ' + path_to_index + ' ' +
#                path_to_fastq + id + '/' + id + raw_suffix1 + ' ' + path_to_fastq + id + '/' + id + raw_suffix2 +
#                ' | samtools sort -@ ' + thread_num + ' - | samtools view -Sb -@ ' + thread_num + ' - > ' +
#                out_path + id + suffix + '.bam', shell=True)
#     check_call('samtools index -@ ' + thread_num + ' ' + out_path + id + suffix + '.bam', shell=True)


def smalt(id):
    check_call(path_to_smalt + 'smalt map -c 0.9 -x -n ' + thread_num + ' ' + path_to_index + ' ' +
               path_to_reads + id + '/' + id +
               '_tadpole.fasta | samtools sort -@ ' + thread_num + ' - | samtools view -Sb -@ ' + thread_num + ' - > ' +
               out_path + id + suffix + '.bam', shell=True)
    check_call('samtools index -@ ' + thread_num + ' ' + out_path + id + suffix + '.bam', shell=True)


def bwa_mem(id):
    check_call('bwa mem -t ' + thread_num + ' -R \'@RG\\tID:1\\tSM:' + id
               + '\' ' + path_to_index + ' ' + path_to_reads + id + '/' + id + '_tadpole.fasta | samtools sort -@ ' +
               thread_num + ' - | samtools view -bS -@ ' + thread_num + ' - > ' + out_path + id + suffix + '.bam',
               shell=True)
    check_call('samtools index -@ ' + thread_num + ' ' + out_path + id + suffix + '.bam', shell=True)


# def bowtie(id):
#     check_call('bowtie2 --mm --very-sensitive --rg-id ' + id + ' -p ' + thread_num + ' -x ' + path_to_index +
#                ' -U ' + path_to_reads + id + '/' + id + '_merged.fasta | samtools sort -@ ' + thread_num +
#                ' - | samtools view -bS -@ ' + thread_num + ' - > ' + out_path + id + '_h37rv.bam', shell=True)
#     check_call('samtools index -@ ' + thread_num + ' ' + out_path + id + suffix + '.bam', shell=True)


def minimap(id):
    check_call('/export/home/fedonin/minimap2-2.17_x64-linux/minimap2 -x sr -a -2 -t ' + thread_num +
               ' -R \'@RG\\tID:1\\tSM:' + id + '\' ' + data_path + 'h37rv.fasta ' + path_to_reads + id + '/' + 'pilon' +
               suffix + '.fasta | samtools sort -@ ' + thread_num +
               ' - | samtools view -bS -@ ' + thread_num + ' - > ' + out_path + id + suffix + '.bam', shell=True)
    check_call('samtools index -@ ' + thread_num + ' ' + out_path + id + suffix + '.bam', shell=True)


mapper = minimap


if __name__ == '__main__':
    if not exists(out_path):
        makedirs(out_path)
    tasks = Parallel(n_jobs=4)(delayed(mapper)(l.strip()) for l in open(path_to_list, 'r').readlines()
                               if not exists(out_path + l.strip() + '.bam'))
    c = 0
    for task in tasks:
        c += 1
    print(str(c) + ' samples processed')
    # for l in open(path_to_list, 'r').readlines():
    #     id = l.strip()
    #     os.system('%s map -c 0.8 -x -i 800 -n 144 %s %s%s_1.fastq.gz %s%s_2.fastq.gz | samtools view -Sb - > %s%s.bam' %
    #               (path_to_smalt, path_to_index, path_to_fastq, id, path_to_fastq, id, out_path, id))