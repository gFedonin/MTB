from multiprocessing.pool import Pool
from os import makedirs, system
from os.path import exists

from src.core.constants import data_path

sample_list = data_path + 'coll18_supp.samples'
path_to_raw_data = data_path + 'coll18/'
out_path = data_path + 'coll18_fastqc/'

thread_num = 144


def fastqc(sample_id):
    if not exists(out_path + sample_id):
        makedirs(out_path + sample_id)
    system('~/FastQC/fastqc -q -t 2 --extract -o ' + out_path + sample_id + '/ ' + path_to_raw_data + sample_id +
           '_1.fastq.gz ' + path_to_raw_data + sample_id + '_2.fastq.gz > ' + out_path + sample_id + '/log.txt ' +
           '2>> ' + out_path + sample_id + '/log.txt')
    return 1


def run():
    if not exists(out_path):
        makedirs(out_path)
    with Pool(thread_num//2) as p:
        tasks = p.map_async(fastqc, [l.strip() for l in open(sample_list).readlines()])
        c = 0
        for task in tasks.get():
            c += task
        print(str(c))


def parse_summary():
    sample_ids = [l.strip() for l in open(sample_list).readlines()]
    with open(out_path + 'samples_with_adapters.list', 'w') as fout:
        for sample_id in sample_ids:
            adapter_found = False
            if exists(out_path + sample_id + '/' + sample_id + '_1_fastqc/summary.txt'):
                with open(out_path + sample_id + '/' + sample_id + '_1_fastqc/summary.txt') as f:
                    for line in f.readlines():
                        s = line.split('\t')
                        if s[1] == 'Adapter Content' and s[0] != 'PASS':
                            adapter_found = True
            if exists(out_path + sample_id + '/' + sample_id + '_2_fastqc/summary.txt'):
                with open(out_path + sample_id + '/' + sample_id + '_2_fastqc/summary.txt') as f:
                    for line in f.readlines():
                        s = line.split('\t')
                        if s[1] == 'Adapter Content' and s[0] != 'PASS':
                            adapter_found = True
            if adapter_found:
                fout.write(sample_id)
                fout.write('\n')


if __name__ == '__main__':
    parse_summary()