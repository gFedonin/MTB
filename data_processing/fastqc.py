from multiprocessing.pool import Pool
from os import makedirs, system
from os.path import exists

from core.constants import data_path

# sample_list = '/export/home/fedonin/Enc/list.txt'
# sample_list = data_path + 'coll18_supp.samples'
sample_list = data_path + 'all_with_pheno.txt'
# path_to_raw_data = '/export/home/fedonin/Enc/data/raw/'
# path_to_raw_data = data_path + 'coll18/'
# path_to_raw_data = data_path + 'coll18_trimmed/'
# path_to_raw_data = '/export/data/kkuleshov/myc/sra/'
path_to_raw_data = data_path + 'kkuleshov_bbduk/'
# out_path = data_path + 'coll18_fastqc/'
# out_path = data_path + 'coll18_trimmed_fastqc/'
# out_path = data_path + 'kuleshov_fastqc/'
out_path = data_path + 'kuleshov_bbduk_fastqc/'
# out_path = '/export/home/fedonin/Enc/data/fastqc/'
# suffix1 = '_R1_001.fastq.gz'
# suffix1 = '_R1'
suffix1 = '_p1'
# suffix1 = '_1'
# suffix2 = '_R2_001.fastq.gz'
# suffix2 = '_R2'
suffix2 = '_p2'
# suffix2 = '_2'

thread_num = 32
overwrite = False


def fastqc(sample_id):
    if not exists(out_path + sample_id):
        makedirs(out_path + sample_id)
    if overwrite or not exists(out_path + sample_id + '/' + sample_id + suffix1 + '_fastqc/summary.txt') or \
            not exists(out_path + sample_id + '/' + sample_id + suffix2 + '_fastqc/summary.txt'):
        system('~/FastQC/fastqc -q -t 2 --extract -o ' + out_path + sample_id + '/ ' + path_to_raw_data + sample_id +
               '/' + sample_id +
               suffix1 + '.fastq.gz ' + path_to_raw_data + sample_id +
               '/' + sample_id +
               suffix2 + '.fastq.gz > ' + out_path + sample_id + '/log.txt ' +
               '2>> ' + out_path + sample_id + '/log.txt')
        return 1
    else:
        return 0


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
            overrepresented = set()
            sources = set()
            if exists(out_path + sample_id + '/' + sample_id + suffix1 + '_fastqc/summary.txt'):
                with open(out_path + sample_id + '/' + sample_id + suffix1 + '_fastqc/summary.txt') as f:
                    for line in f.readlines():
                        s = line.split('\t')
                        if s[1] in ('Adapter Content', 'Overrepresented sequences') and s[0] != 'PASS':
                            adapter_found = True
                            seq_list = False
                            with open(out_path + sample_id + '/' + sample_id + suffix1 + '_fastqc/fastqc_data.txt') as f1:
                                for l in f1.readlines():
                                    if l[0] == '#':
                                        continue
                                    if 'Overrepresented sequences' in l:
                                        seq_list = True
                                    elif 'END_MODULE' in l:
                                        seq_list = False
                                    elif seq_list:
                                        s = l.strip().split('\t')
                                        overrepresented.add(s[0])
                                        i = s[-1].find('(')
                                        sources.add(s[-1][0:i - 1])
            if exists(out_path + sample_id + '/' + sample_id + suffix2 + '_fastqc/summary.txt'):
                with open(out_path + sample_id + '/' + sample_id + suffix2 + '_fastqc/summary.txt') as f:
                    for line in f.readlines():
                        s = line.split('\t')
                        if s[1] in ('Adapter Content', 'Overrepresented sequences') and s[0] != 'PASS':
                            adapter_found = True
                            seq_list = False
                            with open(out_path + sample_id + '/' + sample_id + suffix2 + '_fastqc/fastqc_data.txt') as f1:
                                for l in f1.readlines():
                                    if l[0] == '#':
                                        continue
                                    if 'Overrepresented sequences' in l:
                                        seq_list = True
                                    elif 'END_MODULE' in l:
                                        seq_list = False
                                    elif seq_list:
                                        s = l.strip().split('\t')
                                        overrepresented.add(s[0])
                                        i = s[-1].find('(')
                                        sources.add(s[-1][0:i - 1])
            if adapter_found:
                fout.write(sample_id + '\t')
                fout.write(';'.join(overrepresented))
                fout.write('\t')
                fout.write(';'.join(sources))
                fout.write('\n')
    with open(out_path + 'samples_with_low_quality_warn_or_fail.list', 'w') as fout:
        for sample_id in sample_ids:
            bad_quality = False
            if exists(out_path + sample_id + '/' + sample_id + suffix1 + '_fastqc/summary.txt'):
                with open(out_path + sample_id + '/' + sample_id + suffix1 + '_fastqc/summary.txt') as f:
                    for line in f.readlines():
                        s = line.split('\t')
                        if s[1] == 'Per base sequence quality' and s[0] != 'PASS':
                            bad_quality = True
            if exists(out_path + sample_id + '/' + sample_id + suffix2 + '_fastqc/summary.txt'):
                with open(out_path + sample_id + '/' + sample_id + suffix2 + '_fastqc/summary.txt') as f:
                    for line in f.readlines():
                        s = line.split('\t')
                        if s[1] == 'Per base sequence quality' and s[0] != 'PASS':
                            bad_quality = True
            if bad_quality:
                fout.write(sample_id)
                fout.write('\n')
    with open(out_path + 'samples_with_low_quality_fail.list', 'w') as fout:
        for sample_id in sample_ids:
            bad_quality = False
            if exists(out_path + sample_id + '/' + sample_id + suffix1 + '_fastqc/summary.txt'):
                with open(out_path + sample_id + '/' + sample_id + suffix1 + '_fastqc/summary.txt') as f:
                    for line in f.readlines():
                        s = line.split('\t')
                        if s[1] == 'Per base sequence quality' and s[0] == 'FAIL':
                            bad_quality = True
            if exists(out_path + sample_id + '/' + sample_id + suffix2 + '_fastqc/summary.txt'):
                with open(out_path + sample_id + '/' + sample_id + suffix2 + '_fastqc/summary.txt') as f:
                    for line in f.readlines():
                        s = line.split('\t')
                        if s[1] == 'Per base sequence quality' and s[0] == 'FAIL':
                            bad_quality = True
            if bad_quality:
                fout.write(sample_id)
                fout.write('\n')


if __name__ == '__main__':
    run()
    parse_summary()
