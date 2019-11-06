from os import system, makedirs, listdir
from os.path import exists
from subprocess import check_call

from Bio import SeqIO
from sklearn.externals.joblib import Parallel, delayed

from core.constants import data_path
from core.data_reading import path_to_ref
from testing.skesa_to_bam import parse_delta, print_bam

path_to_spades = '/export/home/fedonin/SPAdes-3.13.1-Linux/bin/'
path_to_mummer = '/export/home/fedonin/mummer4/bin/'
path_to_samples = '/export/data/kkuleshov/myc/sra/'
# path_to_samples = data_path + 'coll18/coll18_bbduk_q30/'
# path_to_list = data_path + 'list5.txt'
path_to_list = data_path + 'all_with_pheno.txt'
# path_to_list = data_path + 'coll18/coll1.txt'
spades_out_path = data_path + 'spades/'
# out_path = data_path + 'coll18/spades/'
# log_path = data_path + 'coll18/spades_logs/'
mummer_out_path = data_path + 'spades_mummer/dnadiff/'
mummer_log_path = data_path + 'spades_mummer/log/'
bam_out_path = data_path + 'spades_bam/'


thread_num = '144'
max_memory = '32'


def run_spades():
    if not exists(spades_out_path):
        makedirs(spades_out_path)
    sample_ids = [l.strip() for l in open(path_to_list).readlines()]
    # tasks = Parallel(n_jobs=thread_num, batch_size=len(sample_ids)//thread_num + 1)(delayed(main)(
    #     ['-t', str(2*core_num//thread_num), '-m', str(max_memory), '--careful', '-1', path_to_samples + sample_id +
    #      '/' + sample_id + '_R1.fastq.gz', '-2', path_to_samples + sample_id + '/' +
    #     sample_id + '_R2.fastq.gz', '-o', out_path + sample_id]) for sample_id in sample_ids)
    # c = 0
    # for task in tasks:
    #     c += task
    # print('sum return %d' % c)
    for sample_id in sample_ids:
        check_call(path_to_spades + 'spades.py -t ' + thread_num + ' -m ' + max_memory + ' --careful -1 ' +
               path_to_samples + sample_id + '/' + sample_id + '_R1.fastq.gz -2 ' + path_to_samples + sample_id + '/' +
               sample_id + '_R2.fastq.gz -o ' + spades_out_path + sample_id, shell=True)


def run_mummer(sample_id):
    if not exists(mummer_out_path + sample_id):
        makedirs(mummer_out_path + sample_id)
    if not exists(mummer_out_path + sample_id + '/' + sample_id + '.snps'):
        if not exists(spades_out_path + sample_id + '/contigs.fasta'):
            return 0
        contigs = [l for l in open(spades_out_path + sample_id + '/contigs.fasta').readlines()]
        if len(contigs) == 0:
            print('no contigs for ' + sample_id)
            return 1
        status = check_call(path_to_mummer + 'dnadiff  --prefix=' + mummer_out_path + sample_id + '/' + sample_id + ' ' +
                        path_to_ref + ' ' + spades_out_path + sample_id + '/contigs.fasta > ' + mummer_log_path +
                        sample_id + '.dnadiff.log 2>> ' + mummer_log_path + sample_id + '.dnadiff.log', shell=True)
        if status != 0:
            print('mummer problem with ' + sample_id)
        return 1
    return 0


def mass_mummer():
    if not exists(mummer_out_path):
        makedirs(mummer_out_path)
    if not exists(mummer_log_path):
        makedirs(mummer_log_path)
    # if not exists(mummer_out_path_filtered):
    #     makedirs(mummer_out_path_filtered)
    tasks = Parallel(n_jobs=-1)(delayed(run_mummer)(sample_id) for sample_id in listdir(spades_out_path))
    c = 0
    for task in tasks:
        c += task
    print('%d samples processed' % c)


def parse_mummer(sample_id):
    if not exists(spades_out_path + sample_id + '/contigs.fasta'):
        return 0
    contigs = [record for record in SeqIO.parse(spades_out_path + sample_id + '/contigs.fasta', 'fasta')]
    if len(contigs) == 0:
        return 0
    contig_to_mapping = parse_delta(sample_id, mummer_out_path)
    print_bam(sample_id, contigs, contig_to_mapping, bam_out_path)
    return 1


def mass_mummer_to_bam():
    if not exists(bam_out_path):
        makedirs(bam_out_path)
    c = 0
    for sample_id in listdir(mummer_out_path):
        c += parse_mummer(sample_id)
    print("%d samples processed" % c)


if __name__ == '__main__':
    # run_spades()
    mass_mummer()
    mass_mummer_to_bam()
