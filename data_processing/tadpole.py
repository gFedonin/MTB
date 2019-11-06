from os import makedirs
from os.path import exists
from subprocess import check_call

from core.constants import data_path

path_to_bbmap = '/export/home/fedonin/bbmap/'

path_to_ids = data_path + 'debug1.list'
# path_to_ids = data_path + 'all_with_pheno.txt'
path_to_reads = '/export/data/kkuleshov/myc/sra/'
out_path = data_path + 'tadpole_extend100/'

thread_num = '32'


def assemble():
    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    for sample_id in sample_ids:
        if not exists(out_path + sample_id):
            makedirs(out_path + sample_id)
        not_ready = True
        # if exists(out_path + sample_id + '/' + sample_id + '_tadpole.fasta'):
        #     not_ready = open(out_path + sample_id + '/' + sample_id + '_tadpole.fasta').readline() == ''
        if not_ready and exists(path_to_reads + sample_id + '/' + sample_id + '_R1.fastq.gz'):
            check_call(path_to_bbmap + 'tadpole.sh k=93 mincountseed=1 mincontig=0 mincoverage=1 overwrite=true -t=' + thread_num + ' in=' + path_to_reads + sample_id + '/' + sample_id +
                       '_R1.fastq.gz,' + path_to_reads + sample_id + '/' + sample_id +
                       '_R2.fastq.gz out=' + out_path + sample_id + '/' + sample_id + '_tadpole.fasta > ' + out_path +
                       sample_id + '/' + sample_id + '_tadpole1.log 2> ' + out_path + sample_id + '/' + sample_id +
                       '_tadpole2.log', shell=True, cwd=out_path + sample_id)


def extend():
    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    for sample_id in sample_ids:
        if not exists(out_path + sample_id):
            makedirs(out_path + sample_id)
        not_ready = True
        # if exists(out_path + sample_id + '/' + sample_id + '_tadpole.fasta'):
        #     not_ready = open(out_path + sample_id + '/' + sample_id + '_tadpole.fasta').readline() == ''
        if not_ready and exists(path_to_reads + sample_id + '/' + sample_id + '_R1.fastq.gz'):
            check_call(path_to_bbmap + 'tadpole.sh -Xmx20g k=93 mode=extend el=100 er=100 overwrite=true -t=' + thread_num + ' in=' + path_to_reads + sample_id + '/' + sample_id +
                       '_R1.fastq.gz in2=' + path_to_reads + sample_id + '/' + sample_id +
                       '_R2.fastq.gz out=' + out_path + sample_id + '/' + sample_id + '_R1.fastq.gz out2=' + out_path +
                       sample_id + '/' + sample_id + '_R2.fastq.gz > ' + out_path +
                       sample_id + '/' + sample_id + '_tadpole1.log 2> ' + out_path + sample_id + '/' + sample_id +
                       '_tadpole2.log', shell=True, cwd=out_path + sample_id)


if __name__ == '__main__':
    extend()
