from os import makedirs
from os.path import exists
from subprocess import check_call

from core.constants import data_path

path_to_disco = '/export/home/fedonin/Disco/'

path_to_ids = data_path + 'debug1.list'
path_to_reads = '/export/data/kkuleshov/myc/sra/'
out_path = data_path + 'disco/'


if __name__ == '__main__':
    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    for sample_id in sample_ids:
        if not exists(out_path + sample_id):
            makedirs(out_path + sample_id)
        check_call(path_to_disco + 'runAssembly.sh -in1 ' + path_to_reads + sample_id + '/' + sample_id +
                   '_R1.fastq.gz -in2 ' + path_to_reads + sample_id + '/' + sample_id +
                   '_R2.fastq.gz -d ' + out_path + sample_id, shell=True, cwd=out_path + sample_id)