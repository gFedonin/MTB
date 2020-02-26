from os import makedirs, system

from os.path import exists
from subprocess import call

from sklearn.externals.joblib import Parallel, delayed

from core.constants import data_path

data_set = 'walker18'
# path_to_bams = '/export/data/kchukreev/data/bam_files/gatk/realigned_bams_with_known/'
# path_to_bams = data_path + 'bwa_mem_rev_rm_dup/'
path_to_bams = data_path + data_set + '/bwa_mem_rev_rm_dup/'
# path_to_ids = data_path + 'Full_subset.txt'
# path_to_ids = data_path + 'all_with_pheno.txt'
# path_to_ids = data_path + data_set + '/' + data_set + '_trimmed.list'
path_to_ids = data_path + data_set + '/' + 'trim1234.list'
# out_path = data_path + 'coverages_rm_dup/'
out_path = data_path + data_set + '/coverages/'

suffix = ''


def genomecov(sample_id):
    call('genomeCoverageBed -d -ibam ' + path_to_bams + sample_id + suffix + '.bam > ' + out_path +
         sample_id + '.depth', shell=True)
    print(sample_id)
    return 1


if __name__ == '__main__':

    if not exists(out_path):
        makedirs(out_path)

    tasks = Parallel(n_jobs=-1)(delayed(genomecov)(line.strip())
                                for line in open(path_to_ids).readlines() if not exists(out_path + line.strip() + '.depth'))

    c = 0
    for task in tasks:
        c += task
    print(str(c))
