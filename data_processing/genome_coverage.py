from os import makedirs, system

from os.path import exists
from sklearn.externals.joblib import Parallel, delayed

from src.core.constants import data_path

path_to_bams = '/export/data/kchukreev/data/bam_files/gatk/realigned_bams_with_known/'
path_to_ids = data_path + 'Full_subset.txt'
out_path = data_path + 'coverages/'


def genomecov(sample_id):
    if exists(out_path + sample_id + '.depth'):
        return 0
    system('genomeCoverageBed -d -ibam ' + path_to_bams + sample_id + '_h37rv.bam > ' + out_path + sample_id + '.depth')
    print(sample_id)
    return 1


if __name__ == '__main__':

    if not exists(out_path):
        makedirs(out_path)

    tasks = Parallel(n_jobs=-1)(delayed(genomecov)(line.strip()) for line in open(path_to_ids).readlines())

    c = 0
    for task in tasks:
        c += task
    print(str(c))
