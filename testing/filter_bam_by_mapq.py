from os import makedirs
from os.path import exists
from subprocess import check_call

from numpy import mean, std
from sklearn.externals.joblib import Parallel, delayed

import pysam as ps
from core.constants import data_path

# path_to_ids = data_path + 'debug1.list'
path_to_ids = data_path + 'all_with_pheno.txt'
# path_to_bam = '/export/data/kchukreev/data/bam_files/gatk/realigned_bams_with_known/'
# path_to_bam = data_path + 'minimap/'
# path_to_bam = data_path + 'bowtie_rev_l15n1/'
# path_to_bam = data_path + 'smalt_k10s1c0.9/'
# path_to_bam = data_path + 'bwa_mem_rev/'
# path_to_bam = data_path + 'skesa_minimap/'
path_to_bam = data_path + 'skesa_bwa_mem/'
# out_path = data_path + 'realigned_bams_with_known_mapq/'
# out_path = data_path + 'minimap_mapq/'
# out_path = data_path + 'bowtie_rev_l15n1_mapq/'
# out_path = data_path + 'smalt_k10s1c0.9_mapq/'
# out_path = data_path + 'bwa_mem_rev_mapq/'
# out_path = data_path + 'skesa_minimap_mapq/'
out_path = data_path + 'skesa_bwa_mapq/'

std_mult = 2
# suffix = '_h37rv'
suffix_in = ''
suffix_out = ''


def filter_bam(sample_id):
    samfile = ps.AlignmentFile(path_to_bam + sample_id + suffix_in + '.bam', "rb")
    # samfile = ps.AlignmentFile(path_to_bam + sample_id + '/' + sample_id + '_notmerged_minimap.bam', "rb")
    mapq = []
    alns = []
    for read in samfile.fetch():
        alns.append(read)
        mapq.append(read.mapping_quality)
    mean_mapq = mean(mapq)
    std_mapq = std(mapq)
    with ps.AlignmentFile(out_path + sample_id + suffix_out + '.bam', "wb", header=samfile.header) as outf:
        for read in alns:
            if read.mapping_quality >= mean_mapq - std_mult*std_mapq:
                outf.write(read)
    check_call('samtools index ' + out_path + sample_id + suffix_out + '.bam', shell=True)
    return 1


if __name__ == '__main__':
    if not exists(out_path):
        makedirs(out_path)
    tasks = Parallel(n_jobs=-1)(delayed(filter_bam)(sample_id.strip()) for sample_id in open(path_to_ids).readlines())
    c = 0
    for task in tasks:
        c += task
    print('%d bam files filtered\n' % c)
