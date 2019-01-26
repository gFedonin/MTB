import pysam as ps
from scipy.stats import ttest_ind
from sklearn.externals.joblib import Parallel, delayed

from src.core.annotations import read_annotations
from src.core.constants import upstream_length

path_to_ids = '/export/data/fedonin/MTB/data/dr_covered_with_pheno_and_snp_new.txt'
path_to_bam = '/export/data/kchukreev/data/bam_files/gatk/realigned_bams_with_known/'


gene_filter = 'PGRS'

out_path = '../../res/sample_pgrs_mapq_stat.csv'


def parse_bam(sample_id, cds_list, threads=1):
    samfile = ps.AlignmentFile(path_to_bam + sample_id + '_h37rv.bam', "rb", threads=threads)
    ref_name = samfile.get_reference_name(0)
    targer_genes = []
    nontarget_genes = []
    for cds in cds_list:
        av_mapq = 0
        read_num = 0
        for read in samfile.fetch(ref_name, cds.start, cds.end):
            av_mapq += read.mapping_quality
            read_num += 1
        if read_num > 0:
            if gene_filter in cds.name:
                targer_genes.append(av_mapq/read_num)
            else:
                nontarget_genes.append(av_mapq/read_num)
    samfile.close()
    return sample_id, ttest_ind(targer_genes, nontarget_genes)[1]


if __name__ == '__main__':
    sample_ids = [name.strip() for name in open(path_to_ids).readlines()]
    cds_list = read_annotations(upstream_length)
    # tasks = Parallel(n_jobs=-1)(delayed(parse_bam)(sample_id, cds_list) for sample_id in sample_ids)
    # with open(out_path, 'w') as f:
    #     for sample_id, p_val in tasks:
    #         f.write(sample_id + '\t' + str(p_val) + '\n')
    with open(out_path, 'w') as f:
        for sample_id in sample_ids:
            id, p_val = parse_bam(sample_id, cds_list, 144)
            f.write(sample_id + '\t' + str(p_val) + '\n')


