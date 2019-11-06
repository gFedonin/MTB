import gzip
from os import makedirs
from os.path import exists
from subprocess import check_call, call

from sklearn.externals.joblib import Parallel, delayed

from core.constants import data_path
from core.data_reading import path_to_ref

# path_to_ids = data_path + 'coll18/coll1.txt'
# path_to_ids = data_path + 'debug1.list'
path_to_ids = data_path + 'list_54.txt'
# path_to_ids = data_path + 'all_with_pheno.txt'
path_to_bams = data_path + 'bwa_mem_rev/'
# path_to_bams = data_path + 'coll18/bwa_mem_rev/'
out_path = data_path + 'pilon/'
out_path_bam = data_path + 'pilon_minimap/'
out_path_variants = data_path + 'snps/pilon/raw_variants/'

min_alt_frac = 0.75
min_dp = 10
min_mq = 30
min_bq = 20

thread_num = '8'


def run_pilon(sample_id):
    if not exists(out_path + sample_id):
        makedirs(out_path + sample_id)
    check_call('java -Xmx4G -jar ~/pilon-1.23.jar --threads ' + thread_num +
               ' --variant --vcf --genome ' + path_to_ref +
               ' --frags ' + path_to_bams + sample_id + '_h37rv.bam --outdir ' + out_path +
               sample_id + ' > ' + out_path + sample_id + '/log.txt 2> ' + out_path + sample_id + '/log.txt',
               shell=True, cwd=out_path + sample_id)
    # with open(out_path + sample_id + '/' + sample_id + '.vcf', 'w') as f:
    #     for l in open(out_path + sample_id + '/pilon.vcf').readlines():
    #     # for l in gzip.open(out_path + sample_id + '/pilon.vcf.gz', 'rt').readlines():
    #         if l[0] == '#':
    #             f.write(l)
    #             continue
    #         s = l.strip().split('\t')
    #         if s[4] == '.' and s[6] == 'PASS':
    #             continue
    #         f.write(l)
    # check_call('rm ' + out_path + sample_id + '/pilon.vcf', shell=True)
    call('gzip ' + out_path + sample_id + '/pilon.vcf', shell=True, cwd=out_path + sample_id)
    return 1


def pilon_all():
    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    tasks = Parallel(n_jobs=8)(delayed(run_pilon)(sample_id) for sample_id in sample_ids
                               if exists(path_to_bams + sample_id + '_h37rv.bam') and
                               not exists(out_path + sample_id + 'pilon.fasta'))
    c = 0
    for task in tasks:
        c += task
    print('%d samples processed' % c)


def parse_stats(stat_str):
    stats = {}
    for s in stat_str.split(';'):
        i = s.find('=')
        if i == -1:
            stats[s] = None
        else:
            stats[s[:i]] = s[i + 1]
    return stats


def parse_pilon_vcf():
    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    for sample_id in sample_ids:
        with open(out_path_variants + sample_id + '.variants', 'w') as f:
            for l in gzip.open(out_path + sample_id + '/pilon.vcf.gz', 'rt').readlines():
                if l[0] == '#':
                    continue
                s = l.strip().split('\t')
                pos = int(s[1])
                ref = s[3]
                alt = s[4]
                if s[6] == 'PASS' or s[6] == 'Amb':
                    stats = parse_stats(s[-3])
                    if 'IMPRECISE' in stats:
                        continue
                    if int(stats['DP']) < min_dp:
                        continue
                    if float(stats['AF']) < min_alt_frac:
                        continue
                    if int(stats['MQ']) < min_mq:
                        continue
                    if int(stats['BQ']) < min_bq:
                        continue
                    if int(stats['XC']) > int(stats['DP']):
                        continue


def minimap(id):
    if not exists(out_path + id + '.bam'):
        check_call('/export/home/fedonin/minimap2-2.17_x64-linux/minimap2 -x asm5 -a -t ' + thread_num +
                    ' ' + data_path + 'h37rv.fasta ' + out_path + id +
                   '/pilon.fasta | samtools sort -@ ' + thread_num + ' - | samtools view -bS -@ ' + thread_num + ' - > ' +
                   out_path_bam + id + '.bam', shell=True)
        check_call('samtools index -@ ' + thread_num + ' ' + out_path_bam + id + '.bam', shell=True)


def map_all_contigs():
    if not exists(out_path):
        makedirs(out_path)
    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    for sample_id in sample_ids:
        minimap(sample_id)


if __name__ == '__main__':
    pilon_all()
    # map_all_contigs()
