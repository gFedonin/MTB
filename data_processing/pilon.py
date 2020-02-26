import gzip
from os import makedirs
from os.path import exists
from subprocess import check_call, call

from Bio import pairwise2
from sklearn.externals.joblib import Parallel, delayed

from core.constants import data_path
from core.data_reading import path_to_ref

# path_to_ids = data_path + 'coll18/coll1.txt'
from data_processing.parse_gatk_variants import in_filter_interval, get_intervals_to_filter_out

# path_to_ids = data_path + 'debug4.list'
# path_to_ids = data_path + 'list1.txt'
path_to_ids = data_path + 'all_with_pheno.txt'
path_to_bams = data_path + 'bwa_mem_rev/'
# path_to_bams = data_path + 'coll18/bwa_mem_rev/'
out_path = data_path + 'pilon/'
out_path_bam = data_path + 'pilon_minimap/'
out_path_variants = data_path + 'snps/pilon/raw_variants/'
out_path_reduced_vcf = data_path + 'snps/pilon/reduced_vcf/'

min_alt_frac = 0.75
min_dp = 10
min_mq = 30
min_bq = 20

thread_num = '8'


def run_pilon(sample_id):
    if not exists(out_path + sample_id):
        makedirs(out_path + sample_id)
    call('java -Xmx8G -jar ~/pilon-1.23.jar --threads ' + thread_num +
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
    call('pigz ' + out_path + sample_id + '/pilon.vcf', shell=True, cwd=out_path + sample_id)
    return 1


def pilon_all():
    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    tasks = Parallel(n_jobs=16)(delayed(run_pilon)(sample_id) for sample_id in sample_ids
                               if exists(path_to_bams + sample_id + '_h37rv.bam') and
                               not exists(out_path + sample_id + '/pilon.fasta'))
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
            stats[s[:i]] = s[i + 1:]
    return stats


def parse_pilon_vcf(sample_id, filter_intervals):
    filter_intervals_starts = [x[0] for x in filter_intervals]
    if not exists(out_path + sample_id + '/pilon.vcf.gz'):
        return 0
    if exists(out_path_variants + sample_id + '.variants'):
        with open(out_path_variants + sample_id + '.variants') as f:
            if f.readline() != '':
                return 0
    print(sample_id)
    with open(out_path_variants + sample_id + '.variants', 'w') as f:
        with open(out_path_reduced_vcf + sample_id + '.vcf', 'w') as fvcf:
            for l in gzip.open(out_path + sample_id + '/pilon.vcf.gz', 'rt').readlines():
                if l[0] == '#':
                    continue
                s = l.strip().split('\t')
                pos = s[1]
                ref = s[3]
                alt = s[4]
                if s[6] == 'PASS' or s[6] == 'Amb':
                    if alt == '.':
                        continue
                    stats = parse_stats(s[-3])
                    if 'IMPRECISE' in stats:
                        continue
                    stat_dp = stats.get('DP')
                    if stat_dp is not None and int(stat_dp) < min_dp:
                        continue
                    stat = stats.get('AF')
                    if stat is not None and float(stat) < min_alt_frac:
                        continue
                    stat = stats.get('MQ')
                    if stat is not None and int(stat) < min_mq:
                        continue
                    stat = stats.get('BQ')
                    if stat is not None and int(stat) < min_bq:
                        continue
                    stat = stats.get('XC')
                    if stat is not None and stat_dp is not None and (int(stat) > int(stat_dp)):
                        continue
                    if in_filter_interval(int(pos), filter_intervals_starts, filter_intervals):
                        continue
                    f.write(pos + '\t' + ref + '\t' + alt + '\n')
                    fvcf.write(l)
    return 1


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


def parse_complex_event(pos, old, new):
    res = []
    if new == '*':
        for i in range(len(old)):
            res.append((pos + i, '-', 'del'))
        return res
    if old == '*':
        res.append((pos, new, 'ins'))
        return res
    aln = pairwise2.align.globalms(old, new, 2, -1, -1, -.1)[0]
    aln_old = aln[0]
    aln_new = aln[1]
    insert_pos = 0
    insert_s = 0
    insert = False
    p = 0
    for i in range(len(aln_old)):
        c1 = aln_old[i]
        c2 = aln_new[i]
        if c1 != c2:
            if c1 == '-':
                # ins
                if not insert:
                    insert_pos = p
                    insert_s = i
                    insert = True
            elif c2 == '-':
                #del
                if insert:
                    res.append((pos + insert_pos, aln_new[insert_s: i], 'ins'))
                    insert = False
                res.append((pos + p, c2, 'del'))
                p += 1
            else:
                if insert:
                    res.append((pos + insert_pos, aln_new[insert_s: i], 'ins'))
                    insert = False
                res.append((pos + p, c2, 'snp'))
                p += 1
        else:
            if insert:
                res.append((pos + insert_pos, aln_new[insert_s: i], 'ins'))
                insert = False
            p += 1
    if insert:
        res.append((pos + insert_pos, aln_new[insert_s:], 'ins'))
    return res


def parse_all_vcf():
    if not exists(out_path_variants):
        makedirs(out_path_variants)
    if not exists(out_path_reduced_vcf):
        makedirs(out_path_reduced_vcf)
    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    filter_intervals = get_intervals_to_filter_out()
    tasks = Parallel(n_jobs=-1)(delayed(parse_pilon_vcf)(sample_id, filter_intervals) for sample_id in sample_ids)
    c = 0
    for task in tasks:
        c += task
    print('%d samples proccessed' % c)


if __name__ == '__main__':
    # pilon_all()
    # map_all_contigs()
    parse_all_vcf()
