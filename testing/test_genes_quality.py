from os.path import exists

from scipy.stats import ttest_ind
from sklearn.externals.joblib import Parallel, delayed
import numpy as np

import pysam as ps

from src.core.annotations import read_annotations, CDSType, localize_all_snps
from src.core.constants import data_path, upstream_length, ref_len
from src.core.data_reading import read_variants

path_to_var = data_path + 'snps/annotated_with_DR_with_indel_with_pheno_and_snp_mc10_long_del/'
path_to_raw_var = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_mc10/'
path_to_ids = data_path + 'dr_covered_with_pheno_and_snp_new.txt'#'list3.txt'
path_to_bam = '/export/data/kchukreev/data/bam_files/gatk/realigned_bams_with_known/'
out_path = '../../res/bam_stat/'

gene_group_label = 'PGRS'

proper_pairs_threshold = 0.5
long_isnerts_threshold = 0.5


def print_snp_density_stat(sample_to_variants, cds_list):
    gene_to_var_num = {}
    for name, snp_list in sample_to_variants.items():
        for var in snp_list:
            s = var.split('\t')
            if s[0] != 'non_cds' and s[0] != 'upstream':
                c = gene_to_var_num.get(s[1])
                if c is not None:
                    gene_to_var_num[s[1]] = c + 1
                else:
                    gene_to_var_num[s[1]] = 1
    labeled_counts = []
    unlabeled_counts = []
    for cds in cds_list:
        if cds.type != CDSType.upstream:
            if gene_group_label in cds.name:
                c = gene_to_var_num.get(cds.name)
                if c is None:
                    labeled_counts.append(0)
                else:
                    labeled_counts.append(c/(cds.end - cds.start))
            else:
                c = gene_to_var_num.get(cds.name)
                if c is None:
                    unlabeled_counts.append(0)
                else:
                    unlabeled_counts.append(c/(cds.end - cds.start))
    stat, p_value = ttest_ind(labeled_counts, unlabeled_counts)
    print("mean_var_density %s = %f non_%s = %f, TTest stat = %f, p_value = %f" % (gene_group_label,
            np.mean(labeled_counts), gene_group_label, np.mean(unlabeled_counts), stat, p_value))


def convert_var_list(sample_to_variants):
    variant_to_samples = {}
    for sample_id, var_list in sample_to_variants.items():
        for var in var_list:
            sample_list = variant_to_samples.get(var)
            if sample_list is None:
                sample_list = [sample_id]
                variant_to_samples[var] = sample_list
            else:
                sample_list.append(sample_id)
    return variant_to_samples


class Read:
    name: str
    start: int
    end: int
    is_reverse: bool
    is_first: bool

    def __init__(self, read):
        self.name = read.query_name
        self.start = read.reference_start
        self.end = read.reference_end
        self.is_first = read.is_read1


class VarStat:

    def __init__(self, cov=0, proper_pairs=0, long_inserts=0):
        self.proper_pairs = proper_pairs
        self.cov = cov
        self.insert_lengths = []
        self.long_inserts = long_inserts


def parse_bam_pos_wise_and_all(sample_id, var_list, threads=1):
    var_to_stat = {}
    if exists(out_path + sample_id + '.stat'):
        if len(open(out_path + sample_id + '.stat', 'r').readlines()) - 1 == len(var_list):
            return
    samfile = ps.AlignmentFile(path_to_bam + sample_id + '_h37rv.bam', "rb", threads=threads)
    ref_name = samfile.get_reference_name(0)
    var_pos = []
    read_pairs = {}
    read_to_vars = {}
    total_proper_pairs = 0
    for var in var_list:
        s = var.split('\t')
        pos = int(s[0]) - 1
        var_pos.append(pos)
        var_stat = VarStat()
        var_to_stat[var] = var_stat
        for read in samfile.fetch(ref_name, pos, pos + 1):
            var_stat.cov += 1
            if read.is_proper_pair:
                var_stat.proper_pairs += 1
                total_proper_pairs += 1
                if read.query_name not in read_pairs:
                    read_pairs[read.query_name] = Read(read)
                    read_to_vars[read.query_name] = [var]
                else:
                    read_to_vars[read.query_name].append(var)
    mean = 0
    std = 0
    if total_proper_pairs > 0:
        for read in samfile.fetch():
            p_read = read_pairs.get(read.query_name)
            if p_read is not None:
                if p_read.is_first == read.is_read1:
                    continue
                for var in read_to_vars[read.query_name]:
                    var_stat = var_to_stat[var]
                    if read.is_reverse:
                        var_stat.insert_lengths.append(read.reference_end - p_read.start)
                    else:
                        var_stat.insert_lengths.append(p_read.end - read.reference_start)
        c = 0
        for var, var_stat in var_to_stat.items():
            mean += sum(var_stat.insert_lengths)
            std += sum(x*x for x in var_stat.insert_lengths)
            c += len(var_stat.insert_lengths)
        mean /= c
        std = np.sqrt(std/c - mean*mean)
    samfile.close()
    with open(out_path + sample_id + '.stat', 'w') as f:
        f.write('var\tcov\tproper_pairs\tlong_inserts\n')
        for var, var_stat in var_to_stat.items():
            long_inserts = sum(abs(x - mean) > 2 * std for x in var_stat.insert_lengths)
            f.write('%s\t%d\t%d\t%d\n' % (var, var_stat.cov, var_stat.proper_pairs, long_inserts))
    # print(file + ' processing done')


def parse_bam_files_parallel(sample_to_variants):
    # sample_to_varStat = {}
    tasks = Parallel(n_jobs=-1)(delayed(parse_bam_pos_wise_and_all)(sample_id, var_list)
                                for sample_id, var_list in sample_to_variants.items())#, backend="threading"
    c = 0
    for task in tasks:
        c += 1
    return c
    # for sample_id, var_to_stat in tasks:
    #     sample_to_varStat[sample_id] = var_to_stat
    # return sample_to_varStat


def get_bam_stat(sample_to_variants):
    variant_to_bam_stat = {}
    # sample_to_var_stat = read_bam_files_parallel(sample_to_variants)
    sample_to_var_stat = {}
    for sample_id, var_list in sample_to_variants.items():
        var_to_stat = {}
        sample_to_var_stat[sample_id] = var_to_stat
        with open(out_path + sample_id + '.stat', 'r') as f:
            for line in f.readlines()[1:]:
                var_stat = VarStat()
                s = line.strip().split('\t')
                var_stat.cov = int(s[-3])
                var_stat.proper_pairs = int(s[-2])
                var_stat.long_inserts = int(s[-1])
                var_to_stat['\t'.join(s[:-3])] = var_stat
        if len(var_list) != len(var_to_stat.keys()):
            print(sample_id + ' var_num = ' + str(len(var_list)) + ' found ' + str(len(var_to_stat.keys())))

    variant_to_samples = convert_var_list(sample_to_variants)
    for var, sample_list in variant_to_samples.items():
        cov = 0
        proper_pairs = 0
        long_inserts = 0
        for sample in sample_list:
            var_stat = sample_to_var_stat[sample][var]
            cov += var_stat.cov
            proper_pairs += var_stat.proper_pairs
            long_inserts += var_stat.long_inserts
        variant_to_bam_stat[var] = VarStat(cov, proper_pairs, long_inserts)
    return variant_to_bam_stat


def print_bam_stat(sample_to_variants, path):
    variant_to_bam_stat = get_bam_stat(sample_to_variants)
    with open(path, 'w') as f:
        f.write('var\tcov\tproper\tlong_inserts\n')
        for var, stat in variant_to_bam_stat.items():
            f.write('%s\t%d\t%d\t%d\n' % (var, stat.cov, stat.proper_pairs, stat.long_inserts))


def print_bam_stat_by_sample(sample_to_variants, path):
    with open(path, 'w') as fout:
        fout.write('sample_id\tpos\talt\ttype\tcov\tproper\tlong_inserts\tbad_pairs_freq\tlong_pairs_freq\n')
        for sample_id, var_list in sample_to_variants.items():
            with open(out_path + sample_id + '.stat', 'r') as fin:
                for line in fin.readlines()[1:]:
                    line = line.strip()
                    s = line.split('\t')
                    cov = int(s[-3])
                    proper_pairs = int(s[-2])
                    long_inserts = int(s[-1])
                    bad_pairs_freq = (cov - proper_pairs)/cov
                    if proper_pairs > 0:
                        long_pairs_freq = long_inserts/proper_pairs
                    else:
                        long_pairs_freq = 0
                    if bad_pairs_freq > proper_pairs_threshold or long_pairs_freq > long_isnerts_threshold:
                        fout.write('%s\t%s\t%1.2f\t%1.2f\n' % (sample_id, line, bad_pairs_freq, long_pairs_freq))


def print_bam_stat_multi_sample(sample_to_variants, path):
    sample_to_var_stat = {}
    for sample_id, var_list in sample_to_variants.items():
        var_to_stat = {}
        sample_to_var_stat[sample_id] = var_to_stat
        with open(out_path + sample_id + '.stat', 'r') as f:
            for line in f.readlines()[1:]:
                var_stat = VarStat()
                s = line.strip().split('\t')
                var_stat.cov = int(s[-3])
                var_stat.proper_pairs = int(s[-2])
                var_stat.long_inserts = int(s[-1])
                var_to_stat['\t'.join(s[:-3])] = var_stat
        if len(var_list) != len(var_to_stat.keys()):
            print(sample_id + ' var_num = ' + str(len(var_list)) + ' found ' + str(len(var_to_stat.keys())))
    with open(path, 'w') as fout:
        variant_to_samples = convert_var_list(sample_to_variants)
        fout.write('pos\talt\ttype\tbad_pairs_freq_share\tlong_pairs_freq_share\n')
        for var, sample_list in variant_to_samples.items():
            bad_pairs_freq_share = 0
            long_pairs_freq_share = 0
            for sample in sample_list:
                var_stat = sample_to_var_stat[sample][var]
                bad_pairs_freq = (var_stat.cov - var_stat.proper_pairs) / var_stat.cov
                if var_stat.proper_pairs > 0:
                    long_pairs_freq = var_stat.long_inserts / var_stat.proper_pairs
                else:
                    long_pairs_freq = 0
                if bad_pairs_freq > proper_pairs_threshold:
                    bad_pairs_freq_share += 1
                if long_pairs_freq > long_isnerts_threshold:
                    long_pairs_freq_share += 1
            fout.write('%s\t%1.2f\t%1.2f\n' % (var, bad_pairs_freq_share/len(sample_list), long_pairs_freq_share/len(sample_list)))


class GeneStat:

    def __init__(self):
        self.mut_num: int = 0
        self.cov: int = 0
        self.proper_pairs: int = 0
        self.long_inserts: int = 0


def gene_bam_stat(sample_to_variants, path):
    variant_to_bam_stat = get_bam_stat(sample_to_variants)
    cds_list = read_annotations(upstream_length)
    pos_to_stat = {}
    for var, stat in variant_to_bam_stat.items():
        pos = int(var.split('\t')[0])
        if pos not in pos_to_stat:
            pos_to_stat[pos] = [stat]
        else:
            pos_to_stat[pos].append(stat)
    print('total pos %d' % len(pos_to_stat))
    pos_to_cds = localize_all_snps(list(pos_to_stat.keys()), cds_list)
    print('localized %d pos' % len(pos_to_cds))
    gene_to_stat = {cds.name: GeneStat() for cds in cds_list}
    print('total genes %d' % len(gene_to_stat))
    for pos, stat_list in pos_to_stat.items():
        cds = pos_to_cds.get(pos)
        if cds is not None:
            g_stat = gene_to_stat[cds.name]
            for stat in stat_list:
                if stat.cov == 0:
                    print('AAA!')
                g_stat.mut_num += 1
                g_stat.cov += stat.cov
                g_stat.proper_pairs += stat.proper_pairs
                g_stat.long_inserts += stat.long_inserts
    labeled_counts = []
    unlabeled_counts = []
    with open(path, 'w') as fout:
        fout.write('gene name\tmut_num\timproper_pairs_freq\tlong_inserts_freq\n')
        for gene_name, g_stat in gene_to_stat.items():
            if g_stat.cov != 0:
                improper_pairs_freq = (g_stat.cov - g_stat.proper_pairs)/g_stat.cov
                long_inserts_freq = g_stat.long_inserts/g_stat.proper_pairs
                fout.write('%s\t%d\t%1.2f\t%1.2f\n' % (gene_name, g_stat.mut_num, improper_pairs_freq, long_inserts_freq))
                if gene_group_label in gene_name:
                    labeled_counts.append(improper_pairs_freq)
                else:
                    unlabeled_counts.append(improper_pairs_freq)
    stat, p_value = ttest_ind(labeled_counts, unlabeled_counts)
    print("improper pairs freq %s = %f non_%s = %f, TTest stat = %f, p_value = %f" % (gene_group_label,
            np.mean(labeled_counts), gene_group_label, np.mean(unlabeled_counts), stat, p_value))
    labeled_counts = []
    unlabeled_counts = []
    for gene_name, g_stat in gene_to_stat.items():
        if g_stat.cov > 0:
            if gene_group_label in gene_name:
                if g_stat.proper_pairs > 0:
                    labeled_counts.append(g_stat.long_inserts/g_stat.proper_pairs)
            else:
                if g_stat.proper_pairs > 0:
                    unlabeled_counts.append(g_stat.long_inserts/g_stat.proper_pairs)
    stat, p_value = ttest_ind(labeled_counts, unlabeled_counts)
    print("long insertion lengths freq %s = %f non_%s = %f, TTest stat = %f, p_value = %f" % (gene_group_label,
            np.mean(labeled_counts), gene_group_label, np.mean(unlabeled_counts), stat, p_value))


if __name__ == '__main__':
    with open(path_to_ids, 'r') as f:
        sample_ids = [name.strip() for name in f.readlines()]
    tasks = Parallel(n_jobs=-1)(delayed(read_variants)(path_to_raw_var, sample_id) for sample_id in sample_ids)
    sample_to_variants = {}
    for name, snp_list in tasks:
        sample_to_variants[name] = snp_list
    print('read_var ok!')
    # cds_list = read_annotations(upstream_length)
    # print('read_ann ok!')
    # print_snp_density_stat(sample_to_variants, cds_list)
    # parse_bam_files_parallel(sample_to_variants)
    # print_bam_stat(sample_to_variants, '../../res/gene_read_pairs_stat.csv')
    # print_bam_stat_by_sample(sample_to_variants, '../../res/snp_read_pairs_stat.csv')
    # print_bam_stat_multi_sample(sample_to_variants, '../../res/snp_read_pairs_stat_sample_share.csv')
    gene_bam_stat(sample_to_variants, '../../res/gene_read_pairs_stat.csv')
