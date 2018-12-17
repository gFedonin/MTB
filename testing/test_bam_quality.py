from os.path import exists

from scipy.stats import ttest_ind
from sklearn.externals.joblib import Parallel, delayed
import numpy as np

import pysam as ps

from src.core.constants import data_path
from src.core.data_reading import read_variants

path_to_raw_var = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_mc10/'
path_to_ids = data_path + 'dr_covered_with_pheno_and_snp_new.txt'
path_to_bam = '/export/data/kchukreev/data/bam_files/gatk/realigned_bams_with_known/'
out_path = '../../res/bam_stat/'

proper_pairs_threshold = 0.5
long_isnerts_threshold = 0.5


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


def parse_bam(sample_id, var_list, threads=1):
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
    tasks = Parallel(n_jobs=-1)(delayed(parse_bam)(sample_id, var_list)
                                for sample_id, var_list in sample_to_variants.items())#, backend="threading"
    c = 0
    for task in tasks:
        c += 1
    return c


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


def print_bam_stat_by_variant(sample_to_variants, path):
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


if __name__ == '__main__':
    with open(path_to_ids, 'r') as f:
        sample_ids = [name.strip() for name in f.readlines()]
    tasks = Parallel(n_jobs=-1)(delayed(read_variants)(path_to_raw_var, sample_id) for sample_id in sample_ids)
    sample_to_variants = {}
    for name, snp_list in tasks:
        sample_to_variants[name] = snp_list
    print('read_var ok!')
    # parse_bam_files_parallel(sample_to_variants)
    # print_bam_stat(sample_to_variants, '../../res/gene_read_pairs_stat.csv')
    # print_bam_stat_by_sample(sample_to_variants, '../../res/snp_read_pairs_stat.csv')
    # print_bam_stat_multi_sample(sample_to_variants, '../../res/snp_read_pairs_stat_sample_share.csv')
