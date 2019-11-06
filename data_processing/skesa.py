import gzip
from os import system, listdir, makedirs
from os.path import exists

from sklearn.externals.joblib import Parallel, delayed

from core.data_reading import path_to_ref
from core.constants import data_path

path_to_skesa = '/export/home/fedonin/SKESA/skesa.centos6.10'
# path_to_mummer = '/export/home/fedonin/MUMmer3.23/'
path_to_mummer = '/export/home/fedonin/mummer4/bin/'
# path_to_samples = '/export/data/kkuleshov/myc/sra/'
# path_to_samples = data_path + 'coll18/coll18_bbduk_q30/'
path_to_samples = data_path + 'kkuleshov_bbduk_q30/'
# path_to_samples = data_path + 'walker18/walker18_bbduk_q30/'
path_to_samples_reformat = data_path + 'kkuleshov_bbduk_q30_fixed/'
# path_to_samples_reformat = data_path + 'walker18/walker18_bbduk_q30_fixed/'
# path_to_samples_reformat = data_path + 'coll18/coll18_bbduk_q30_fixed/'
# path_to_list = data_path + 'coll18/coll5.txt'
# path_to_list = data_path + 'walker18/walker2018.list'
# path_to_list = data_path + 'coll18/coll18_supp.samples'
path_to_list = data_path + 'all_with_pheno.txt'
# path_to_list = data_path + 'debug1.list'
skesa_out_path = data_path + 'skesa_bbduk_q30/'
# skesa_out_path = data_path + 'coll18/skesa_bbduk_q30/'
# skesa_out_path = data_path + 'walker18/skesa_bbduk_q30/'
# skesa_out_path = data_path + 'test/'
# skesa_out_path = data_path + 'skesa_coll_bbduk/'
skesa_log_path = data_path + 'skesa_bbduk_q30_log/'
# skesa_log_path = data_path + 'walker18/skesa_bbduk_q30_log/'
# skesa_log_path = data_path + 'coll18/skesa_bbduk_q30_log/'
# skesa_log_path = data_path + 'test_log/'
# mummer_out_path = data_path + 'mummer/delta/'
# mummer_out_path = data_path + 'walker18/mummer/dnadiff_bbduk_q30/'
# mummer_out_path = data_path + 'coll18/mummer/dnadiff_bbduk_q30/'
mummer_out_path = data_path + 'mummer/dnadiff_bbduk_q30/'

# mummer_log_path = data_path + 'mummer/nucmer_log/'
# mummer_log_path = data_path + 'walker18/mummer/log_bbduk_q30/'
# mummer_log_path = data_path + 'coll18/mummer/log_bbduk_q30/'
mummer_log_path = data_path + 'mummer/log_bbduk_q30/'

# mummer_out_path_filtered = data_path + 'mummer/delta_filtered/'
# mummer_out_path_snps = data_path + 'walker18/mummer/snps_bbduk_q30/'
# mummer_out_path_snps = data_path + 'coll18/mummer/snps_bbduk_q30/'
mummer_out_path_snps = data_path + 'mummer/snps_bbduk_q30/'

# final_raw_snps_path = data_path + 'walker18/snps/skesa_mummer_raw_ld_mum4_q30/'
# final_raw_snps_path = data_path + 'coll18/snps/skesa_mummer_raw_ld_mum4_q30/'
final_raw_snps_path = data_path + 'snps/skesa_mummer_raw_ld_mum4_q30/'



thread_num = '144'


def run(sample_id):
    if not exists(skesa_out_path + sample_id + '.skesa.fa'):
        system(path_to_skesa + ' --fastq ' + path_to_samples_reformat + sample_id + '/' + sample_id + '_p1.fastq.gz,' +
               path_to_samples_reformat + sample_id + '/' + sample_id + '_p2.fastq.gz --cores ' + thread_num +
               ' --memory 100 > ' + skesa_out_path + sample_id +
               '.skesa.fa 2> ' + skesa_log_path + sample_id + '.skesa.log')
        return 1
    return 0


def run_skesa():
    if not exists(skesa_out_path):
        makedirs(skesa_out_path)
    if not exists(skesa_log_path):
        makedirs(skesa_log_path)
    tasks = Parallel(n_jobs=5)(delayed(run)(l.strip()) for l in open(path_to_list).readlines())
    c = 0
    for task in tasks:
        c += task
    print('%d samples processed' % c)


# def run_mummer(sample_id):
#     if not exists(mummer_out_path_filtered + sample_id + '.snps'):
#         system(path_to_mummer + 'nucmer --prefix=' + mummer_out_path + sample_id + ' ' + path_to_ref + ' ' +
#                skesa_out_path + sample_id + '.skesa.fa > ' + mummer_log_path + sample_id + '.nucmer.log 2>> ' +
#                mummer_log_path + sample_id + '.nucmer.log')
#         system(path_to_mummer + 'delta-filter -r -q ' + mummer_out_path + sample_id + '.delta > ' +
#                mummer_out_path_filtered + sample_id + '.filter')
#         system(path_to_mummer + 'show-snps -Clr ' + mummer_out_path_filtered + sample_id + '.filter > ' +
#                mummer_out_path_snps + sample_id + '.snps')
#         return 1
#     return 0


def run_mummer(sample_id):
    if not exists(mummer_out_path + sample_id):
        makedirs(mummer_out_path + sample_id)
    if not exists(mummer_out_path + sample_id + '/' + sample_id + '.snps'):
        contigs = [l for l in open(skesa_out_path + sample_id + '.skesa.fa').readlines()]
        if len(contigs) == 0:
            print('no contigs for ' + sample_id)
            return 1
        status = system(path_to_mummer + 'dnadiff  --prefix=' + mummer_out_path + sample_id + '/' + sample_id + ' ' + path_to_ref +
               ' ' + skesa_out_path + sample_id + '.skesa.fa > ' + mummer_log_path + sample_id + '.dnadiff.log 2>> ' +
               mummer_log_path + sample_id + '.dnadiff.log')
        if status != 0:
            print('mummer problem with ' + sample_id)
        return 1
    return 0


def mass_mummer():
    if not exists(mummer_out_path):
        makedirs(mummer_out_path)
    if not exists(mummer_log_path):
        makedirs(mummer_log_path)
    # if not exists(mummer_out_path_filtered):
    #     makedirs(mummer_out_path_filtered)
    if not exists(mummer_out_path_snps):
        makedirs(mummer_out_path_snps)
    tasks = Parallel(n_jobs=-1)(delayed(run_mummer)(l[:-9]) for l in listdir(skesa_out_path))
    c = 0
    for task in tasks:
        c += task
    print('%d samples processed' % c)


# def parse_mummer_snps(sample_id):
#     with open(final_raw_snps_path + sample_id + '.variants', 'w') as f:
#         is_ins = False
#         ins_str = []
#         ins_pos = None
#         is_del = False
#         del_len = 0
#         del_pos = None
#         for l in open(mummer_out_path_snps + sample_id + '.snps').readlines()[5:]:
#             s = l.strip().split()
#             if s[1] == '.':
#                 if is_del:
#                     f.write('%s\t%d\tdel\n' % (del_pos, del_len))
#                     is_del = False
#                     del_len = 0
#                     del_pos = None
#                 is_ins = True
#                 ins_pos = s[0]
#                 ins_str.append(s[2])
#             elif s[2] == '.':
#                 if is_ins:
#                     f.write('%s\t%s\tins\n' % (ins_pos, ''.join(ins_str)))
#                     is_ins = False
#                     ins_str.clear()
#                     ins_pos = None
#                 is_del = True
#                 del_pos = s[0]
#                 del_len += 1
#             else:
#                 if is_del:
#                     f.write('%s\t%d\tdel\n' % (del_pos, del_len))
#                     is_del = False
#                     del_len = 0
#                     del_pos = None
#                 elif is_ins:
#                     f.write('%s\t%s\tins\n' % (ins_pos, ''.join(ins_str)))
#                     is_ins = False
#                     ins_str.clear()
#                     ins_pos = None
#                 f.write('%s\t%s\tsnp\n' % (s[0], s[2]))
#     return 1


def parse_mummer_snps(sample_id):
    if not exists(mummer_out_path + sample_id + '/' + sample_id + '.snps'):
        print('no snps for ' + sample_id)
        return 1
    with open(final_raw_snps_path + sample_id + '.variants', 'w') as f:
        is_ins = False
        ins_str = []
        ins_pos = None
        is_del = False
        del_len = 0
        del_pos = None
        for l in open(mummer_out_path + sample_id + '/' + sample_id + '.snps').readlines()[5:]:
            s = l.strip().split()
            if s[1] == '.':
                if is_del:
                    f.write('%s\t%d\tdel\n' % (del_pos, del_len))
                    is_del = False
                    del_len = 0
                    del_pos = None
                is_ins = True
                ins_pos = s[0]
                ins_str.append(s[2])
            elif s[2] == '.':
                if is_ins:
                    f.write('%s\t%s\tins\n' % (ins_pos, ''.join(ins_str)))
                    is_ins = False
                    ins_str.clear()
                    ins_pos = None
                is_del = True
                del_pos = s[0]
                del_len += 1
            else:
                if is_del:
                    f.write('%s\t%d\tdel\n' % (del_pos, del_len))
                    is_del = False
                    del_len = 0
                    del_pos = None
                elif is_ins:
                    f.write('%s\t%s\tins\n' % (ins_pos, ''.join(ins_str)))
                    is_ins = False
                    ins_str.clear()
                    ins_pos = None
                f.write('%s\t%s\tsnp\n' % (s[0], s[2]))
    return 1


# def parse_all_mummer_snps():
#     if not exists(final_raw_snps_path):
#         makedirs(final_raw_snps_path)
#     tasks = Parallel(n_jobs=-1)(delayed(parse_mummer_snps)(l[:-5]) for l in listdir(mummer_out_path_snps))
#     c = 0
#     for task in tasks:
#         c += task
#     print('%d samples processed' % c)


def parse_all_mummer_snps():
    if not exists(final_raw_snps_path):
        makedirs(final_raw_snps_path)
    tasks = Parallel(n_jobs=-1)(delayed(parse_mummer_snps)(l) for l in listdir(mummer_out_path))
    c = 0
    for task in tasks:
        c += task
    print('%d samples processed' % c)


def convert_names(sample_id):
    if not exists(path_to_samples_reformat + sample_id):
        makedirs(path_to_samples_reformat + sample_id)
    if not exists(path_to_samples + sample_id + '/' + sample_id + '_p1.fastq.gz') or not exists(path_to_samples + sample_id + 
        '/' + sample_id + '_p2.fastq.gz'):
        return 0
    if exists(path_to_samples_reformat + sample_id + '/' + sample_id + '_p1.fastq.gz'):
        return 0
    lines = [l for l in gzip.open(path_to_samples + sample_id + '/' + sample_id + '_p1.fastq.gz', 'rt').readlines()]
    read_num = len(lines)//4
    with gzip.open(path_to_samples_reformat + sample_id + '/' + sample_id + '_p1.fastq.gz', 'wt', 2) as f:
        for i in range(read_num):
            name = lines[4*i][1:]
            f.write(lines[4*i])
            f.write(lines[4*i + 1])
            f.write('+' + name)
            f.write(lines[4*i + 3])
    if exists(path_to_samples_reformat + sample_id + '/' + sample_id + '_p2.fastq.gz'):
        return 0
    lines = [l for l in gzip.open(path_to_samples + sample_id + '/' + sample_id + '_p2.fastq.gz', 'rt').readlines()]
    read_num = len(lines)//4
    with gzip.open(path_to_samples_reformat + sample_id + '/' + sample_id + '_p2.fastq.gz', 'wt', 2) as f:
        for i in range(read_num):
            name = lines[4*i][1:]
            f.write(lines[4*i])
            f.write(lines[4*i + 1])
            f.write('+' + name)
            f.write(lines[4*i + 3])
    return 1


def convert_all_names():
    tasks = Parallel(n_jobs=-1)(delayed(convert_names)(l.strip()) for l in open(path_to_list).readlines())
    c = 0
    for task in tasks:
        c += task
    print('%d samples processed' % c)


if __name__ == '__main__':
    # convert_all_names()
    # run_skesa()
    mass_mummer()
    parse_all_mummer_snps()

