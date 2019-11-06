from os import makedirs, listdir
from os.path import exists
from sklearn.externals.joblib import Parallel, delayed
from src.core.constants import data_path
from src.core.data_reading import read_all_variants

path_to_pheno = data_path + 'pheno_mc5_mega/'
path_to_var1 = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_no_win_qual_mqm_std3_mqm30_no_highcov/'
path_to_var2 = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_no_win_qual_mqm_std3_mqm30_no_highcov_bam_filtered/'
path_to_ids = data_path + 'snps/intersect.list'
path_to_subsets = data_path + 'subsets/'

out_path = '../../res/filter_comparison/'

# snp_count_threshold = 3

thread_num = 32

variant = '4355978\tCGA\tins'


def count_mutations(sample_to_vars):
    var_to_count = {}
    for sample_id, var_list in sample_to_vars.items():
        for var in var_list:
            if var in var_to_count:
                var_to_count[var] += 1
            else:
                var_to_count[var] = 1
    return var_to_count


def compare_on_intersection():
    if not exists(out_path):
        makedirs(out_path)
    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    sample_to_vars1 = read_all_variants(path_to_var1, sample_ids, thread_num=thread_num)
    var_to_count1 = count_mutations(sample_to_vars1)
    sample_to_vars2 = read_all_variants(path_to_var2, sample_ids, thread_num=thread_num)
    var_to_count2 = count_mutations(sample_to_vars2)
    with open(out_path + 'different.csv', 'w') as fout:
        for var, c1 in var_to_count1.items():
            c2 = var_to_count2.get(var)
            if c2 is None:
                fout.write('%s\t%d\t%d\n' % (var, c1, 0))
            else:
                if c1 != c2:
                    fout.write('%s\t%d\t%d\n' % (var, c1, c2))
        for var, c in var_to_count2.items():
            if var not in var_to_count1:
                fout.write('%s\t%d\t%d\n' % (var, 0, c))


def compare_on_all():
    if not exists(out_path):
        makedirs(out_path)
    sample_ids1 = [fName[:-9] for fName in listdir(path_to_var1) if fName.endswith('.variants')]
    sample_ids2 = [fName[:-9] for fName in listdir(path_to_var2) if fName.endswith('.variants')]
    sample_to_vars1 = read_all_variants(path_to_var1, sample_ids1, thread_num=thread_num)
    var_to_count1 = count_mutations(sample_to_vars1)
    sample_to_vars2 = read_all_variants(path_to_var2, sample_ids2, thread_num=thread_num)
    var_to_count2 = count_mutations(sample_to_vars2)
    with open(out_path + 'different.csv', 'w') as fout:
        for var, c1 in var_to_count1.items():
            c2 = var_to_count2.get(var)
            if c2 is None:
                fout.write('%s\t%d\t%d\n' % (var, c1, 0))
            else:
                if c1 != c2:
                    fout.write('%s\t%d\t%d\n' % (var, c1, c2))
        for var, c in var_to_count2.items():
            if var not in var_to_count1:
                fout.write('%s\t%d\t%d\n' % (var, 0, c))


def check_sample(sample_id, vars1, vars2):
    if variant in vars1:
        if variant not in vars2:
            return sample_id
    else:
        if variant in vars2:
            return sample_id
    return None


def print_list_of_different_samples():
    if not exists(out_path):
        makedirs(out_path)
    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    sample_to_vars1 = read_all_variants(path_to_var1, sample_ids, thread_num=thread_num)
    sample_to_vars2 = read_all_variants(path_to_var2, sample_ids, thread_num=thread_num)
    tasks = Parallel(n_jobs=thread_num, batch_size=len(sample_ids)//thread_num + 1)(
        delayed(check_sample)(sample_id, sample_to_vars1[sample_id], sample_to_vars2[sample_id]) for sample_id in sample_ids)
    with open(out_path + variant.replace('\t', '') + '.sample_list', 'w') as f:
        for sample_id in tasks:
            if sample_id is not None:
                f.write(sample_id + '\n')


if __name__ == '__main__':
    # compare_on_all()
    # compare_on_intersection()
	print_list_of_different_samples()