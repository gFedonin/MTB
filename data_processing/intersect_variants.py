from os import makedirs
from os.path import exists

from core.constants import data_path
from core.data_reading import read_variants
from sklearn.externals.joblib import Parallel, delayed

path_to_ids = data_path + 'snps/gatk_pilon_old_intersection.list'
path1 = data_path + 'snps/pilon/annotated_pg_NWds10_filtered_test/'
path2 = data_path + 'snps/gatk_before_cortex/annotated_pg_NWds10_mq40_keep_complex_filtered_test/'
# out_path = data_path + 'snps/gatk_and_pilon/intersection/'
out_path = data_path + 'snps/gatk_and_pilon/unification/'

thread_num = 32


def read_all_variants(path_to_var, sample_ids):
    sample_to_variants = {}
    tasks = Parallel(n_jobs=thread_num, batch_size=len(sample_ids)//thread_num + 1)(delayed(read_variants)(path_to_var, sample_id)
                                    for sample_id in sample_ids if exists(path_to_var + sample_id + '.variants'))
    for name, snp_list in tasks:
        sample_to_variants[name] = snp_list
    return sample_to_variants


def intersect():
    if not exists(out_path):
        makedirs(out_path)
    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    sample_to_variants1 = read_all_variants(path1, sample_ids)
    sample_to_variants2 = read_all_variants(path2, sample_ids)
    for sample_id in sample_ids:
        vars1 = set(sample_to_variants1[sample_id])
        vars2 = sample_to_variants2[sample_id]
        with open(out_path + sample_id + '.variants', 'w') as f:
            for var in vars2:
                if var in vars1:
                    f.write(var + '\n')


def unite():
    if not exists(out_path):
        makedirs(out_path)
    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    sample_to_variants1 = read_all_variants(path1, sample_ids)
    sample_to_variants2 = read_all_variants(path2, sample_ids)
    for sample_id in sample_ids:
        vars = set(sample_to_variants1[sample_id])
        vars.update(sample_to_variants2[sample_id])
        with open(out_path + sample_id + '.variants', 'w') as f:
            for var in vars:
                f.write(var + '\n')


if __name__ == '__main__':
    # intersect()
    unite()