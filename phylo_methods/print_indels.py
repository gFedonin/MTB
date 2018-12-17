import os

from os.path import exists
from sklearn.externals.joblib import Parallel, delayed

from src.core.constants import data_path
from src.core.data_reading import read_variants

path_to_var = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_mc10/'
path_to_ids = data_path + 'dr_covered_with_pheno_and_snp.txt'
out_path = data_path + 'indels/'


def filter_indels(sample_id, variants):
    res = []
    del_len = 0
    del_start_pos = 0
    for var in variants:
        s = var.split('\t')
        if s[-1] == 'del':
            pos = int(s[0])
            if del_len == 0:
                del_start_pos = pos
                del_len = 1
            else:
                if del_start_pos + del_len == pos:
                    del_len += 1
                else:
                    res.append('del_' + str(del_start_pos) + '_' + str(del_len))
                    del_start_pos = pos
                    del_len = 1
        else:
            if del_len > 0:
                res.append('del_' + str(del_start_pos) + '_' + str(del_len))
                del_len = 0
            if s[-1] == 'ins':
                res.append('ins_' + s[0] + '_' + s[1])
    if del_len > 0:
        res.append('del_' + str(del_start_pos) + '_' + str(del_len))
    return sample_id, res


def main():
    if not exists(out_path):
        os.makedirs(out_path)
    with open(path_to_ids, 'r') as f:
        sample_ids = [name.strip() for name in f.readlines()]
    print(str(len(sample_ids)) + ' samples')

    indels_to_samples = {}
    tasks = Parallel(n_jobs=-1)(delayed(read_variants)(path_to_var, sample_id) for sample_id in sample_ids)
    tasks = Parallel(n_jobs=-1)(delayed(filter_indels)(sample_id, l) for sample_id, l in tasks)
    for name, l in tasks:
        for var in l:
            sample_set = indels_to_samples.get(var)
            if sample_set is None:
                sample_set = set()
                indels_to_samples[var] = sample_set
            sample_set.add(name)
    for var, sample_set in indels_to_samples.items():
        with open(out_path + var + '.samples', 'w') as f:
            for sample_id in sample_ids:
                f.write(sample_id)
                f.write('\t')
                if sample_id in sample_set:
                    f.write('1\n')
                else:
                    f.write('0\n')


if __name__ == '__main__':
    main()