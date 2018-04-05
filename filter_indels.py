from os.path import exists

from os import makedirs
from sklearn.externals.joblib import Parallel, delayed

path_to_ids = './data/all_with_pheno_and_snp.txt'
path_to_snps = './data/snps/realigned_vcfs/'
out_path = './data/snps/realigned_vcfs_indels/'


def filter_indels(sample_id):
    filtered = []
    with open(path_to_snps + sample_id + '_h37rv.vcf', 'r') as f:
        for line in f.readlines():
            if line[0] != '#':
                if 'TYPE=ins' in line or 'TYPE=del' in line:
                    filtered.append(line)
            else:
                filtered.append(line)

    with open(out_path + sample_id + '_h37rv.vcf', 'w') as f:
        f.write(''.join(filtered))
    return 0


def main():
    if not exists(out_path):
        makedirs(out_path)
    sample_ids = [sample_id.strip() for sample_id in open(path_to_ids, 'r').readlines()]
    for sample_id in sample_ids:
        # print(sample_id)
        filter_indels(sample_id)
    # tasks = Parallel(n_jobs=-1)(delayed(filter_indels)(sample_id) for sample_id in sample_ids)
    # c = 0
    # for task in tasks:
    #     c += task
    # print(c)


if __name__ == '__main__':
    main()
