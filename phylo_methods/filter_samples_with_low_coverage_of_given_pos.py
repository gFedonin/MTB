from os import listdir

from src.core.constants import data_path

pos = '2938363'

path_to_cov = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_win_cov_m3_genes/'
out_path = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_win_cov_m3_genes/bad_samples.list'


def print_bad_samples():
    with open(out_path, 'w') as f:
        f.write('sample_id\tpos_is_undercovered\tpos_is_overcovered\twin_is_undercovered\t'
                'win_is_overcovered\twin_has_big_std_cov\n')
        for file in listdir(path_to_cov):
            if '.cov' in file:
                for line in open(path_to_cov + file).readlines()[1:]:
                    s = line.strip().split('\t')
                    if s[0] == pos:
                        if s[1] == '0':
                            f.write(file[:-4] + '\t' + '\t'.join(s[2:]) + '\n')
                        break


if __name__ == '__main__':
    print_bad_samples()
