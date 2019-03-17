from os import listdir

from src.core.constants import data_path

out_path_win_cov = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_no_win_cov_m3/'

if __name__ == '__main__':
    for file_name in listdir(out_path_win_cov):
        if '.cov' in file_name:
            with open(out_path_win_cov + file_name) as f:
                lines = f.readlines()
            with open(out_path_win_cov + file_name, 'w') as f:
                l = lines[0]
                i = l.rfind('v')
                if i == -1:
                    print(file_name)
                    break
                else:
                    f.write(l[0: i + 1] + '\n')
                    f.write(l[i + 1:])
                    for line in lines[1:]:
                        f.write(line)