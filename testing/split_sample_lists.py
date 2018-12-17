from src.core.constants import data_path

path_to_list = data_path + 'dr_covered_with_pheno_and_snp_new.txt'

split_num = 4

lines = open(path_to_list, 'r').readlines()
batch_size = len(lines)//split_num
for i in range(split_num):
    with open(data_path + 'list' + str(i) + '.list', 'w') as f:
        f.write(''.join(lines[i*batch_size: (i + 1)*batch_size]))
for i in range(len(lines) % split_num):
    with open(data_path + 'list' + str(i) + '.list', 'a') as f:
        f.write(lines[-i - 1])
