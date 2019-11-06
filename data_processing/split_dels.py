from os import listdir, makedirs
from os.path import exists

from core.constants import data_path

path_to_variants = data_path + 'snps/skesa_mummer_raw_ld_mum4/'
out_path = data_path + 'snps/skesa_mummer_raw_mum4/'


if __name__ == '__main__':
    if not exists(out_path):
        makedirs(out_path)
    for fname in listdir(path_to_variants):
        if fname.endswith('.variants'):
            with open(out_path + fname, 'w') as f:
                for l in open(path_to_variants + fname).readlines():
                    s = l.strip().split('\t')
                    if s[-1] == 'del':
                        del_len = int(s[1])
                        if del_len == 1:
                            f.write(l)
                        else:
                            pos = int(s[0])
                            for i in range(del_len):
                                f.write(str(pos + i) + '\t1\tdel\n')
                    else:
                        f.write(l)