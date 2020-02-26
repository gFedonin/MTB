from os import listdir, makedirs
from os.path import exists

path_to_coeffs = '../../res/ml_gatk_vs_pilon_vs_old/ml_log_mc1_gatk_annotated_long_del_pg_NWds10_mq40_keep_complex_filtered/'
out_path = '../../res/ml_gatk_vs_pilon_vs_old/gatk/'

path1 = '../../res/ml_gatk_vs_pilon_vs_old/gatk/'
path2 = '../../res/ml_gatk_vs_pilon_vs_old/pilon/'
out_path1 = '../../res/ml_gatk_vs_pilon_vs_old/gatk_unique/'
out_path2 = '../../res/ml_gatk_vs_pilon_vs_old/pilon_unique/'
out_path3 = '../../res/ml_gatk_vs_pilon_vs_old/gatk_pilon/'


def filter_features():
    if not exists(out_path):
        makedirs(out_path)
    for fname in listdir(path_to_coeffs):
        if fname.endswith('.log'):
            filtered_list = []
            for l in open(path_to_coeffs + fname).readlines()[1:]:
                s = l.strip().split('\t')
                if float(s[-1]) > 0:
                    filtered_list.append(l)
            filtered_list.sort()
            with open(out_path + fname, 'w') as f:
                for l in filtered_list:
                    f.write(l)


def pick_unique():
    if not exists(out_path1):
        makedirs(out_path1)
    if not exists(out_path2):
        makedirs(out_path2)
    if not exists(out_path3):
        makedirs(out_path3)
    fname_to_features1 = {}
    for fname in listdir(path1):
        features = set()
        fname_to_features1[fname] = features
        for l in open(path1 + fname).readlines():
            s = l.strip().split('\t')
            features.add('\t'.join(s[:-1]))
    fname_to_features2 = {}
    for fname in listdir(path2):
        features = set()
        fname_to_features2[fname] = features
        for l in open(path2 + fname).readlines():
            s = l.strip().split('\t')
            features.add('\t'.join(s[:-1]))
    for fname, features1 in fname_to_features1.items():
        f1 = open(out_path1 + fname, 'w')
        f3 = open(out_path3 + fname, 'w')
        features2 = fname_to_features2[fname]
        for f in features1:
            if f in features2:
                f3.write(f + '\n')
            else:
                f1.write(f + '\n')
        f1.close()
        f3.close()
    for fname, features2 in fname_to_features2.items():
        with open(out_path2 + fname, 'w') as f2:
            features1 = fname_to_features1[fname]
            for f in features2:
                if f not in features1:
                    f2.write(f + '\n')


if __name__ == '__main__':
    pick_unique()
