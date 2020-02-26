from core.constants import data_path

data_set = 'farhat19'
path_to_full_list = data_path + data_set + '/' + data_set + '.list'
path_to_filtered_list_samea = data_path + data_set + '/' + data_set + '_converted_new.samples'
metadata_path = data_path + data_set + '/metadata/'

path_to_pheno_raw = data_path + data_set + '/' + data_set + '.pheno'
path_to_pheno_filtered = data_path + data_set + '/' + data_set + '_converted.pheno'

batch_size = 100


def filter_samea():
    filter_set = set(l.strip() for l in open(path_to_filtered_list_samea).readlines())
    with open(path_to_pheno_filtered, 'w') as f:
        lines = open(path_to_pheno_raw).readlines()
        f.write(lines[0])
        for l in lines[1:]:
            s = l.strip().split()
            if s[0] in filter_set:
                f.write(l)


def filter_and_convert():
    filter_set = set(l.strip() for l in open(path_to_filtered_list_samea).readlines())
    meta_to_samea = []
    samples = [l.strip() for l in open(path_to_full_list).readlines()]
    for i in range(0, len(samples), batch_size):
        for l in open(metadata_path + str(i) + '.csv').readlines()[1:-1]:
            s = l.split(',')
            meta_to_samea.append((l, s[25]))
    with open(path_to_pheno_filtered, 'w') as f:
        lines = open(path_to_pheno_raw).readlines()
        f.write(lines[0])
        for l in lines[1:]:
            i = l.find('\t')
            name = l[:i]
            for m, samea in meta_to_samea:
                if name in m:
                    if samea in filter_set:
                        f.write(samea + l[i:])


if __name__ == '__main__':
    # filter_samea()
    filter_and_convert()
