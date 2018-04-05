path_to_subsets = './data/subsets/'
path_to_pheno = './data/pheno/'


def read_subset(name):
    list = set()
    with open(path_to_subsets + name + '.txt', 'r') as f:
        for line in f.readlines():
            list.add(line[:-1])
    return name, list


def read_pheno(drug):
    sample_id_to_pheno = []
    with open(path_to_pheno + drug + '.pheno', 'r') as f:
        for line in f.readlines():
            i = line.find('\t')
            sample_id_to_pheno.append((line[:i], int(line[i + 1:-1])))
    return drug, sample_id_to_pheno