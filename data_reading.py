from constants import drug_names


def read_snp_list(path_to_snp_list, convert_to_int=False):
    list = []
    with open(path_to_snp_list, 'r') as f:
        for line in f.readlines():
            if convert_to_int:
                list.append(int(line.strip()))
            else:
                list.append(line.strip())
    return list


def read_subset(path_to_subsets, name):
    subset = set()
    with open(path_to_subsets + name + '.txt', 'r') as f:
        for line in f.readlines():
            subset.add(line.strip())
    return name, subset


def read_pheno(path_to_pheno, drug):
    sample_id_to_pheno = []
    with open(path_to_pheno + drug + '.pheno', 'r') as f:
        for line in f.readlines():
            line = line.strip()
            i = line.find('\t')
            sample_id_to_pheno.append((line[:i], int(line[i + 1:])))
    return drug, sample_id_to_pheno


def read_dict(path_to_dict, name):
    drug_to_mut_list = {}
    with open(path_to_dict + name + '.txt', 'r') as f:
        for line in f.readlines():
            line = line.strip()
            j = line.rfind('\t')
            drug = line[j + 1:]
            if drug in drug_names.keys():
                full_name = drug_names[drug]
                if type(full_name) is list:
                    for n in full_name:
                        if n in drug_to_mut_list.keys():
                            drug_to_mut_list[n].append(line[: j])
                        else:
                            drug_to_mut_list[n] = [line[: j]]
                else:
                    if full_name in drug_to_mut_list.keys():
                        drug_to_mut_list[full_name].append(line[: j])
                    else:
                        drug_to_mut_list[full_name] = [line[: j]]
            # else:
            #     print(drug + ' not found')
    return name, drug_to_mut_list


def read_variants(path_to_snp, sample_id, filter_set=None):
    snp_list = []
    with open(path_to_snp + sample_id + '.variants', 'r') as f:
        for line in f.readlines():
            l = line.strip()
            if l != '':
                if filter_set is not None:
                    if l in filter_set:
                        snp_list.append(l)
                else:
                    snp_list.append(l)
    return sample_id, snp_list