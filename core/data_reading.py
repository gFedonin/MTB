from Bio import SeqIO

from src.core.constants import drug_names, data_path

path_to_ref = data_path + 'h37rv.fasta'


def read_h37rv():
    fasta_sequences = SeqIO.parse(open(path_to_ref), 'fasta')
    sequence = ''
    for fasta in fasta_sequences:
        sequence = str(fasta.seq)
    return sequence.upper()


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


def read_dict(path_to_dict: 'str', name: 'str')->'tuple[str, dict[str, list[str]]]':
    """
    Reads given dictionary and returns given dictionary name and map from drug name to mutation lists

    :param path_to_dict: str path to folder with dictionary files
    :param name: str dictionary name, no extension
    :return: tuple[str, dict[str,list[str]]] (name, drug_name to muy_list dictionary
    """
    drug_to_mut_list = {}
    with open(path_to_dict + name + '.txt', 'r') as f:
        for line in f.readlines():
            line = line.strip()
            j = line.rfind('\t')
            drug = line[j + 1:]
            mut = line[: j]
            if drug in drug_names.keys():
                full_name = drug_names[drug]
                if type(full_name) is list:
                    for n in full_name:
                        mut_list = drug_to_mut_list.get(n)
                        if mut_list is not None:
                            mut_list.append(mut)
                        else:
                            drug_to_mut_list[n] = [mut]
                else:
                    mut_list = drug_to_mut_list.get(full_name)
                    if mut_list is not None:
                        mut_list.append(mut)
                    else:
                        drug_to_mut_list[full_name] = [mut]
            # else:
            #     print(drug + ' not found')
    return name, drug_to_mut_list


def read_variants(path_to_snp, sample_id, filter_set=None, keep_type=True):
    snp_list = []
    with open(path_to_snp + sample_id + '.variants', 'r') as f:
        for line in f.readlines():
            l = line.strip()
            if not keep_type:
                i = l.rfind('\t')
                l = l[:i]
            if l != '':
                if filter_set is not None:
                    if l in filter_set:
                        snp_list.append(l)
                else:
                    snp_list.append(l)
    return sample_id, snp_list