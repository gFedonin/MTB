import os
from os.path import exists

from src.core.annotations import read_annotations
from src.core.constants import data_path, upstream_length
from src.core.data_reading import read_dict

path_to_features = '../../res/ml_log_mc1-noncds_long_del/'
path_to_dict = data_path + 'dictionaries_indels/'
dict_name = 'Walker_dictionary_broken_genes'

out_path = '../../res/selected_features_stat/ml_log_mc1-noncds_long_del/'


def filter_selected():
    for (dirpath, dirnames, filenames) in os.walk(path_to_features):
        for filename in filenames:
            i = filename.find('.log')
            if i != -1:
                drug = filename[0: i]
                with open(path_to_features + drug + '.selected', 'w') as f:
                    for line in open(path_to_features + filename, 'r').readlines()[1:]:
                        s = line.strip().split('\t')
                        if float(s[-1]) != 0:
                            f.write(line)


def find_missing_from_dict():
    if not exists(out_path):
        os.makedirs(out_path)
    name, drug_to_mut_list = read_dict(path_to_dict, dict_name)
    for drug, mut_list in drug_to_mut_list.items():
        with open(out_path + drug + '.stat', 'w') as f:
            broken = set()
            changed = set()
            all = set()
            total = 0
            for line in open(path_to_features + drug+ '.selected').readlines():#
                s = line.strip().split('\t')
                total += 1
                all.add('\t'.join(s[:-1]))
                if s[-2] == 'broken':
                    broken.add(s[1])
                elif s[-2] == 'changed':
                    changed.add(s[1])

            in_broken = 0
            in_changed = 0
            point = 0
            missing = 0
            missing_set = []
            for m in mut_list:
                if m in all:
                    point += 1
                else:
                    s = m.split('\t')
                    if s[1] in broken:
                        in_broken += 1
                    elif s[1] in changed:
                        in_changed += 1
                    else:
                        missing += 1
                        missing_set.append(m)
            f.write('total selected %d\n' % (total))
            f.write('broken %d\n' % (len(broken)))
            f.write('changed %d\n' % (len(changed)))
            f.write('Walker total %d\n' % (len(mut_list)))
            f.write('Walker point hit %d\n' % (point))
            f.write('Walker changed hit %d\n' % (in_changed))
            f.write('Walker broken hit %d\n' % (in_broken))
            f.write('Walker missing %d\n' % (missing))
            f.write('missing list:\n')
            f.write('\n'.join(missing_set))


def print_annotation_for_selected():
    cds_list = read_annotations(upstream_length, filter_by_gene_len=False)
    name_to_cds = {}
    for cds in cds_list:
        name_to_cds[cds.name] = cds
    for (dirpath, dirnames, filenames) in os.walk(path_to_features):
        for filename in filenames:
            i = filename.find('.selected')
            if i != -1:
                drug = filename[0: i]
                with open(out_path + drug + '.gene_annotation', 'w') as f:
                    genes = set()
                    for line in open(path_to_features + filename, 'r').readlines():
                        s = line.strip().split('\t')
                        if s[0] != 'non_cds':
                            genes.add(s[1])
                    f.write('name\tpseudogene\thypothetical\tin_proteom\tproduct\n')
                    for gene in genes:
                        cds = name_to_cds[gene]
                        f.write(gene + '\t')
                        if cds.is_pseudogene:
                            f.write('1\t')
                        else:
                            f.write('0\t')
                        if cds.is_hypothetical:
                            f.write('1\t')
                        else:
                            f.write('0\t')
                        if cds.exists_in_proteom:
                            f.write('1\t')
                        else:
                            f.write('0\t')
                        if cds.product is None:
                            f.write('-\n')
                        else:
                            f.write(cds.product)
                            f.write('\n')


if __name__ == '__main__':
    filter_selected()
    find_missing_from_dict()
    print_annotation_for_selected()