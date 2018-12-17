from ete3 import Tree
from os.path import exists
from sklearn.externals.joblib import Parallel, delayed
import os

# path_to_aln = './data/ancestors_merged.fasta'
from src.core.constants import data_path

path_to_node_map = data_path + 'ancestors_mc10_mega/0/split_0-14954_nodeMap.txt'
path_to_tree = data_path + 'tree_with_pheno_and_snp_mc5_mega_rooted.nw'
path_to_pheno = data_path + 'pheno_mc5_mega_mix/'
path_to_anc_probs = data_path + 'anc_probs_mc5_mega_MP_mix/'
path_to_anc_indels_probs = data_path + 'indels_MP/'
path_to_sample_ids = data_path + 'dr_covered_with_pheno_and_snp.txt'
path_to_indels = data_path + 'indels/'

out_path = data_path + 'reconstructed_mc10_mega_MP_mix/'
use_indels = False


def parse_node_map():
    node_id_to_name = {}
    name_to_node_id = {}
    children_to_node_id = {}
    with open(path_to_node_map, 'r') as f:
        f.readline()
        for line in f.readlines():
            s = line.strip().split()
            node_id_to_name[s[1]] = s[0]
            name_to_node_id[s[0]] = s[1]
            if s[2] != '-':
                children_to_node_id[s[2]] = s[1]
                children_to_node_id[s[3]] = s[1]
    return name_to_node_id, node_id_to_name, children_to_node_id


def read_tree():
    name_to_node_id, node_id_to_name, children_to_node_id = parse_node_map()
    t = Tree(path_to_tree)
    for node in t.traverse("postorder"):
      if node.is_leaf():
          node.add_feature('node_id', name_to_node_id[node.name])
      else:
        node.add_feature('node_id', children_to_node_id[node.children[0].node_id])
        node.name = node_id_to_name[node.node_id]
    return t


def print_anc_pheno_by_drug(tree, drug):
    # t = tree.copy()
    samples = []
    with open(path_to_pheno + drug + '.pheno', 'r') as f:
        for line in f.readlines():
            samples.append(line.split('\t')[0])
    tree.prune(samples)
    anc_pheno = []
    with open(path_to_anc_probs + drug + '.probs', 'r') as f:
        for line in f.readlines()[1:]:
            s = line.strip().split('\t')
            if float(s[0]) > float(s[1]):
                anc_pheno.append('0')
            else:
                anc_pheno.append('1')
    if not exists(out_path + drug):
        os.makedirs(out_path + drug)
    # node_to_pheno = {}
    with open(out_path + drug + '/anc.pheno', 'w') as f:
        i = 0
        for node in tree.traverse("preorder"):
            if not node.is_leaf():
                # node_to_pheno[node.name] = anc_pheno[i]
                f.write(node.name + '\t' + anc_pheno[i] + '\n')
                i += 1
    # t.write(format=1, outfile=out_path + drug + "/tree.nw")
    with open(out_path + drug + '/parents.csv', 'w') as f:
        f.write(tree.get_tree_root().name + '\n')
        for node in tree.iter_descendants():
            f.write(node.name + '\t' + node.up.name + '\t' + str(node.dist) + '\n')
    return 0


# def print_anc_pheno_for_MDR_XDR(tree, drug_to_anc_pheno):
#     t = tree.copy()
#     samples = []
#     with open(path_to_pheno + 'MDR.pheno', 'r') as f:
#         for line in f.readlines():
#             samples.append(line.split('\t')[0])
#     t.prune(samples)
#     with open(out_path + 'MDR/anc.pheno', 'w') as f:
#         for node in t.traverse("preorder"):
#             if not node.is_leaf():
#                 isoniazid = drug_to_anc_pheno['Isoniazid'][node.name]
#                 rifampicin = drug_to_anc_pheno['Rifampicin'][node.name]
#                 if isoniazid == '1'and rifampicin == '1':
#                     f.write(node.name + '\t1\n')
#                 else:
#                     f.write(node.name + '\t0\n')


def print_anc_indels(tree):
    sample_ids = [line.strip() for line in open(path_to_sample_ids, 'r').readlines()]
    tree.prune(sample_ids)
    anc_indels = []
    indel_names = []
    for (dirpath, dirnames, filenames) in os.walk(path_to_anc_indels_probs):
        for filename in filenames:
            with open(path_to_anc_indels_probs + filename, 'r') as f:
                indels = []
                f.readline()
                for line in f.readlines():
                    s = line.strip().split('\t')
                    if float(s[0]) > float(s[1]):
                        indels.append('0')
                    else:
                        indels.append('1')
                indel_names.append(filename[:-6])
                anc_indels.append(indels)
    leaf_indels = {}
    for (dirpath, dirnames, filenames) in os.walk(path_to_indels):
        for filename in filenames:
            with open(path_to_indels + filename, 'r') as f:
                indels = {}
                for line in f.readlines():
                    s = line.strip().split('\t')
                    indels[s[0]] = s[1]
                leaf_indels[filename[:-8]] = indels
    with open(out_path + 'anc_indels.fasta', 'w') as f:
        i = 0
        for node in tree.traverse("preorder"):
            f.write('>' + node.name + '\n')
            if not node.is_leaf():
                for indels in anc_indels:
                    f.write(indels[i])
                f.write('\n')
                i += 1
            else:
                for indel_name in indel_names:
                    f.write(leaf_indels[indel_name][node.name])
                f.write('\n')
    with open(out_path + 'indel_list.txt', 'w') as f:
        for name in indel_names:
            f.write(name)
            f.write('\n')


def main():
    if not exists(out_path):
        os.makedirs(out_path)
    big_tree = read_tree()
    drugs = []
    for (dirpath, dirnames, filenames) in os.walk(path_to_pheno):
        for filename in filenames:
            drugs.append(filename[:-6])
            # prune(big_tree, filename[:-6])
    tasks = Parallel(n_jobs=-1)(delayed(print_anc_pheno_by_drug)(big_tree, drug) for drug in drugs)
    if use_indels:
        print_anc_indels(big_tree)
    c = 0
    for task in tasks:
        c += task
    print(c)


if __name__ == '__main__':
    main()

