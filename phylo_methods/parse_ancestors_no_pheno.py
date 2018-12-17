from ete3 import Tree
from os.path import exists
from sklearn.externals.joblib import Parallel, delayed
import os

# path_to_aln = './data/ancestors_merged.fasta'
from src.core.constants import data_path

path_to_node_map = data_path + 'ancestors_mc10_mega/0/split_0-14954_nodeMap.txt'
path_to_tree = data_path + 'tree_with_pheno_and_snp_mc5_mega_rooted.nw'
path_to_anc_indels_probs = data_path + 'indels_MP/'
path_to_sample_ids = data_path + 'dr_covered_with_pheno_and_snp.txt'
path_to_indels = data_path + 'indels/'

out_path = data_path + 'reconstructed_mc10_mega_MP/'


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
    sample_ids = [line.strip() for line in open(path_to_sample_ids, 'r').readlines()]
    t.prune(sample_ids)
    return t


def print_anc_pheno_by_drug(tree):
    with open(out_path + 'full_tree_parents.csv', 'w') as f:
        f.write(tree.get_tree_root().name + '\n')
        for node in tree.iter_descendants():
            f.write(node.name + '\t' + node.up.name + '\t' + str(node.dist) + '\n')
    return 0


def main():
    if not exists(out_path):
        os.makedirs(out_path)
    big_tree = read_tree()
    print_anc_pheno_by_drug(big_tree)


if __name__ == '__main__':
    main()

