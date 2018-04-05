from ete3 import Tree
from os.path import exists
from sklearn.externals.joblib import Parallel, delayed
import os

# path_to_aln = './data/ancestors_merged.fasta'
path_to_node_map = './data/ancestors_sparse_mc5/0/split_0-7037_nodeMap.txt'
path_to_tree = './data/tree_dr_covered_with_pheno_and_snp_mc5_rooted.nw'
path_to_pheno = './data/pheno_mc5_Walker/'
path_to_anc_probs = './data/anc_probs_mc5_Walker_MP/'
out_path = './data/reconstructed_mc5_Walker_MP/'

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
                if s[2] < s[3]:
                    children_to_node_id[s[2] + '_' + s[3]] = s[1]
                else:
                    children_to_node_id[s[3] + '_' + s[2]] = s[1]
    return name_to_node_id, node_id_to_name, children_to_node_id


def read_tree(name_to_node_id, node_id_to_name, children_to_node_id):
    t = Tree(path_to_tree)
    for node in t.traverse("postorder"):
      if node.is_leaf():
          node.add_feature('node_id', name_to_node_id[node.name])
      else:
        id0 = node.children[0].node_id
        id1 = node.children[1].node_id
        if id0 < id1:
            node.add_feature('node_id', children_to_node_id[id0 + '_' + id1])
        else:
            node.add_feature('node_id', children_to_node_id[id1 + '_' + id0])
        node.name = node_id_to_name[node.node_id]
    return t


def prune(tree, drug):
    t = tree.copy()
    samples = []
    with open(path_to_pheno + drug + '.pheno', 'r') as f:
        for line in f.readlines():
            samples.append(line.split('\t')[0])
    t.prune(samples)
    anc_pheno = []
    with open(path_to_anc_probs + drug + '.probs', 'r') as f:
        f.readline()
        for line in f.readlines():
            s = line.split('\t')
            if float(s[0]) > float(s[1]):
                anc_pheno.append('0')
            else:
                anc_pheno.append('1')
    if not os.path.exists(out_path + drug):
        os.makedirs(out_path + drug)
    with open(out_path + drug + '/anc.pheno', 'w') as f:
        i = 0
        for node in t.traverse("preorder"):
            if not node.is_leaf():
                f.write(node.name + '\t' + anc_pheno[i] + '\n')
                i += 1
    # t.write(format=1, outfile=out_path + drug + "/tree.nw")
    with open(out_path + drug + '/parents.csv', 'w') as f:
        f.write(t.get_tree_root().name + '\n')
        for node in t.iter_descendants():
            f.write(node.name + '\t' + node.up.name + '\t' + str(node.dist) + '\n')
    return 0


def main():
    if not exists(out_path):
        os.makedirs(out_path)

    name_to_node_id, node_id_to_name, children_to_node_id = parse_node_map()
    big_tree = read_tree(name_to_node_id, node_id_to_name, children_to_node_id)
    drugs = []
    for (dirpath, dirnames, filenames) in os.walk(path_to_pheno):
        for filename in filenames:
            drugs.append(filename[:-6])
            # prune(big_tree, filename[:-6])
    tasks = Parallel(n_jobs=-1)(delayed(prune)(big_tree, drug) for drug in drugs)
    c = 0
    for task in tasks:
        c += task
    print(c)

main()

