from ete3 import Tree
import os

from os.path import exists
from sklearn.externals.joblib import Parallel, delayed

path_to_tree = './data/tree_dr_covered_with_pheno_and_snp_mc5_rooted.nw'
path_to_tree_breaker_output = './data/tree_breaker_mc5_Walker/'
out_path = './data/tree_features_mc5_Walker/'

min_posterior = 0.5


def gen_features_for_drug(filename, big_tree):
    features = []
    drug = filename.split('.')[0]
    with open(path_to_tree_breaker_output + filename, 'r') as f:
        marked_tree = Tree(f.readlines()[-1], format=1)
        for node in filter(lambda n: float(n.name.split('|')[-1].split('=')[1][:-1]) > min_posterior, marked_tree.traverse()):
            nodes = []
            for n in node.get_leaves():
                nodes.append(n.name.split('{')[0])
            # nodes = list(filter(lambda n: n.name in nodes, big_tree.get_leaves()))
            nodes = big_tree.get_common_ancestor(nodes).get_leaves()
            features.append(','.join([n.name for n in nodes]))
            # features.append(','.join(nodes))
    with open(out_path + drug + '.features', 'w') as f:
        for node in features:
            f.write(node + '\n')
    return 1


def main():
    if not exists(out_path):
        os.mkdir(out_path)
    big_tree = Tree(path_to_tree)
    big_tree.set_outgroup('canetti')

    for (dirpath, dirnames, filenames) in os.walk(path_to_tree_breaker_output):
        tasks = Parallel(n_jobs=-1)(delayed(gen_features_for_drug)(filename, big_tree) for filename in filenames)
        # for filename in filenames:
        #     features = []
        #     drug = filename.split('.')[0]
        #     with open(path_to_tree_breaker_output + filename, 'r') as f:
        #         marked_tree = Tree(f.readlines()[-1])
        #         for node in filter(lambda n: n.posterior > min_posterior, marked_tree.traverse()):
        #             features.append(big_tree.get_common_ancestor(node.get_leaves()).name)
        #     with open(out_path + drug + '.features', 'w') as f:
        #         for node in features:
        #             f.write(node + '\n')
        sum = 0
        for c in tasks:
            sum += c
        print('generated features for ' + str(sum) + ' drugs')

if __name__ == '__main__':
    main()
