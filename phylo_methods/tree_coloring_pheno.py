import os
import random
from math import floor

from os.path import exists

from ete3 import Tree

from src.core.constants import data_path

drug = 'Ethambutol'
forward_tag = '[&!color=#ff0033]'
backward_tag = '[&!color=#00ff00]'

path_to_fig_tree = data_path + 'fig_tree.txt'
path_to_pheno_and_trees = data_path + 'reconstructed_mc10_mega_MP_mix/'
path_to_pheno = data_path + 'pheno_mc5_mega_mix/'

out_path = data_path + 'colored_trees_pheno_RS/'

color_only_changes = False
prune = False
random_fraction = 0.1

def read_parents(drug):
    parents = []
    with open(path_to_pheno_and_trees + drug + '/parents.csv', 'r') as f:
        root = f.readline().strip()
        for line in f.readlines():
            s = line.strip().split('\t')
            parents.append((s[0], s[1], s[2]))
    return root, parents


def read_pheno(drug):

    sample_to_pheno = {}
    with open(path_to_pheno + drug + '.pheno', 'r') as f:
        for line in f.readlines():
            s = line.strip().split('\t')
            if s[1] == '1':
                sample_to_pheno[s[0]] = 'R'
            else:
                sample_to_pheno[s[0]] = 'S'
    with open(path_to_pheno_and_trees + drug + '/anc.pheno', 'r') as f:
        for line in f.readlines():
            s = line.strip().split('\t')
            if s[1] == '1':
                sample_to_pheno[s[0]] = 'R'
            else:
                sample_to_pheno[s[0]] = 'S'
    return sample_to_pheno


def format_mut_lists(sample_to_pheno, parents):
    forward = set()
    backward = set()
    for node_id, parent_id, dist in parents:
        pheno = sample_to_pheno[node_id]
        parent_pheno = sample_to_pheno[parent_id]
        if color_only_changes and pheno == parent_pheno:
            continue
        if pheno == 'R':
            forward.add(node_id)
        else:
            backward.add(node_id)
    return forward, backward


def main():

    if not exists(out_path):
        os.makedirs(out_path)

    root, parents = read_parents(drug)

    sample_to_pheno = read_pheno(drug)

    forward, backward = format_mut_lists(sample_to_pheno, parents)
    print('forward ' + str(len(forward)) + ' backward ' + str(len(backward)))

    tree = Tree()
    tree.name = root
    node_name_to_node = {root: tree}
    # if root in forward:
    #     tree.name += forward_tag
    # if root in backward:
    #     tree.name += backward_tag
    for node_id, parent_id, dist in parents:
        parent = node_name_to_node[parent_id]
        node = Tree()
        node_name_to_node[node_id] = node
        node.name = node_id
        # if node_id in forward:
        #     node.name += forward_tag
        # if node_id in backward:
        #     node.name += backward_tag
        parent.add_child(node, node_id, dist)
    if prune:
        selected_ids = random.sample(sample_to_pheno.keys(), floor(random_fraction*len(sample_to_pheno)))
        tree.prune(selected_ids)
    newick = tree.write(format=1)
    for sample_id, pheno in sample_to_pheno.items():
        if sample_id in forward:
            newick = newick.replace(sample_id, sample_id + '_' + pheno + forward_tag)
        elif sample_id in backward:
            newick = newick.replace(sample_id, sample_id + '_' + pheno + backward_tag)
        else:
            newick = newick.replace(sample_id, sample_id + '_' + pheno)
    with open(out_path + drug + '.nex', 'w') as f:
        f.write('#NEXUS\n')
        f.write('begin taxa;\n')
        f.write('\tdimensions ntax=' + str(len(sample_to_pheno)) + ';\n')
        f.write('\ttaxlabels\n')
        for sample_id, pheno in sample_to_pheno.items():
            f.write('\t' + sample_id + '_' + pheno + '\n')
        f.write(';\n')
        f.write('end;\n\n')
        f.write('begin trees;\n')
        f.write('\ttree tree_1 = [&R] ' + newick)
        f.write('\nend;\n\n')
        for line in open(path_to_fig_tree, 'r').readlines():
            f.write(line)
    # tree.write(format=1, outfile=out_path + drug + '.nw')


if __name__ == '__main__':
    main()
