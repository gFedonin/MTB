from ete3 import Tree

path_to_tree = '../../data/tree_with_pheno_and_snp_mc5_rooted.nw'
out_path = '../../data/tree_with_pheno_and_snp_mc5.branch_lengths'

if __name__ == '__main__':
    t = Tree(path_to_tree)
    branch_lengths = []
    for node in t.traverse():
        branch_lengths.append(node.dist)
    with open(out_path, 'w') as f:
        for l in branch_lengths:
            f.write(str(l) + '\n')