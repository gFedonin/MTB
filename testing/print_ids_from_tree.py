from src.core.constants import data_path
from ete3 import Tree

path_to_tree = data_path + 'tree_with_pheno_and_snp_mc5_mega_rooted.nw'
out_path = data_path + 'tree_with_pheno_and_snp_mc5_mega_rooted.ids'


if __name__ == '__main__':
    t = Tree(path_to_tree)
    ids = [node.name for node in t.traverse() if node.is_leaf()]
    ids.sort()
    with open(out_path, 'w') as f:
        for id in ids:
            f.write(id)
            f.write('\n')