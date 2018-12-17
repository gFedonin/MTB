from ete3 import Tree

from src.core.constants import data_path

path_to_tree = data_path + 'tree_with_pheno_and_snp_mc5_mega.nw'
out_path = data_path + 'tree_with_pheno_and_snp_mc5_mega_rooted.nw'

t = Tree(path_to_tree)
t.set_outgroup('canetti')
t.write(outfile=out_path)
