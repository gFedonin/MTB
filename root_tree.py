from ete3 import Tree

path_to_tree = './data/tree_with_pheno_and_snp_mc5.nw'
out_path = './data/tree_with_pheno_and_snp_mc5_rooted.nw'

t = Tree(path_to_tree)
t.set_outgroup('canetti')
t.write(outfile=out_path)
