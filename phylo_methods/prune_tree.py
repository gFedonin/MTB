from ete3 import Tree

from core.constants import data_path

path_to_tree = data_path + 'tree_with_pheno_and_snp_mc5_mega_rooted.nw'
path_to_sample_ids = data_path + 'dr_covered_with_pheno_and_snp.txt'
out_path = data_path + 'tree_with_pheno_and_snp_mc5_mega_rooted_pruned.nw'


if __name__ == '__main__':
    t = Tree(path_to_tree)
    sample_ids = [line.strip() for line in open(path_to_sample_ids, 'r').readlines()]
    t.prune(sample_ids)
    t.write(format=1, outfile=out_path)
