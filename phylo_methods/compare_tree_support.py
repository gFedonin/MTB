from ete3 import Tree
from numpy import mean

path_to_tree1 = '/export/data/fedonin/MTB/data/casali14_snp/casali_snps.fasta.contree'
path_to_tree2 = '/export/data/fedonin/MTB/data/casali14_with_indels/casali.partitions.contree'


if __name__ == '__main__':
    print(mean([n.support for n in Tree(path_to_tree1).traverse()]))
    print(mean([n.support for n in Tree(path_to_tree2).traverse()]))
