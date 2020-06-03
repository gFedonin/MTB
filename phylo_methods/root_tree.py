from ete3 import Tree

from core.constants import data_path

# path_to_tree = data_path + 'tree_with_pheno_and_snp_mc5_mega.nw'
# path_to_tree = data_path + '/iqtree_with_indels/combined_mq40_keep_complex_std_names_filtered_no_DR.partitions.treefile'
path_to_tree = data_path + '/iqtree/combined_mq40_keep_complex_std_names_filtered_no_DR.partitions.treefile'
# out_path = data_path + 'tree_with_pheno_and_snp_mc5_mega_rooted.nw'
# out_path = data_path + '/iqtree_with_indels/combined_mq40_keep_complex_std_names_filtered_no_DR_partitions_rooted.nw'
out_path = data_path + '/iqtree/combined_mq40_keep_complex_std_names_filtered_no_DR_partitions_rooted.nw'

t = Tree(path_to_tree)
t.set_outgroup('canetti')
t.write(outfile=out_path)
