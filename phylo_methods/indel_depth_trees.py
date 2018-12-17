from os import makedirs

from os.path import exists

from Bio import SeqIO
from ete3 import Tree

from src.core.constants import data_path
from src.phylo_methods.parse_ancestors import read_tree, parse_node_map

path_to_indel_aln = data_path + 'reconstructed_mc10_mega_MP/anc_indels.fasta'
path_to_indel_list = data_path + 'reconstructed_mc10_mega_MP/indel_list.txt'
path_to_sample_ids = data_path + 'dr_covered_with_pheno_and_snp_new.txt'
path_to_fig_tree = data_path + 'fig_tree.txt'

out_path = '../../res/indel_trees_MP/'


threshold = 50
red_color = '[&!color=#ff0033]'


def main():
    if not exists(out_path):
        makedirs(out_path)
    name_to_node_id, node_id_to_name, children_to_node_id = parse_node_map()
    big_tree = read_tree(name_to_node_id, node_id_to_name, children_to_node_id)
    sample_ids = [sample_id.strip() for sample_id in open(path_to_sample_ids, 'r').readlines()]
    big_tree.prune(sample_ids)
    sample_id_to_indels = {}

    indel_list = [s.strip() for s in open(path_to_indel_list, 'r').readlines()]
    indel_num = len(indel_list)

    sample_sequences = SeqIO.parse(open(path_to_indel_aln, 'r'), 'fasta')
    for seq in sample_sequences:
        sample_id_to_indels[seq.name] = seq.seq
    for node in big_tree.traverse():
        node.add_feature('indels', sample_id_to_indels[node.name])

    for node in big_tree.traverse('postorder'):
        node.add_feature('min_len', 0)
        node.add_feature('max_len', 0)
        if not node.is_leaf():
            node.min_len = min(node.children[0].min_len, node.children[1].min_len) + 1
            node.max_len = max(node.children[0].max_len, node.children[1].max_len) + 1
    for node in big_tree.traverse('preorder'):
        for i in range(indel_num):
            if node.indels[i] == '1' and node.max_len >= threshold:
                newick = node.write(format=1)
                node_num = 0
                for n in node.traverse():
                    if n.indels[i] == '1':
                        newick = newick.replace(n.name, n.name + red_color)
                    node_num += 1

                with open(out_path + node.name + '_' + indel_list[i] + '.nex', 'w') as f:
                    f.write('#NEXUS\n')
                    f.write('begin taxa;\n')
                    f.write('\tdimensions ntax=' + str(node_num) + ';\n')
                    f.write('\ttaxlabels\n')
                    for n in node.traverse():
                        f.write('\t' + n.name + '\n')
                    f.write(';\n')
                    f.write('end;\n\n')
                    f.write('begin trees;\n')
                    f.write('\ttree tree_1 = [&R] ' + newick)
                    f.write('\nend;\n\n')
                    for line in open(path_to_fig_tree, 'r').readlines():
                        f.write(line)


if __name__ == '__main__':
    main()