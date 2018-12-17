from Bio import SeqIO
from ete3 import Tree
import matplotlib.pyplot as plt

from src.core.constants import data_path
from src.phylo_methods.parse_ancestors import read_tree, parse_node_map

path_to_indel_aln = data_path + 'reconstructed_mc10_mega_MP/anc_indels.fasta'
out_path_min = '../../res/indel_hist_min_MP.png'
out_path_max = '../../res/indel_hist_max_MP.png'
out_path_min_path = '../../res/indel_hist_min_path_MP.png'
out_path_max_path = '../../res/indel_hist_max_path_MP.png'
path_to_sample_ids = data_path + 'dr_covered_with_pheno_and_snp_new.txt'


def main():

    name_to_node_id, node_id_to_name, children_to_node_id = parse_node_map()
    big_tree = read_tree(name_to_node_id, node_id_to_name, children_to_node_id)
    sample_ids = [sample_id.strip() for sample_id in open(path_to_sample_ids, 'r').readlines()]
    big_tree.prune(sample_ids)
    min_depths = []
    max_depths = []
    min_indel_path_lens = []
    max_indel_path_lens = []
    sample_id_to_indels = {}
    sample_sequences = SeqIO.parse(open(path_to_indel_aln, 'r'), 'fasta')
    indel_num = 0
    for seq in sample_sequences:
        sample_id_to_indels[seq.name] = seq.seq
        indel_num = len(seq.seq)
    for node in big_tree.traverse():
        node.add_feature('indels', sample_id_to_indels[node.name])

    for node in big_tree.traverse('postorder'):
        node.add_feature('min_len', 0)
        node.add_feature('max_len', 0)
        node.add_feature('min_indel_path_len', {})
        node.add_feature('max_indel_path_len', {})
        if not node.is_leaf():
            node.min_len = min(node.children[0].min_len, node.children[1].min_len) + 1
            node.max_len = max(node.children[0].max_len, node.children[1].max_len) + 1
            for i in range(indel_num):
                if node.indels[i] == '1':
                    ch0_len = node.children[0].min_indel_path_len.get(i)
                    ch1_len = node.children[1].min_indel_path_len.get(i)
                    if ch0_len is None:
                        if ch1_len is None:
                            node.min_indel_path_len[i] = 0
                            node.max_indel_path_len[i] = 0
                        else:
                            node.min_indel_path_len[i] = ch1_len + 1
                            node.max_indel_path_len[i] = node.children[1].max_indel_path_len[i] + 1
                    else:
                        if ch1_len is None:
                            node.min_indel_path_len[i] = ch0_len + 1
                            node.max_indel_path_len[i] = node.children[0].max_indel_path_len[i] + 1
                        else:
                            node.min_indel_path_len[i] = min(ch0_len, ch1_len) + 1
                            node.max_indel_path_len[i] = \
                                max(node.children[0].max_indel_path_len[i], node.children[1].max_indel_path_len[i]) + 1
        else:
            for i in range(indel_num):
                if node.indels[i] == '1':
                    node.min_indel_path_len[i] = 0
                    node.max_indel_path_len[i] = 0
    for node in big_tree.traverse('postorder'):
        for i in range(indel_num):
            if node.indels[i] == '1':
                min_depths.append(node.min_len)
                max_depths.append(node.max_len)
                min_indel_path_lens.append(node.min_indel_path_len[i])
                max_indel_path_lens.append(node.max_indel_path_len[i])
    plt.switch_backend('agg')
    plt.title('Histogram of min depths')
    plt.xlabel('Min branch len to the tip')
    plt.ylabel('Percent of nodes')
    n, bins, patches = plt.hist(min_depths, 50, density=True, facecolor='g', alpha=0.75)
    plt.savefig(out_path_min)
    plt.clf()
    plt.title('Histogram of max depths')
    plt.xlabel('Max branch len to the tip')
    plt.ylabel('Percent of nodes')
    n, bins, patches = plt.hist(max_depths, 50, density=True, facecolor='g', alpha=0.75)
    plt.savefig(out_path_max)
    plt.title('Histogram of min path lens')
    plt.xlabel('Min path len to the tip')
    plt.ylabel('Percent of nodes')
    n, bins, patches = plt.hist(min_indel_path_lens, 50, density=True, facecolor='g', alpha=0.75)
    plt.savefig(out_path_min_path)
    plt.clf()
    plt.title('Histogram of max path lens')
    plt.xlabel('Max path len to the tip')
    plt.ylabel('Percent of nodes')
    n, bins, patches = plt.hist(max_indel_path_lens, 50, density=True, facecolor='g', alpha=0.75)
    plt.savefig(out_path_max_path)


if __name__ == '__main__':
    main()