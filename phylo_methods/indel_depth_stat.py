from os import makedirs

from os.path import exists

from Bio import SeqIO
from ete3 import Tree
from sklearn.externals.joblib import Parallel, delayed

from src.core.constants import data_path
from src.phylo_methods.parse_ancestors import read_tree, parse_node_map

path_to_indel_aln = data_path + 'reconstructed_mc10_mega_ARD/anc_indels.fasta'
path_to_indel_list = data_path + 'reconstructed_mc10_mega_ARD/indel_list.txt'
path_to_sample_ids = data_path + 'dr_covered_with_pheno_and_snp.txt'

out_path = '../../res/indel_stat_ARD.txt'


def count_errors(indel_n_list, node_indels, ch0_indels, ch1_indels):
    count = []
    rcount = []
    for i in indel_n_list:
        if node_indels[i] == '1':
            if ch0_indels[i] == '0' and ch1_indels[i] == '0':
                count.append(i)
        else:
            if ch0_indels[i] == '1' and ch1_indels[i] == '1':
                rcount.append(i)

    return count, rcount


thread_num = 144


def count_err_in_all_indels(node, indel_num):
    count = []
    rcount = []
    indel_n_lists = []
    for i in range(thread_num):
        indel_n_lists.append([])
    for i in range(indel_num):
        indel_n_lists[i%thread_num].append(i)
    tasks = Parallel(n_jobs=thread_num)(delayed(count_errors)(indel_n_list, node.indels, node.children[0].indels,
                                                              node.children[1].indels) for indel_n_list in indel_n_lists)
    for cnt, rcnt in tasks:
        count.extend(cnt)
        rcount.extend(rcnt)
    return count, rcount


def main():
    big_tree = read_tree()
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

    indel_n_to_count = {}
    indel_n_to_rcount = {}
    for node in big_tree.traverse():
        if not node.is_leaf():
            for i in range(indel_num):
                if node.indels[i] == '1':
                    if node.children[0].indels[i] == '0' and node.children[1].indels[i] == '0':
                        c = indel_n_to_count.get(i)
                        if c is None:
                            indel_n_to_count[i] = 1
                        else:
                            indel_n_to_count[i] = c + 1
                else:
                    if node.children[0].indels[i] == '1' and node.children[1].indels[i] == '1':
                        c = indel_n_to_rcount.get(i)
                        if c is None:
                            indel_n_to_rcount[i] = 1
                        else:
                            indel_n_to_rcount[i] = c + 1
    with open(out_path, 'w') as f:
        for i in range(indel_num):
            c1 = indel_n_to_count.get(i)
            c2 = indel_n_to_rcount.get(i)
            if c1 is None:
                if c2 is not None:
                    f.write(indel_list[i])
                    f.write('\t0\t')
                    f.write(str(c2))
                    f.write('\n')
            else:
                if c2 is None:
                    f.write(indel_list[i])
                    f.write('\t')
                    f.write(str(c1))
                    f.write('\t0\n')
                else:
                    f.write(indel_list[i])
                    f.write('\t')
                    f.write(str(c1))
                    f.write('\t')
                    f.write(str(c2))
                    f.write('\n')


if __name__ == '__main__':
    main()