import os

path_to_splits = './data/split_aln_sparse/'
path_to_mega = './data/megacc'
path_to_config = './data/ancestral_seqs_ML_nucleotide.mao'
path_to_tree = './data/big_tree_filtered2.newick'
out_path = './data/ancestors_sparse/'

split_num = 144


def main():
    for i in range(split_num):
        out = out_path + str(i) + '/'
        os.system('mkdir ' + out)
        os.system('nohup ' + path_to_mega + ' -a ' + path_to_config + ' -d ' + path_to_splits + 'split_' + str(i) +
                  '.fasta -f Fasta -t ' + path_to_tree + ' -o ' + out + ' > ' + out + 'out.txt 2>> ' + out + 'out.txt &')


if __name__ == '__main__':
    main()