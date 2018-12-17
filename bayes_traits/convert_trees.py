from Bio import Phylo
import os

path_to_trees = './data/trees/'
out_path = './data/trees_nexus/'

def newick_to_nexus():
    for (dirpath, dirnames, filenames) in os.walk(path_to_trees):
        for filename in filenames:
            Phylo.convert(path_to_trees + filename, 'newick', out_path + filename[:-3] + '.nexus', 'nexus')


def reformat_nexus():
    for (dirpath, dirnames, filenames) in os.walk(out_path):
        for filename in filenames:
            if 'nexus' not in filename:
                continue
            with open(out_path + filename, 'r') as f:
                f.readline()
                f.readline()
                f.readline()
                labels = f.readline()[:-2].split()[1:]
                f.readline()
                f.readline()
                tree = f.readline().strip()
                tree = tree.replace('Tree', 'tree')
                tree = tree.replace(';', ');')
                tree = tree.replace('tree1=', 'tree1 = (')
            with open(out_path + filename, 'w') as f:
                f.write('#NEXUS\n')
                f.write('begin trees;\n')
                f.write('\ttranslate\n')
                for i in range(len(labels)-1):
                    f.write('\t\t' + str(i + 1) + ' ' + labels[i] + ',\n')
                    tree = tree.replace(labels[i], str(i + 1))
                f.write('\t\t' + str(len(labels)) + ' ' + labels[-1] + ';\n')
                tree = tree.replace(labels[-1], str(len(labels)))
                f.write('\t\t' + tree + '\n')
                f.write('end;\n')


if __name__ == '__main__':
    newick_to_nexus()
    reformat_nexus()