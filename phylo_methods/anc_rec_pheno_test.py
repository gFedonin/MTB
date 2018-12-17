from ete3 import Tree

path_to_tree = '../../R/test/tree.nw'
path_to_pheno = '../../R/test/pheno.txt'
path_to_anc = '../../R/test/anc.probs'
path_to_fig_tree = '../../data/fig_tree.txt'

out_path = '../../R/test/colored_tree.nex'

red_color = '[&!color=#ff0033]'

def main():
    t = Tree(path_to_tree)
    anc_pheno = []
    with open(path_to_anc, 'r') as f:
        f.readline()
        for line in f.readlines():
            s = line.strip().split('\t')
            if float(s[0]) > float(s[1]):
                anc_pheno.append(0)
            else:
                anc_pheno.append(1)
    sample_pheno = {}
    with open(path_to_pheno, 'r') as f:
        for line in f.readlines():
            s = line.strip().split('\t')
            sample_pheno[s[0]] = int(s[1])
    i = 0
    for node in t.traverse('levelorder'):
        node.add_feature('pheno', 0)
        if node.is_leaf():
            node.pheno = sample_pheno[node.name]
        else:
            node.name = 'node_' + str(i)
            node.pheno = anc_pheno[i]
            i += 1
    nw = t.write(format=1)
    node_num = 0
    for node in t.traverse():
        if node.pheno == 1:
            nw = nw.replace(node.name, node.name + red_color)
        node_num += 1

    with open(out_path, 'w') as f:
        f.write('#NEXUS\n')
        f.write('begin taxa;\n')
        f.write('\tdimensions ntax=' + str(node_num) + ';\n')
        f.write('\ttaxlabels\n')
        for n in t.traverse():
            f.write('\t' + n.name + '\n')
        f.write(';\n')
        f.write('end;\n\n')
        f.write('begin trees;\n')
        f.write('\ttree tree_1 = [&R] ' + nw)
        f.write('\nend;\n\n')
        for line in open(path_to_fig_tree, 'r').readlines():
            f.write(line)


if __name__ == '__main__':
    main()