from src.core.annotations import read_annotations
from src.core.constants import data_path, ref_len, upstream_length
from ete3 import Tree

path_to_xparr = data_path + 'xparr/mc10_mega_MP/'
out_path = data_path + 'colored_trees_xparr/'
path_to_pheno = data_path + 'pheno_mc5_mega/'
path_to_pheno_and_trees = data_path + 'reconstructed_mc10_mega_MP/'
drug = 'Rifampicin'
mut = 'Gene\trpoB\t450'
red = '&!color=#ff0033'
violet = '&!color=#cc00ff'
green = '&!color=#00ff00'
blue = '&!color=#0066ff'
path_to_fig_tree = data_path + 'fig_tree.txt'


def get_target_pos():
    s = mut.split('\t')
    target_gene = s[1]
    target_pos = int(s[2])
    cds_list = read_annotations(upstream_length)
    for cds in cds_list:
        if cds.name == target_gene:
            if cds.strand == 1:
                return cds.start + (target_pos - 1)*3
            else:
                return cds.end - (target_pos - 1)*3 - 2


def read_tree():
    mut_with_pheno_change = set()
    mut_with_no_pheno_change = set()
    no_mut_with_pheno_change = set()
    node_to_nonsyn = {}
    target_pos = get_target_pos()
    prot_pos = mut.split('\t')[2]
    print('target pos = ' + str(target_pos))
    fin = open(path_to_xparr + drug + '.xparr')
    fin.readline()
    tree = Tree()
    tree.name = fin.readline().strip()
    node_name_to_node = {tree.name: tree}
    c = 0
    for line in fin.readlines():
        s = line.strip().split('\t')
        node = Tree()
        node.name = s[0]
        node_name_to_node[node.name] = node
        parent = node_name_to_node[s[1]]
        node.up = parent
        parent.add_child(node, s[0], s[2])
        pheno_change = ''
        if len(s) >= 5:
            pheno_change = s[4]
        mut_found = False
        if len(s) == 6:
            nonsyn = s[5].split(';')
            snps = 0
            for m in nonsyn:
                if not m[0].isdigit():
                    pos_str = m[1:-1]
                    snps += 1
                    pos = int(pos_str)
                    if pos == target_pos:
                        c += 1
                        if pheno_change == '':
                            mut_with_no_pheno_change.add(node.name)
                        else:
                            mut_with_pheno_change.add(node.name)
                        mut_found = True
                        node_to_nonsyn[node.name] = m[0] + prot_pos + m[-1]
                        break
        if not mut_found and pheno_change != '':
            no_mut_with_pheno_change.add(node.name)
    print('target pos found ' + str(c))
    return tree, mut_with_pheno_change, mut_with_no_pheno_change, no_mut_with_pheno_change, node_to_nonsyn


def print_tree(tree, sample_to_pheno, mut_with_pheno_change, mut_with_no_pheno_change, no_mut_with_pheno_change, node_to_nonsyn):
    newick = tree.write(format=1)
    for node_id, pheno in sample_to_pheno.items():
        if node_id in mut_with_pheno_change:
            newick = newick.replace(node_id, node_id + '_' + pheno + '[' + green + ',mut=' + node_to_nonsyn[node_id] + ']')
        elif node_id in mut_with_no_pheno_change:
            newick = newick.replace(node_id, node_id + '_' + pheno + '[' + red + ',mut=' + node_to_nonsyn[node_id] + ']')
        elif node_id in no_mut_with_pheno_change:
            newick = newick.replace(node_id, node_id + '_' + pheno + '[' + blue + ']')
        else:
            newick = newick.replace(node_id, node_id + '_' + pheno)
    with open(out_path + drug + '.nex', 'w') as f:
        f.write('#NEXUS\n')
        f.write('begin taxa;\n')
        f.write('\tdimensions ntax=' + str(len(sample_to_pheno)) + ';\n')
        f.write('\ttaxlabels\n')
        for node_id, pheno in sample_to_pheno.items():
            f.write('\t' + node_id + '_' + pheno + '\n')
        f.write(';\n')
        f.write('end;\n\n')
        f.write('begin trees;\n')
        f.write('\ttree tree_1 = [&R] ' + newick)
        f.write('\nend;\n\n')
        for line in open(path_to_fig_tree, 'r').readlines():
            f.write(line)


def read_pheno(drug):
    sample_to_pheno = {}
    with open(path_to_pheno + drug + '.pheno', 'r') as f:
        for line in f.readlines():
            s = line.strip().split('\t')
            if s[1] == '1':
                sample_to_pheno[s[0]] = 'R'
            else:
                sample_to_pheno[s[0]] = 'S'
    with open(path_to_pheno_and_trees + drug + '/anc.pheno', 'r') as f:
        for line in f.readlines():
            s = line.strip().split('\t')
            if s[1] == '1':
                sample_to_pheno[s[0]] = 'R'
            else:
                sample_to_pheno[s[0]] = 'S'
    return sample_to_pheno


if __name__ == '__main__':
    tree, mut_with_pheno_change, mut_with_no_pheno_change, no_mut_with_pheno_change, node_to_nonsyn = read_tree()
    sample_to_pheno = read_pheno(drug)
    print_tree(tree, sample_to_pheno, mut_with_pheno_change, mut_with_no_pheno_change, no_mut_with_pheno_change, node_to_nonsyn)