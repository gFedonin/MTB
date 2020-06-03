path_to_tree = 'd:\\MyProjects\\Neverov\\MTB\\readme\\data\\Casali\\casali_fast_tree.nw'
path_to_map = 'd:\\MyProjects\\Neverov\\MTB\\readme\\data\\Casali\\PRJEB2138_AssemblyDetails.txt'
out_path = 'd:\\MyProjects\\Neverov\\MTB\\readme\\data\\Casali\\casali_fast_tree_converted.nw'


if __name__ == '__main__':
    tree = open(path_to_tree).readline()
    for l in open(path_to_map).readlines():
        if l[0] == '#':
            continue
        s = l.strip().split('\t')
        tree = tree.replace(s[3], s[4])
    with open(out_path, 'w') as f:
        f.write(tree)

