path_to_tree = '../../data/xparr/mc10_mega_MP_RR_vars_vs_vars_filtered_frd0.1_pvalue0.05_pheno_array_10drugs.xparr'
node_name = 'Node_5437'
out_path = '../../data/xparr/mc10_mega_MP_RR_vars_vs_vars_filtered_frd0.1_pvalue0.05_pheno_array_10drugs_bad.xparr'


def count_resistant_in_subtree():
    parents = [node_name]
    child = 0
    total = 0
    resistant = 0
    for line in open(path_to_tree).readlines()[2:]:
        s = line.strip().split()
        if s[1] == parents[0]:
            if s[0][0] == 'S':
                total += 1
                if s[-1][0] == '1':
                    resistant += 1
            else:
                parents.append(s[0])
            if child == 0:
                child = 1
            else:
                child = 0
                parents.pop(0)
                if len(parents) == 0:
                    break
    print('%d from %d' % (resistant, total))


def print_subtree():
    parents = [node_name]
    child = 0
    with open(out_path, 'w') as fout:
        for line in open(path_to_tree).readlines()[2:]:
            s = line.strip().split()
            if s[1] == parents[0]:
                fout.write(line)
                if s[0][0] != 'S':
                    parents.append(s[0])
                if child == 0:
                    child = 1
                else:
                    child = 0
                    parents.pop(0)
                    if len(parents) == 0:
                        break



if __name__ == '__main__':
    print_subtree()