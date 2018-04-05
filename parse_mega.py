import os
import numpy as np
from sklearn.externals.joblib import Parallel, delayed

path_to_mega_output = './data/ancestors_mc5/'
out_path = './data/temp/'

split_num = 144


def parse_single(split):
    node_id_to_name = {}
    node_id_to_children = {}
    node_id_to_seq = {}
    for (dirpath, dirnames, filenames) in os.walk(path_to_mega_output + str(split)):
        for filename in filenames:
            if 'nodeMap' in filename:
                with open(path_to_mega_output + str(split) + '/' + filename, 'r') as f:
                    f.readline()
                    for line in f.readlines():
                        s = line.strip().split()
                        node_id_to_name[s[1]] = s[0]
                        if s[2] != '-':
                            node_id_to_children[s[1]] = (s[2], s[3])
            elif '.csv' in filename:
                with open(path_to_mega_output + str(split) + '/' + filename, 'r') as f:
                    first = f.readline().strip()
                    s = first.split('Node_')
                    nodes = s[1:]
                    node_seq = []
                    for i in range(len(nodes)):
                        node_seq.append([])
                    pos = ''
                    for line in f.readlines():
                        s = line[:-1].split()
                        p = s[0][0:-2]
                        n = s[0][-1]
                        vals = np.zeros(len(nodes))
                        for i in range(1, len(s)):
                            vals[i - 1] = float(s[i])
                        if p != pos:
                            if pos != '':
                                for i in range(len(nodes)):
                                    node_seq[i].append(max_n[i])
                            pos = p
                            max_prob = np.empty(len(nodes), dtype=float)
                            max_n = np.empty(len(nodes), dtype=str)
                            for i in range(1, len(s)):
                                max_prob[i - 1] = float(s[i])
                                max_n[i - 1] = n
                        else:
                            for i in range(1, len(s)):
                                v = float(s[i])
                                if v > max_prob[i - 1]:
                                    max_prob[i - 1] = v
                                    max_n[i - 1] = n
                    for i in range(len(nodes)):
                        node_seq[i].append(max_n[i])
                        node_seq[i].append('\n')
                        node_id_to_seq[nodes[i]] = ''.join(node_seq[i])
    with open(out_path + 'split_' + str(split) + '.fasta', 'w') as f:
        for node_id, name in node_id_to_name.items():
            f.write('>')
            f.write(name)
            f.write('\n')
            f.write(node_id_to_seq[node_id])
    return 0


def main():
    tasks = Parallel(n_jobs=-1)(delayed(parse_single)(split) for split in range(split_num))
    c = 0
    for task in tasks:
        c += task
    print(c)




if __name__ == '__main__':
    main()
