path1 = '../../data/anc_probs_mc5_mega_ARD/Rifampicin.probs'
path2 = '../../data/anc_probs_mc5_mega_MP/Rifampicin.probs'


if __name__ == '__main__':
    anc1 = []
    for line in open(path1, 'r').readlines()[1:]:
        s = line.split('\t')
        if float(s[0]) > float(s[1]):
            anc1.append(0)
        else:
            anc1.append(1)
    anc2 = []
    for line in open(path2, 'r').readlines()[1:]:
        s = line.split('\t')
        if float(s[0]) > float(s[1]):
            anc2.append(0)
        else:
            anc2.append(1)
    c = 0
    for i in range(len(anc1)):
        if anc1[i] != anc2[i]:
            c += 1
    print(str(c))