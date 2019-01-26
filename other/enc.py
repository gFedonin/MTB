in_path = '../../../Enc/data/primers.csv'
out_path = '../../../Enc/data/primers.fasta'

with open(out_path, 'w') as f:
    for line in open(in_path).readlines():
        s = line.split('\t')
        f.write('>' + s[0].replace(' ', '-') + '\n')
        f.write(s[1].replace('-', '').replace(' ', '').upper())