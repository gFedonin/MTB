import os

path = './data/ancestors_sparse/'

split_num = 144

for i in range(split_num):
    os.system('mv ' + path + str(i) + '/aln.fasta ' + path + str(i) + '/split_' + str(i) + '.fasta')