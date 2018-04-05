import os
from math import sqrt

from Bio.SubsMat import MatrixInfo
from os.path import exists
from sklearn.externals.joblib import delayed, Parallel

from constants import aminoacids

path_to_snp = './data/snps/annotated_with_pheno_and_snp_mc5/'
out_path = './data/snps/annotated_transformed_genes_with_pheno_and_snp_mc5/'
out_path_list = out_path + 'all_snp_list.csv'


def create_distance_matrix(similarity_matrix):
    dist_matrix = {}
    aa_to_vec = {}
    for aa1 in aminoacids:
        vec = []
        for aa2 in aminoacids:
            if (aa1, aa2) in similarity_matrix.keys():
                vec.append(similarity_matrix[(aa1, aa2)])
            else:
                vec.append(similarity_matrix[(aa2, aa1)])
        aa_to_vec[aa1] = vec
    for aa1 in aminoacids:
        vec1 = aa_to_vec[aa1]
        for aa2 in aminoacids:
            vec2 = aa_to_vec[aa2]
            dist = 0
            for i in range(len(aminoacids)):
                dist += (vec1[i] - vec2[i])**2
            dist_matrix[(aa1, aa2)] = sqrt(dist)
    return dist_matrix


def print_broken_genes(sample_id):
    broken_genes = []
    with open(path_to_snp + sample_id + '.snp', 'r') as f:
        for line in f.readlines():
            s = line.strip().split('\t')
            if s[0] == 'Gene' and (s[4] == '*' or s[3] == '*'):
                broken_genes.append(s[1])
    with open(out_path + sample_id + '.snp', 'w') as f:
        for gene in broken_genes:
            f.write(gene + '\n')
    return broken_genes


def print_all_broken_genes():
    if not exists(out_path):
        os.mkdir(out_path)
    for (dirpath, dirnames, filenames) in os.walk(path_to_snp):
        broken_genes_list = set()
        tasks = Parallel(n_jobs=-1)(delayed(print_broken_genes)(sample_id) for sample_id in [name[:-4] for name in filenames if '.snp' in name])
        for l in tasks:
            for gene in l:
                broken_genes_list.add(gene)
        with open(out_path_list, 'w') as f:
            for gene in broken_genes_list:
                f.write(gene + '\n')


def print_changed_genes(sample_id, distance_matrix):
    reformatted_snps = []
    with open(path_to_snp + sample_id + '.variants', 'r') as fin:
        with open(out_path + sample_id + '.variants', 'w') as fout:
            for line in fin.readlines():
                s = line.strip().split('\t')
                if s[0] == 'Gene':
                    if s[4] != '*' and s[3] != '*' and s[5] != 'FS':
                        fout.write(s[0] + '\t' + s[1] + '\t' + s[2] + '\t' + str(distance_matrix[(s[3], s[4])]) + '\n')
                    else:
                        fout.write(line)
                    reformatted_snps.append(s[0] + '\t' + s[1] + '\t' + s[2] + '\n')
                else:
                    fout.write(line)
                    reformatted_snps.append(line)
    return reformatted_snps


def print_all_changed_genes():
    if not exists(out_path):
        os.mkdir(out_path)
    matrix = create_distance_matrix(MatrixInfo.blosum62)
    for (dirpath, dirnames, filenames) in os.walk(path_to_snp):
        all_snps_list = set()
        tasks = Parallel(n_jobs=-1)(delayed(print_changed_genes)(sample_id, matrix) for sample_id in [name[:-4] for name in filenames if '.snp' in name])
        for l in tasks:
            for snp in l:
                all_snps_list.add(snp)
        with open(out_path_list, 'w') as f:
            for snp in all_snps_list:
                f.write(snp)


def main():
    print_all_changed_genes()


if __name__ == '__main__':
    main()