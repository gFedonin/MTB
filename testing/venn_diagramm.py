from os import listdir, makedirs
from os.path import exists

from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_circles

from src.core.constants import data_path

path_to_epistat = '../../res/tree_was/10_drug_vs_vars_0.1.csv'
path_to_lr = '../../res/burkina/filtered_features/logreg/intersection/'
path_to_sgimc = '../../res/burkina/filtered_features/SGIMC/intersection/'
out_path = '../../res/venn/'
path_to_dict = data_path + 'dictionaries_indels/Walker_dictionary.txt'


def read_Walker():
    walker_genes = set()
    with open(path_to_dict, 'r') as f:
        for line in f.readlines():
            s = line.split('\t')
            walker_genes.add(s[1])
    return walker_genes


def read_burkina(path):
    walker_genes = read_Walker()
    drug_to_gene_set = {}
    for drug in listdir(path):
        fin = open(path + drug)
        line = fin.readline()
        gene_set = set()
        while line != '':
            name = line.split('\t')[1]
            line = fin.readline()
            all_positive = True
            while line != '\n':
                if float(line) <= 0:
                    all_positive = False
                line = fin.readline()
            if all_positive and name != '-' and 'PGRS' not in name and name not in walker_genes:
                gene_set.add(name)
            line = fin.readline()
        drug_to_gene_set[drug] = gene_set
    return drug_to_gene_set


def read_epistat():
    walker_genes = read_Walker()
    epistat = set()
    for line in open(path_to_epistat).readlines()[1:]:
        name = line.split('\t')[1]
        if name not in walker_genes:
            epistat.add(name)
    return epistat


if __name__ == '__main__':
    if not exists(out_path):
        makedirs(out_path)
    epistat = read_epistat()
    lr_drug_to_gene = read_burkina(path_to_lr)
    sgimc_drug_to_gene = read_burkina(path_to_sgimc)
    for drug, lr_gene_set in lr_drug_to_gene.items():
        sgicm_gene_set = sgimc_drug_to_gene[drug]
        v = venn3([epistat, lr_gene_set, sgicm_gene_set], ('epistat6', 'log_reg', 'SGIMC'))
        plt.title(drug + " Venn diagram")
        plt.savefig(out_path + drug + '.png')
        plt.clf()