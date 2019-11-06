from os import makedirs
from os.path import exists

from core.constants import data_path
from core.data_reading import read_dict
from ete3 import Tree
from phylo_methods.print_XPARR_variants_and_pheno import read_pheno

path_to_pheno = data_path + 'pheno_mc5_mega/'
path_to_pheno_and_trees = data_path + 'reconstructed_mc10_mega_MP/'
path_to_ids = data_path + '10drugs.sample_list'
path_to_dictionaries = data_path + 'dictionaries/'

out_path_pheno = data_path + '10drugs_pheno.csv'
out_path_tree = data_path + '10drugs_with_parent_names.nw'
out_path_walker_genes = data_path + 'all_walker_genes.list'


drug_names = ('Isoniazid', 'Rifampicin', 'Ethambutol', 'Pyrazinamide', 'Streptomycin', 'Moxifloxacin', 'Ofloxacin',
                  'Amikacin', 'Capreomycin', 'Kanamycin')


def create_tree(path_to_pheno_and_trees, drug, sample_ids):
    # parents = []
    t = Tree()
    name_to_node = {}
    with open(path_to_pheno_and_trees + drug + '/parents.csv', 'r') as f:
        t.name = f.readline().strip()
        name_to_node[t.name] = t
        for line in f.readlines():
            s = line.strip().split('\t')
            parent = name_to_node[s[1]]
            child = parent.add_child(name=s[0], dist=float(s[2]))
            name_to_node[s[0]] = child
    t.prune(sample_ids)
    # for node in t.iter_descendants("levelorder"):
    #     parents.append((node.name, node.up.name, str(node.dist)))
    return t


def extract_dr_genes(path_to_dictionaries):
    name, drug_to_mut_list = read_dict(path_to_dictionaries, 'Walker_dictionary')
    drug_to_gene_set = {}
    all_genes = set()
    for drug, mut_list in drug_to_mut_list.items():
        genes = set()
        for mut in mut_list:
            s = mut.split('\t')
            if s[0] != 'rRNA':
                genes.add(s[1])
                all_genes.add(s[1])
        drug_to_gene_set[drug] = genes
    return drug_to_gene_set, all_genes


def print_tree_data():
    sample_ids = [sample_id.strip() for sample_id in open(path_to_ids, 'r').readlines()]
    drug_to_pheno = {}
    drug_to_gene_set, all_genes = extract_dr_genes(path_to_dictionaries)
    with open(out_path_walker_genes, 'w') as f:
        for drug in drug_names:
            sample_to_pheno = read_pheno(path_to_pheno, path_to_pheno_and_trees, drug)
            drug_to_pheno[drug] = sample_to_pheno
            for gene in drug_to_gene_set[drug]:
                f.write('%s\t%s\n' % (drug, gene))        
    tree = create_tree(path_to_pheno_and_trees, drug_names[0], sample_ids)
    with open(out_path_pheno, 'w') as f:
        f.write('NodeID')
        for drug in drug_names:
            f.write('\t' + drug)
        f.write('\n')
        for node in tree.iter_descendants("levelorder"):
            f.write(node.name)
            for drug in drug_names:
                f.write('\t' + drug_to_pheno[drug][node.name])
            f.write('\n')
    tree.write(format=1, outfile=out_path_tree)



if __name__ == "__main__":
    print_tree_data()
