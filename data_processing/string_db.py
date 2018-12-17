from numpy.random.mtrand import choice

from src.core.annotations import read_annotations, CDSType
from src.core.constants import upstream_length, drug_names
from src.core.data_reading import read_dict

mt_db_path = '../../data/Mycobacterium_tuberculosis_H37Rv_txt_v3.txt'
db_path = '../../data/string_db/83332.protein.links.v10.5.txt'
path_to_dict = '../../data/dictionaries_indels/'
path_to_gene_annotations = '../../res/tree_was/annotations/'
drug_list = ['PZA', 'RIF', 'AMI', 'CAP', 'CIP', 'EMB', 'ETH', 'INH', 'KAN', 'MOX', 'OFX', 'PTH', 'STR']#
filter_genes = 'PGRS'
filter_pseudogenes = True
filter_hypothetical = True
filter_walker = True

permutation_num = 100


class Gene:

    def __init__(self, name, locus, function, category, product, walker=None):
        self.name = name
        self.locus = locus
        self.function = function
        self.category = category
        self.product = product
        self.walker = walker


def read_mtb_db():
    cds_list = read_annotations(upstream_length)
    name_to_cds = {cds.name: cds for cds in cds_list if cds.type != CDSType.upstream}
    locus_to_cds = {cds.locus: cds for cds in cds_list if cds.type != CDSType.upstream}
    name_to_gene = {}
    with open(mt_db_path, 'r') as f:
        for line in f.readlines()[1:]:
            s = line.strip().split('\t')
            cds = name_to_cds.get(s[9])
            if cds is None:
                cds = locus_to_cds.get(s[9])
            if cds is None:
                continue
            if filter_genes in s[9]:
                continue
            if filter_pseudogenes and cds.is_pseudogene:
                continue
            if filter_hypothetical and cds.is_hypothetical:
                continue
            if filter_hypothetical and 'hypothetical' in s[-2]:
                continue
            gene = Gene(s[9], s[8], s[10], s[15], s[11])
            name_to_gene[s[9]] = gene
    print('%d lines in mt db' % len(name_to_gene))
    return name_to_gene


def read_string_db()->'dict[str,list[str]]':
    locus_to_neighbors = {}
    for line in open(db_path).readlines()[1:]:
        s = line.strip().split(' ')
        loc1 = s[0].split('.')[1]
        loc2 = s[1].split('.')[1]
        if loc1 in locus_to_neighbors.keys():
            locus_to_neighbors[loc1].append(loc2)
        else:
            locus_to_neighbors[loc1] = [loc2]
    return locus_to_neighbors


def get_neighbors(locus_to_neighbors, locus_list):
    neighbors = set()
    for loc in locus_list:
        if loc in locus_to_neighbors:
            neighbors.update(locus_to_neighbors[loc])
    return neighbors


def read_Walker_dict()->'dict[str, list[str]]':
    """
    Reads Walker dictionary and returns map from drug name to gene list

    :return: dict[str, list[str]]
    """
    name, drug_to_mut_list = read_dict(path_to_dict, 'Walker_dictionary')
    drug_to_gen_list = {}
    for drug, m_list in drug_to_mut_list.items():
        gene_set = set()
        for m in m_list:
            gene_set.add(m.split('\t')[1])
        drug_to_gen_list[drug] = list(gene_set)
    return drug_to_gen_list


def read_gene_annotations():
    """
    Reads .annotations files for all drugs in drug_list and returns map from drug to map from locus name to Gene object

    :return: dict[str, dict[str, Gene]]
    """
    drug_to_genes = {}
    for drug in drug_list:
        locus_to_gene = {}
        drug_to_genes[drug] = locus_to_gene
        for line in open(path_to_gene_annotations + drug.lower() + '.upper.pvalue.pairs.annotations').readlines()[1:]:
            s = line.strip().split('\t')
            if filter_genes in s[0]:
                continue
            if filter_pseudogenes and s[2] == 'True':
                continue
            if filter_hypothetical and s[4] == 'True':
                continue
            if filter_hypothetical and 'hypothetical' in s[-2]:
                continue
            gene = Gene(s[0], s[1], s[-3], s[-2], s[-4], s[-1])
            if filter_walker and gene.walker != 'miss':
                continue
            locus_to_gene[s[1]] = gene
        print('%s %d lines in our annotation' % (drug, len(locus_to_gene)))
    return drug_to_genes


def get_Walker_neighbors()->'dict[str, list[str]]':
    drug_to_gen_list = read_Walker_dict()
    cds_list = read_annotations(upstream_length)
    name_to_cds = {}
    for cds in cds_list:
        if cds.type != CDSType.upstream:
            name_to_cds[cds.name] = cds
            if cds.synonym is not None:
                for syn in cds.synonym:
                    name_to_cds[syn] = cds
    locus_to_cds = {cds.locus: cds for cds in cds_list if cds.locus is not None}
    drug_to_neighbors = {}
    locus_to_neighbors = read_string_db()
    for drug, gene_list in drug_to_gen_list.items():
        locus_names = [name_to_cds[gene_name].locus for gene_name in gene_list]
        drug_to_neighbors[drug] = [locus_to_cds[locus].name for locus in get_neighbors(locus_to_neighbors, locus_names)
                                   if locus in locus_to_cds and locus not in locus_names]
    return drug_to_neighbors


def print_filtered_annotations():
    locus_to_annotation = read_gene_annotations()
    with open('../../res/tree_was/annotations/filtered.annotations', 'w') as f:
        f.write('gene\tlocus\tproduct\tfunction\tfunctional_category\tWalker\n')
        for locus, gene in locus_to_annotation.items():
            f.write('\t'.join([gene.name, gene.locus, gene.product, gene.function, gene.category, gene.walker]) + '\n')


def chi_sqr(func_freqs_all, func_freqs_our):
    res = 0
    for f, freq in func_freqs_all.items():
        if f in func_freqs_our:
            f_our = func_freqs_our[f]
            res += (freq - f_our)*(freq - f_our)
        else:
            res += freq*freq
    return res


def random_genes_stat(func_freqs_all, name_to_gene, n):
    subset = choice(list(name_to_gene.keys()), size=n)
    func_counts_our = {}
    for name in subset:
        gene = name_to_gene[name]
        if gene.category in func_counts_our.keys():
            func_counts_our[gene.category] += 1
        else:
            func_counts_our[gene.category] = 1
    func_freqs_our = {}
    for func, count in func_counts_our.items():
        func_freqs_our[func] = count/n
    return chi_sqr(func_freqs_all, func_freqs_our)


def func_cat_p_value(locus_to_gene, name_to_gene):
    func_counts_all = {}
    for name, gene in name_to_gene.items():
        if gene.category in func_counts_all.keys():
            func_counts_all[gene.category] += 1
        else:
            func_counts_all[gene.category] = 1
    func_freqs_all = {}
    for func, count in func_counts_all.items():
        func_freqs_all[func] = count/len(name_to_gene)
    func_counts_our = {}
    for name, gene in locus_to_gene.items():
        if gene.category in func_counts_our.keys():
            func_counts_our[gene.category] += 1
        else:
            func_counts_our[gene.category] = 1
    func_freqs_our = {}
    for func, count in func_counts_our.items():
        func_freqs_our[func] = count/len(locus_to_gene)
    chi = chi_sqr(func_freqs_all, func_freqs_our)
    c = 0
    for i in range(permutation_num):
        rand_score = random_genes_stat(func_freqs_all, name_to_gene, len(locus_to_gene))
        if chi <= rand_score:
            c += 1
    return c/permutation_num


def count_edges(loci, locus_to_neighbors):
    edge_num = 0
    for locus in loci:
        neighbors = locus_to_neighbors.get(locus)
        if neighbors is not None:
            edge_num += sum(loc in loci for loc in neighbors)
    return edge_num/2


def density_p_value():
    drug_to_genes = read_gene_annotations()
    locus_to_neighbors = read_string_db()
    gene_list = list(locus_to_neighbors.keys())
    for drug, locus_to_gene in drug_to_genes.items():
        genes = list(locus_to_gene.keys())
        edge_num = count_edges(genes, locus_to_neighbors)
        p_value = 0
        for i in range(permutation_num):
            random_genes = choice(gene_list, len(genes))
            if count_edges(random_genes, locus_to_neighbors) >= edge_num:
                p_value += 1
        print('%s %1.2f' % (drug, p_value/permutation_num))


def walker_neighbor_intersection_p_value():
    drug_to_genes = read_gene_annotations()
    drug_to_walker_neighbors = get_Walker_neighbors()
    drug_to_gen_list = read_Walker_dict()
    gene_list = [cds.name for cds in read_annotations(upstream_length)]
    for drug, locus_to_gene in drug_to_genes.items():
        drug = drug_names[drug]
        if drug not in drug_to_walker_neighbors:
            continue
        walker_neighbors = set(drug_to_walker_neighbors[drug])
        gene_set = set(gene.name for gene in locus_to_gene.values())
        stat = len(walker_neighbors.intersection(gene_set))
        p_value = 0
        gene_list_drug = [gene for gene in gene_list if gene not in drug_to_gen_list[drug]]
        for i in range(permutation_num):
            random_genes = set(choice(gene_list_drug, len(gene_set)))
            if len(walker_neighbors.intersection(random_genes)) >= stat:
                p_value += 1
        print('%s %1.2f' % (drug, p_value / permutation_num))


if __name__ == '__main__':
    # print_filtered_annotations()
    # locus_to_gene = read_gene_annotations()
    # name_to_gene = read_mtb_db()
    # print(func_cat_p_value(locus_to_gene, name_to_gene))
    # density_p_value()
    walker_neighbor_intersection_p_value()
