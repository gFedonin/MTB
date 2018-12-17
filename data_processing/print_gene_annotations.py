from src.core.annotations import read_annotations, CDSType
from src.core.constants import upstream_length

drug_list = ['PZA', 'RIF', 'AMI', 'CAP', 'CIP', 'EMB', 'ETH', 'INH', 'KAN', 'MOX', 'OFL', 'PRO', 'STR']#
mt_db_path = '../../data/Mycobacterium_tuberculosis_H37Rv_txt_v3.txt'


def read_mtb_db():
    name_to_data = {}
    with open(mt_db_path, 'r') as f:
        for line in f.readlines()[1:]:
            s = line.strip().split('\t')
            name_to_data[s[9]] = (s[10], s[11], s[15], s[8])
    return name_to_data


class Gene:

    def __init__(self, name, locus, function, category, product, is_pseudogene, exists_in_proteom, is_hypothetical, walker=None):
        self.name = name
        self.locus = locus
        self.function = function
        self.category = category
        self.product = product
        self.walker = walker
        self.is_pseudogene = is_pseudogene
        self.exists_in_proteom = exists_in_proteom
        self.is_hypothetical = is_hypothetical


if __name__ == '__main__':
    cds_list = read_annotations(upstream_length)
    name_to_cds = {}
    for cds in cds_list:
        name_to_cds[cds.name] = cds
    name_to_data = read_mtb_db()
    for drug in drug_list:
        in_path = '../../res/tree_was/converted/' + drug.lower() + '.upper.pvalue.pairs.converted'
        out_path = '../../res/tree_was/annotations/' + drug.lower() + '.upper.pvalue.pairs.annotations'
        with open(out_path, 'w') as out_f:
            out_f.write('gene\tlocus\tis_pseudogene\texists_in_proteom\tis_hypothetical\tproduct\tfunction\tfunctional_category\tWalker\n')
            processed = {}
            with open(in_path, 'r') as in_f:
                for line in in_f.readlines()[1:]:
                    s = line.strip().split('\t')
                    if s[1] != 'non_cds':
                        if s[1] not in processed:
                            cds = name_to_cds[s[1]]
                            data = name_to_data[s[1]]
                            gene = Gene(s[1], data[3], data[0], data[2], cds.product, cds.is_pseudogene,
                                        cds.exists_in_proteom, cds.is_hypothetical, s[-1])
                            processed[s[1]] = gene
                        else:
                            gene = processed[s[1]]
                            if gene.walker is None:
                                gene.walker = s[-1]
            for gene in processed.values():
                out_f.write(
                    '%s\t%s\t%r\t%r\t%r\t%s\t%s\t%s\t%s\n' % (gene.name, gene.locus, gene.is_pseudogene, gene.exists_in_proteom,
                                                              gene.is_hypothetical, gene.product, gene.function, gene.category,
                                                              gene.walker))