from bisect import bisect_left, bisect_right
from enum import Enum

from src.core.constants import data_path, upstream_length

path_to_annotations = data_path + 'AL123456_rev.gff'


class CDSType(Enum):
    Gene = 0
    ncRNA = 1
    rRNA = 2
    tRNA = 3
    miscRNA = 4
    upstream = 5


class CDS:
    # type: CDSType
    # start: int  # 1-based genome coordinate
    # end: int
    # strand: int  # 1 or -1
    # name: str
    # synonym: list
    # is_pseudogene: bool
    # product: str
    # exists_in_proteom: bool
    # is_hypothetical: bool

    def __init__(self, type: CDSType, start: int, end: int, strand: int, name: str, synonym: list=None, desc: str=None,
                 locus: str=None):
        self.name: str = name
        self.start: int = start
        self.end: int = end
        self.strand: int = strand
        self.type: CDSType = type
        self.synonym: list = synonym
        self.locus: str = locus
        if desc is not None:
            self.is_pseudogene: bool = 'pseudogene' in desc
            i = desc.find('product')
            if i != -1:
                j = desc.find('\"', i + 7)
                k = desc.find('\"', j + 1)
                self.product: str = desc[j + 1: k]
            else:
                self.product: str = None
            self.exists_in_proteom: bool = 'identified in proteomics study' in desc
            self.is_hypothetical: bool = 'Hypothetical' in desc or 'hypothetical' in desc
        else:
            # if self.type != CDSType.upstream:
            #     print('no desc for ' + name)
            self.is_pseudogene: bool = False
            self.product: str = None
            self.exists_in_proteom: bool = False
            self.is_hypothetical: bool = False

    def __lt__(self, other):
        return self.start < other.start


def read_annotations(upstream_length, filter_by_gene_len=True):
    genes = {}
    nc_rna = {}
    r_rna = {}
    t_rna = {}
    misc_rna = {}
    upstream = {}
    cds_desc = {}
    with open(path_to_annotations, 'r') as f:
        for line in f.readlines()[1:]:
            s = line.strip().split('\t')
            type = s[2]
            start = int(s[3])
            end = int(s[4])
            if type in ('Repeat_region', 'mobile_element'):
                continue
            strand = 1 if s[6] == '+' else -1
            ns = s[8].split(';')
            name = ns[0].split(' ')[1]
            syn = None
            if 'locus_tag' in ns[-1]:
                locus = ns[-1].split(' ')[-1]
            else:
                locus = None
            if type == 'Gene':
                if len(ns) >= 1:
                    for name_tag in ns[1:]:
                        if 'gene_synonym' in name_tag:
                            if syn is None:
                                syn = []
                            syn.append(name_tag.split(' ')[2])
                genes[name] = (start, end, strand, syn, locus)
                if strand == 1:
                    upstream[name] = (start - upstream_length, start - 1, strand, syn, locus)
                else:
                    upstream[name] = (end + 1, end + upstream_length, strand, syn, locus)
            elif type == 'ncRNA':
                if len(ns) >= 1:
                    for name_tag in ns[1:]:
                        if 'gene_synonym' in name_tag:
                            if syn is None:
                                syn = []
                            syn.append(name_tag.split(' ')[2])
                nc_rna[name] = (start, end, strand, syn, locus)
            elif type == 'rRNA':
                r_rna[name] = (start, end, strand)
            elif type == 'tRNA':
                t_rna[name] = (start, end, strand)
            elif type == 'Misc._RNA':
                misc_rna[name] = (start, end, strand)
            elif type == 'CDS':
                desc = s[8].split(';')
                ns = desc[2].split()
                name = ns[1]
                # print('found desc for ' + name)
                cds_desc[name] = line

    cds_list = []
    for name, value in r_rna.items():
        start, end, strand, syn, locus = genes.pop(name, None)
        cds_list.append(CDS(CDSType.rRNA, start, end, strand, name, syn, cds_desc.get(name), locus))
    for name, value in nc_rna.items():
        if name in genes:
            start, end, strand, syn, locus = genes.pop(name, None)
            cds_list.append(CDS(CDSType.ncRNA, start, end, strand, name, syn, cds_desc.get(name), locus))
        else:
            cds_list.append(CDS(CDSType.ncRNA, value[0], value[1], value[2], name, desc=cds_desc.get(name)))
    for name, value in t_rna.items():
        if name in genes:
            start, end, strand, syn, locus = genes.pop(name, None)
            cds_list.append(CDS(CDSType.tRNA, start, end, strand, name, syn, cds_desc.get(name), locus))
        else:
            cds_list.append(CDS(CDSType.tRNA, value[0], value[1], value[2], name, desc=cds_desc.get(name)))
    for name, value in misc_rna.items():
        if name in genes:
            start, end, strand, syn, locus = genes.pop(name, None)
            cds_list.append(CDS(CDSType.miscRNA, start, end, strand, name, syn, cds_desc.get(name), locus))
        else:
            cds_list.append(CDS(CDSType.miscRNA, value[0], value[1], value[2], name, desc=cds_desc.get(name)))
    wrongly_annotated = set()
    for name, value in genes.items():
        start = value[0]
        end = value[1]
        if not filter_by_gene_len or (end - start + 1) % 3 == 0:
            cds_list.append(CDS(CDSType.Gene, start, end, value[2], name, value[3], cds_desc.get(name), value[4]))
        else:
            wrongly_annotated.add(name)
            print(name + ' len % 3 != 0')
    for name, value in upstream.items():
        if name not in wrongly_annotated:
            cds_list.append(CDS(CDSType.upstream, value[0], value[1], value[2], name, value[3]))
    cds_list.sort()
    return cds_list


def localize_all_variants(pos_list: 'list[int]', cds_list: 'list[CDS]', keep_genes_set: 'set[str]'=None) -> 'dict[int, CDS]':
    """
    Finds CDS for each gene coord from all_snps list.

    :param pos_list: list of genome coords
    :param cds_list: list of CDS
    :param keep_genes_set: a set of CDS names, if non None - the others will be filtered out
    :return: coord to CDS dictionary
    """
    pos_list.sort()
    res = {}
    for cds in cds_list:
        if keep_genes_set is not None and cds.name not in keep_genes_set:
            continue
        i = bisect_left(pos_list, cds.start)
        j = bisect_right(pos_list, cds.end, lo=i)
        for k in range(i, j):
            res[pos_list[k]] = cds
    return res


if __name__ == '__main__':
    cds_list = read_annotations(upstream_length, False)
    for cds in cds_list:
        if cds.is_hypothetical and cds.exists_in_proteom:
            print(cds.name + ' is hypothetical but exists in proteom')
        if not cds.is_hypothetical and not cds.exists_in_proteom and cds.type == CDSType.Gene:
            print(cds.name + ' is not hypothetical but does not exist in proteom')
