from bisect import bisect_left, bisect_right
from enum import Enum


class CDSType(Enum):
    Gene = 0
    ncRNA = 1
    rRNA = 2
    tRNA = 3
    upstream = 4


class CDS:
    type = CDSType.Gene
    start = 0
    end = 0
    strand = 1
    name = ''


    def __init__(self, type, start, end, strand, name):
        self.name = name
        self.start = start
        self.end = end
        self.strand = strand
        self.type = type


    def __lt__(self, other):
        return self.start < other.start


def read_annotations(path_to_annotations, upstream_length):
    genes = {}
    nc_rna = {}
    r_rna = {}
    t_rna = {}
    upstream = {}

    with open(path_to_annotations, 'r') as f:
        for line in f.readlines():
            s = line[:-1].split('\t')
            type = s[2]
            start = int(s[3])
            end = int(s[4])
            if type in ('Repeat_region', 'mobile_element', 'CDS'):
                continue
            strand = 1 if s[6] == '+' else -1
            name = s[8].split(' ')[1]
            if type == 'Gene':
                genes[name] = (start, end, strand)
                if strand == 1:
                    upstream[name] = (start - upstream_length, start, strand)
                else:
                    upstream[name] = (end, end + upstream_length, strand)
            elif type == 'ncRNA':
                nc_rna[name] = (start, end, strand)
            elif type == 'rRNA':
                r_rna[name] = (start, end, strand)
            elif type == 'tRNA':
                t_rna[name] = (start, end, strand)
    for rna_name in r_rna.keys():
        genes.pop(rna_name, None)
    for rna_name in nc_rna.keys():
        genes.pop(rna_name, None)
    for rna_name in t_rna.keys():
        genes.pop(rna_name, None)
    cds = []
    for name, value in r_rna.items():
        cds.append(CDS(CDSType.rRNA, value[0], value[1], value[2], name))
    for name, value in t_rna.items():
        cds.append(CDS(CDSType.tRNA, value[0], value[1], value[2], name))
    for name, value in nc_rna.items():
        cds.append(CDS(CDSType.ncRNA, value[0], value[1], value[2], name))
    for name, value in genes.items():
        cds.append(CDS(CDSType.Gene, value[0], value[1], value[2], name))
    for name, value in upstream.items():
        cds.append(CDS(CDSType.upstream, value[0], value[1], value[2], name))
    cds.sort()
    return cds


def localize_all_snps(all_snps, cds_list):
    res = {}
    for cds in cds_list:
        i = bisect_left(all_snps, cds.start)
        j = bisect_right(all_snps, cds.end, lo=i)
        for k in range(i, j):
            res[all_snps[k]] = cds
    return res
