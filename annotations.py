from bisect import bisect_left, bisect_right
from enum import Enum

from constants import codon_table, complement, codon_table_compl


class CDSType(Enum):
    Gene = 0
    ncRNA = 1
    rRNA = 2
    tRNA = 3
    miscRNA = 4
    upstream = 5


class CDS:
    type = CDSType.Gene
    start = 0
    end = 0
    strand = 1  # 1 or -1
    name = ''
    synonym = ''


    def __init__(self, type, start, end, strand, name, synonym=''):
        self.name = name
        self.start = start
        self.end = end
        self.strand = strand
        self.type = type
        self.synonym = synonym


    def __lt__(self, other):
        return self.start < other.start


def read_annotations(path_to_annotations, upstream_length):
    genes = {}
    nc_rna = {}
    r_rna = {}
    t_rna = {}
    misc_rna = {}
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
            ns = s[8].split(' ')
            name = ns[1]
            if type == 'Gene':
                if len(ns) >= 5 and ns[3] == 'gene_synonym':
                    syn = ns[4]
                else:
                    syn = ''
                genes[name] = (start, end, strand, syn)
                if strand == 1:
                    upstream[name] = (start - upstream_length, start, strand, syn)
                else:
                    upstream[name] = (end, end + upstream_length, strand, syn)
            elif type == 'ncRNA':
                if ns[3] == 'gene_synonym':
                    syn = ns[4]
                else:
                    syn = ''
                nc_rna[name] = (start, end, strand, syn)
            elif type == 'rRNA':
                r_rna[name] = (start, end, strand)
            elif type == 'tRNA':
                t_rna[name] = (start, end, strand)
            elif type == 'Misc._RNA':
                misc_rna[name] = (start, end, strand)
    for rna_name in r_rna.keys():
        genes.pop(rna_name, None)
    for rna_name in nc_rna.keys():
        genes.pop(rna_name, None)
    for rna_name in t_rna.keys():
        genes.pop(rna_name, None)
    for rna_name in misc_rna.keys():
        genes.pop(rna_name, None)
    cds_list = []
    for name, value in r_rna.items():
        cds_list.append(CDS(CDSType.rRNA, value[0], value[1], value[2], name))
    for name, value in t_rna.items():
        cds_list.append(CDS(CDSType.tRNA, value[0], value[1], value[2], name))
    for name, value in misc_rna.items():
        cds_list.append(CDS(CDSType.miscRNA, value[0], value[1], value[2], name))
    for name, value in nc_rna.items():
        cds_list.append(CDS(CDSType.ncRNA, value[0], value[1], value[2], name, value[3]))
    for name, value in genes.items():
        cds_list.append(CDS(CDSType.Gene, value[0], value[1], value[2], name, value[3]))
    for name, value in upstream.items():
        cds_list.append(CDS(CDSType.upstream, value[0], value[1], value[2], name, value[3]))
    cds_list.sort()
    return cds_list


def localize_all_snps(all_snps, cds_list):
    res = {}
    for cds in cds_list:
        i = bisect_left(all_snps, cds.start)
        j = bisect_right(all_snps, cds.end, lo=i)
        for k in range(i, j):
            res[all_snps[k]] = cds
    return res


def fix_gene_annotations(cds_list, ref_seq):
    for cds in cds_list:
        if cds.type != CDSType.Gene:
            continue
        cds_len = cds.end - cds.start + 1
        if cds_len % 3 == 0:
            continue
        gene_seq = ref_seq[cds.start - 1: cds.end]
        if cds.strand == 1:
            i = cds_len // 3
            last_codon = gene_seq[3*(i - 1): 3*i]
            next_codon = ref_seq[cds.start - 1 + 3*i: cds.start - 1 + 3*i + 3]
            if codon_table[last_codon] == '*':
                cds.end -= cds_len % 3
                print('fixed cds ' + cds.name)
            elif codon_table[next_codon] == '*':
                cds.end += 3 - cds_len % 3
                print('fixed cds ' + cds.name)
            else:
                print('can\'t fix the cds ' + cds.name)
        else:
            gene_seq = gene_seq.translate(complement)[::-1]
            i = cds_len // 3
            last_codon = gene_seq[3*(i - 1): 3*i]
            next_codon = ref_seq[cds.start - 1 - 3 + cds_len % 3: cds.start - 1 - 3 + cds_len % 3 + 3]
            if codon_table[last_codon] == '*':
                cds.start += cds_len % 3
                print('fixed cds ' + cds.name)
            elif codon_table_compl[next_codon] == '*':
                cds.start -= 3 - cds_len % 3
                print('fixed cds ' + cds.name)
            else:
                print('can\'t fix the cds ' + cds.name)
