from Bio.Seq import Seq

from src.core.annotations import read_annotations, CDSType
from src.core.constants import upstream_length
from src.core.data_reading import read_h37rv, read_variants

gene_names = ['Rv3216', 'oxyR\'']
sample_ids = ['SAMEA1707255', 'SAMN02469331']

path_to_snps = '../../data/snps/raw_with_DR_with_indel_with_pheno_and_snp_mc10/'


def print_ref_seq_genes_seq():
    ref_seq = read_h37rv()
    cds_list = read_annotations(upstream_length, False)

    for cds in cds_list:
        if cds.type == CDSType.Gene and cds.name in gene_names:
            print(cds.name)
            seq = Seq(ref_seq[cds.start - 1: cds.end])
            if cds.strand == 1:
                print(seq.reverse_complement())
            else:
                print(seq)


def print_sample_its():
    ref_seq = read_h37rv()
    cds_list = read_annotations(upstream_length, True)
    rRNA_genes = []
    for cds in cds_list:
        if cds.type == CDSType.rRNA:
            rRNA_genes.append(cds)
    rRNA_genes.sort()
    prev = rRNA_genes[0]
    its_list = {}
    for i in range(1, len(rRNA_genes)):
        curr = rRNA_genes[i]
        its_list['ITS' + str(i)] = (prev.end + 1, curr.start - 1)
        prev = curr
    rRNA_name_to_gene = {}
    for gene in rRNA_genes:
        rRNA_name_to_gene[gene.name] = (gene.start, gene.end)
    for sample_id in sample_ids:
        print(sample_id)
        sample_id, variants = read_variants(path_to_snps, sample_id)
        its_to_var = {}
        for name in its_list.keys():
            its_to_var[name] = []
        rna_gene_to_var = {}
        for gene in rRNA_genes:
            rna_gene_to_var[gene.name] = []
        for var in variants:
            s = var.split('\t')
            pos = int(s[0])
            for name, (start, end) in its_list.items():
                if start <= pos <= end:
                    its_to_var[name].append(s)
            for gene in rRNA_genes:
                if gene.start <= pos <= gene.end:
                    rna_gene_to_var[gene.name].append(s)
        for name, var_list in its_to_var.items():
            start, end = its_list[name]
            l = []
            last_pos = start
            for pos, alt, type in var_list:
                pos = int(pos)
                l.append(ref_seq[last_pos:pos - 1])
                if type == 'snp':
                    l.append(alt)
                    last_pos = pos
                elif type == 'ins':
                    l.append(ref_seq[pos - 1])
                    l.append(alt)
                    last_pos = pos
                else:
                    last_pos = pos
            l.append(ref_seq[last_pos:end])
            print(name)
            print(Seq(''.join(l)).reverse_complement())
        for name, var_list in rna_gene_to_var.items():
            start, end = rRNA_name_to_gene[name]
            l = []
            last_pos = start
            for pos, alt, type in var_list:
                pos = int(pos)
                l.append(ref_seq[last_pos:pos - 1])
                if type == 'snp':
                    l.append(alt)
                    last_pos = pos
                elif type == 'ins':
                    l.append(ref_seq[pos - 1])
                    l.append(alt)
                    last_pos = pos
                else:
                    last_pos = pos
            l.append(ref_seq[last_pos:end])
            print(name)
            print(Seq(''.join(l)).reverse_complement())


if __name__ == '__main__':
    print_sample_its()
