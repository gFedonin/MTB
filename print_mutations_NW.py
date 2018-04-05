from os import mkdir
from os.path import exists
from Bio import SeqIO
from Bio.Seq import Seq
from sklearn.externals.joblib import Parallel, delayed
from Bio import pairwise2

from annotations import CDSType, read_annotations, localize_all_snps, fix_gene_annotations
from constants import codon_table, codon_table_compl, complement

path_to_ids = './data/all_with_pheno_and_snp.txt'
path_to_snps = './data/snps/raw_with_DR_with_indel_with_pheno_and_snp_mc10/'
out_path = './data/snps/annotated_with_DR_with_indel_with_pheno_and_snp_mc10/'
out_path_snp = out_path + 'all_var_list.csv'
path_to_annotations = './data/AL123456_rev.gff'
path_to_ref = './data/h37rv.fasta'

thread_num = 144
upstream_length = 100


def read_h37rv():
    fasta_sequences = SeqIO.parse(open(path_to_ref), 'fasta')
    sequence = ''
    for fasta in fasta_sequences:
        sequence = str(fasta.seq)
    return sequence.upper()


def read_variants(sample_id):
    snps = []
    with open(path_to_snps + sample_id + '.variants', 'r') as f1:
        lines = f1.readlines()
        for line in lines:
            s = line.strip().split('\t')
            snps.append((int(s[0]), s[1], s[2]))
    return sample_id, snps


def get_aminoacids_sense(ref_seq, nucleotide_pos, snps, i):
    var_num = len(snps)
    pos, alt, v_type = snps[i]
    if nucleotide_pos == 0:
        aa0 = codon_table[ref_seq[pos - 1:pos + 2]]
        if i != var_num - 1:
            pos1, alt1, v_type1 = snps[i + 1]
            if pos1 == pos + 1:
                if i != var_num - 2 and snps[i + 2][0] == pos + 2:
                    aa1 = codon_table[alt + alt1 + snps[i + 2][1]]
                    return aa0, aa1, i + 2
                else:
                    aa1 = codon_table[alt + alt1 + ref_seq[pos + 1]]
                    return aa0, aa1, i + 1
            elif pos1 == pos + 2:
                aa1 = codon_table[alt + ref_seq[pos] + alt1]
                return aa0, aa1, i + 1
        aa1 = codon_table[alt + ref_seq[pos] + ref_seq[pos + 1]]
        return aa0, aa1, i
    elif nucleotide_pos == 1:
        aa0 = codon_table[ref_seq[pos - 2:pos + 1]]
        if i != var_num - 1 and snps[i + 1][0] == pos + 1:
            aa1 = codon_table[ref_seq[pos - 2] + alt + snps[i + 1][1]]
            return aa0, aa1, i + 1
        else:
            aa1 = codon_table[ref_seq[pos - 2] + alt + ref_seq[pos]]
            return aa0, aa1, i
    else:
        aa0 = codon_table[ref_seq[pos - 3:pos]]
        aa1 = codon_table[ref_seq[pos - 3:pos - 1] + alt]
        return aa0, aa1, i


def get_aminoacids_antisense(ref_seq, nucleotide_pos, snps, i):
    snp_num = len(snps)
    pos, alt, v_type = snps[i]
    if nucleotide_pos == 2:
        aa0 = codon_table_compl[ref_seq[pos - 1:pos + 2]]
        if i != snp_num - 1:
            pos1, alt1, v_type1 = snps[i + 1]
            if pos1 == pos + 1:
                if i != snp_num - 2 and snps[i + 2][0] == pos + 2:
                    aa1 = codon_table_compl[alt + alt1 + snps[i + 2][1]]
                    return aa0, aa1, i + 2
                else:
                    aa1 = codon_table_compl[alt + alt1 + ref_seq[pos + 1]]
                    return aa0, aa1, i + 1
            elif pos1 == pos + 2:
                aa1 = codon_table_compl[alt + ref_seq[pos] + alt1]
                return aa0, aa1, i + 1
        aa1 = codon_table_compl[alt + ref_seq[pos:pos + 2]]
        return aa0, aa1, i
    elif nucleotide_pos == 1:
        aa0 = codon_table_compl[ref_seq[pos - 2:pos + 1]]
        if i != snp_num - 1 and snps[i + 1][0] == pos + 1:
            aa1 = codon_table_compl[ref_seq[pos - 2] + alt + snps[i + 1][1]]
            return aa0, aa1, i + 1
        aa1 = codon_table_compl[ref_seq[pos - 2] + alt + ref_seq[pos]]
        return aa0, aa1, i
    else:
        aa0 = codon_table_compl[ref_seq[pos - 3:pos]]
        aa1 = codon_table_compl[ref_seq[pos - 3:pos - 1] + alt]
        return aa0, aa1, i


def format_variant(sample_id, variants, ref_seq, ref_seq_compl, snp_to_cds):

    res = []
    ref_seq_len = len(ref_seq)
    genes_to_variants = {}
    with open(out_path + sample_id + '.variants', 'w') as f:
        for pos, alt, v_type in variants:
            try:
                cds = snp_to_cds[pos]
            except:
                cds = None

            if cds is not None:
                if cds.type == CDSType.Gene:
                    var_list = genes_to_variants.get(cds.name)
                    if var_list is None:
                        var_list = [(pos, alt, v_type)]
                        genes_to_variants[cds.name] = var_list
                    else:
                        var_list.append((pos, alt, v_type))
                else:
                    if cds.strand == 1:
                        if cds.type == CDSType.upstream:
                            if v_type == 'snp' or v_type == 'del':
                                res.append('\t'.join(
                                    (cds.type.name, cds.name, str(pos - cds.end - 1), ref_seq[pos - 1], alt, v_type)))
                            else:
                                # ins
                                res.append('\t'.join(
                                    (cds.type.name, cds.name, str(pos - cds.end - 1), '-', alt, v_type)))
                        else:
                            if v_type == 'snp' or v_type == 'del':
                                res.append('\t'.join(
                                    (cds.type.name, cds.name, str(pos - cds.start + 1), ref_seq[pos - 1], alt, v_type)))
                            else:
                                res.append('\t'.join(
                                    (cds.type.name, cds.name, str(pos - cds.start + 1), '-', alt, v_type)))
                    else:  # if strand is '-'
                        if cds.type == CDSType.upstream:
                            if v_type == 'snp' or v_type == 'del':
                                res.append('\t'.join(
                                    (cds.type.name, cds.name, str(cds.start - pos - 1), ref_seq_compl[pos - 1], complement[alt], v_type)))
                            else:
                                # ins
                                if len(alt) == 1:
                                    res.append('\t'.join(
                                        (cds.type.name, cds.name, str(cds.start - pos - 1), '-', complement[alt], v_type)))
                                else:
                                    res.append('\t'.join(
                                        (cds.type.name, cds.name, str(cds.start - pos - 1), '-', str(Seq(alt).reverse_complement()), v_type)))
                        else:
                            if v_type == 'snp' or v_type == 'del':
                                res.append('\t'.join(
                                    (cds.type.name, cds.name, str(cds.end - (pos - 1)), ref_seq_compl[ref_seq_len - pos], complement[alt], v_type)))
                            else:
                                if len(alt) == 1:
                                    res.append('\t'.join(
                                        (cds.type.name, cds.name, str(cds.end - (pos - 1)), '-', complement[alt], v_type)))
                                else:
                                    res.append('\t'.join(
                                        (cds.type.name, cds.name, str(cds.end - (pos - 1)), '-', str(Seq(alt).reverse_complement()),
                                         v_type)))
            else:
                coord = ref_seq_len - pos
                if v_type == 'snp' or v_type == 'del':
                    res.append('\t'.join(('non_cds', ' ', str(coord + 1), ref_seq_compl[coord], complement[alt], v_type)))
                else:
                    # ins
                    if len(alt) == 1:
                        res.append(
                            '\t'.join(('non_cds', ' ', str(coord + 1), '-', complement[alt], v_type)))
                    else:
                        res.append(
                            '\t'.join(('non_cds', ' ', str(coord + 1), '-', str(Seq(alt).reverse_complement()), v_type)))
        processed_genes = set()
        for gene_name, var_list in genes_to_variants.items():
            snp_only = True
            for pos, alt, v_type in var_list:
                if v_type != 'snp':
                    snp_only = False
                    break
                else:
                    if len(alt) != 1:
                        print('wrong v_type!')
            if snp_only:
                processed_genes.add(gene_name)
                j = 0
                gene = snp_to_cds[var_list[0][0]]
                if gene.strand == 1:
                    while j < len(var_list):
                        pos, alt, v_type = var_list[j]
                        protein_pos = (pos - gene.start) // 3 + 1
                        nucleotide_pos = (pos - gene.start) % 3
                        # NEIGHBOURS ANALYSIS
                        aa_ref, aa_alt, j = get_aminoacids_sense(ref_seq, nucleotide_pos, var_list, j)

                        if aa_ref != aa_alt:
                            res.append('\t'.join(
                                (gene.type.name, gene.name, str(protein_pos), aa_ref, aa_alt, v_type)))
                        j += 1
                else:
                    while j < len(var_list):
                        pos, alt, v_type = var_list[j]
                        protein_pos = (gene.end - pos) // 3 + 1
                        nucleotide_pos = (gene.end - pos) % 3
                        # NEIGHBOURS ANALYSIS
                        aa_ref, aa_alt, j = get_aminoacids_antisense(ref_seq, nucleotide_pos, var_list, j)

                        if aa_ref != aa_alt:
                            res.append('\t'.join(
                                (gene.type.name, gene.name, str(protein_pos), aa_ref, aa_alt, v_type)))
                        j += 1
        for gene_name, var_list in genes_to_variants.items():
            if gene_name in processed_genes:
                continue
            gene = snp_to_cds[var_list[0][0]]
            seq_new = []
            new_len = 0
            s = gene.start - 1
            for pos, alt, v_type in var_list:
                if v_type == 'snp':
                    seq_new.append(ref_seq[s: pos - 1])
                    new_len += pos - s
                    seq_new.append(alt)
                    s = pos
                elif v_type == 'ins':
                    seq_new.append(ref_seq[s: pos])
                    seq_new.append(alt)
                    new_len += pos - s + len(alt)
                    s = pos
                else:
                    # del
                    if s != pos - 1:
                        seq_new.append(ref_seq[s: pos - 1])
                        new_len += pos - 1 - s
                    s = pos
            seq_new.append(ref_seq[s:gene.end])
            new_len += gene.end - s
            if new_len % 3 != 0:
                res.append('\t'.join(
                    (gene.type.name, gene.name, '-', '-', '-', 'FS')))
                continue
            seq_old = Seq(ref_seq[gene.start - 1: gene.end])
            if len(seq_old) % 3 != 0:
                print(gene_name + ' wrong annotation')
                continue
            seq_new = Seq(''.join(seq_new))
            if new_len != len(seq_new):
                print('we have a problem!')
            if gene.strand == 1:
                translated_old = seq_old.translate()
                translated_new = seq_new.translate()
            else:
                translated_old = seq_old.reverse_complement().translate()
                translated_new = seq_new.reverse_complement().translate()
            # alignment = pairwise2.align.globaldx(str(translated_old), str(translated_new), matlist.blosum62)[0]
            alignment = pairwise2.align.globalxx(str(translated_old), str(translated_new))[0]
            aln_old = alignment[0]
            aln_new = alignment[1]
            insert_pos = 0
            insert_s = 0
            insert = False
            p = 0
            for i in range(len(aln_old)):
                c1 = aln_old[i]
                c2 = aln_new[i]
                if c1 != c2:
                    if c1 == '-':
                        # ins
                        if not insert:
                            insert_pos = p
                            insert_s = i
                    elif c2 == '-':
                        #del
                        if insert:
                            res.append('\t'.join(
                                (gene.type.name, gene.name, str(insert_pos), '-', aln_new[insert_s: i], 'ins')))
                            insert = False
                        res.append('\t'.join(
                            (gene.type.name, gene.name, str(p), c1, c2, 'del')))
                        p += 1
                    else:
                        if insert:
                            res.append('\t'.join(
                                (gene.type.name, gene.name, str(insert_pos), '-', aln_new[insert_s: i], 'ins')))
                            insert = False
                        res.append('\t'.join(
                            (gene.type.name, gene.name, str(p), c1, c2, 'snp')))
                        p += 1
                else:
                    if insert:
                        res.append('\t'.join(
                            (gene.type.name, gene.name, str(insert_pos), '-', aln_new[insert_s: i], 'ins')))
                        insert = False
                    p += 1
        if len(res) > 0:
            f.write('\n'.join(res))
            f.write('\n')
    return res


def format_sample_set(sample_ids, sample_to_snps, ref_seq, ref_seq_compl, snp_to_cds):
    res = {}
    for sample_id in sample_ids:
        res[sample_id] = format_variant(sample_id, sample_to_snps[sample_id], ref_seq, ref_seq_compl, snp_to_cds)
    return res


def format_all_samples(sample_ids, sample_to_snps, ref_seq, ref_seq_compl, snp_to_cds):
    sample_lists = {}
    for i in range(thread_num):
        sample_lists[i] = []
    for i in range(len(sample_ids)):
        sample_lists[i%thread_num].append(sample_ids[i])
    tasks = Parallel(n_jobs=thread_num)(
        delayed(format_sample_set)(sample_list, sample_to_snps, ref_seq, ref_seq_compl, snp_to_cds)
        for sample_list in sample_lists.values()
    )
    res = {}
    for task in tasks:
        for sample_id, var_list in task.items():
            res[sample_id] = var_list
    return res


def main():
    if not exists(out_path):
        mkdir(out_path)
    sample_to_snps = {}
    all_snp_pos = set()
    h37rv = read_h37rv()
    h37rv_compl = str(Seq(h37rv).reverse_complement())

    sample_ids = [sample_id[:-1] for sample_id in open(path_to_ids, 'r').readlines()]

    cds_list = read_annotations(path_to_annotations, upstream_length)
    fix_gene_annotations(cds_list, h37rv)

    tasks = Parallel(n_jobs=-1)(
        delayed(read_variants)(sample_id) for sample_id in sample_ids)
    for sample_id, variants in tasks:
        sample_to_snps[sample_id] = variants
        for snp_pos, alt, v_type in variants:
            all_snp_pos.add(snp_pos)
    all_snps = list(all_snp_pos)
    all_snps.sort()
    print('done with variant reading')

    snp_to_cds = localize_all_snps(all_snps, cds_list)

    formatted_snps = format_all_samples(sample_ids, sample_to_snps, h37rv, h37rv_compl, snp_to_cds)
    print('done with variant format')

    all_formatted_snps = set()
    for f_snps in formatted_snps:
        for formatted_snp in f_snps:
            all_formatted_snps.add(formatted_snp)

    with open(out_path_snp, 'w') as f:
        f.write('\n'.join(all_formatted_snps))
        f.write('\n')
    print('printed all snps')


if __name__ == '__main__':
    main()
