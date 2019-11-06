import multiprocessing as mp
from os import listdir, mkdir, makedirs
from os.path import exists

from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SubsMat.MatrixInfo import blosum62
from sklearn.externals.joblib import Parallel, delayed

from core.annotations import (CDSType, localize_all_variants,
                                  read_annotations)
from core.constants import (codon_table, codon_table_compl, complement,
                                data_path, upstream_length)
from core.data_reading import read_h37rv

path_to_variants = data_path + 'snps/gatk_before_cortex/raw_variants/'
out_path = data_path + 'snps/gatk_before_cortex/annotated_long_del_pg_NWds10/'
# path_to_snps = data_path + 'snps/skesa_mummer_raw_mum4/'
# out_path = data_path + 'snps/skesa_mum4_annotated_long_del_pg_NWds10/'
out_path_snp = out_path + 'all_var_list.csv'
debug_path = data_path + 'snps/gatk_before_cortex/dead_genes/'
# debug_path = data_path + 'snps/skesa_mum4_dead_genes/'

thread_num = 32


split_dels = False
translate_pseudogene = False
add_codons_to_frameshift = 0.1
shortening_threshold = 0.9


def read_variants(sample_id):
    snps = []
    with open(path_to_variants + sample_id + '.variants', 'r') as f1:
        for line in f1.readlines():
            if line[0] == '#':
                continue
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


def append_dels(res, deletions, ref_seq, ref_seq_compl):
    ref_seq_len = len(ref_seq)
    for cds, del_list in deletions.items():
        if cds == 'non_cds':
            for pos, l in del_list:
                coord = ref_seq_len - pos - l + 1
                res.append(
                    '\t'.join(('non_cds', '-', str(coord + 1), ref_seq_compl[coord: coord + l], '-', 'del')))
            continue
        if cds.strand == 1:
            if cds.type == CDSType.upstream:
                for pos, l in del_list:
                    res.append('\t'.join(
                        (cds.type.name, cds.name, str(pos - cds.end - 1), ref_seq[pos - 1: pos - 1 + l], '-',
                         'del')))
            else:
                for pos, l in del_list:
                    res.append('\t'.join(
                        (cds.type.name, cds.name, str(pos - cds.start + 1), ref_seq[pos - 1: pos - 1 + l], '-',
                         'del')))
        else:
            if cds.type == CDSType.upstream:
                for pos, l in del_list:
                    res.append('\t'.join(
                        (cds.type.name, cds.name, str(cds.start - pos - l),
                         ref_seq_compl[ref_seq_len - pos - l + 1: ref_seq_len - pos + 1], '-', 'del')))
            else:
                for pos, l in del_list:
                    res.append('\t'.join(
                        (cds.type.name, cds.name, str(cds.end - (pos + l - 2)),
                         ref_seq_compl[ref_seq_len - pos - l + 1: ref_seq_len - pos + 1], '-', 'del')))


def add_del(deletions, cds, pos):
    del_list = deletions.get(cds)
    if del_list is None:
        del_list = [[pos, 1]]
        deletions[cds] = del_list
    else:
        found = False
        for l in del_list:
            if pos == l[0] + l[1]:
                l[1] += 1
                found = True
                break
        if not found:
            del_list.append([pos, 1])


def process_non_cds(pos, alt, v_type, ref_seq_len, ref_seq_compl, deletions, res):
    coord = ref_seq_len - pos
    if v_type == 'snp':
        res.append('\t'.join(('non_cds', '-', str(coord + 1), ref_seq_compl[coord], complement[alt], v_type)))
    elif v_type == 'del':
        if split_dels:
            res.append(
                '\t'.join(('non_cds', '-', str(coord + 1), ref_seq_compl[coord], complement[alt], v_type)))
        else:
            add_del(deletions, 'non_cds', pos)
    else:
        # ins
        if len(alt) == 1:
            res.append(
                '\t'.join(('non_cds', '-', str(coord + 1), '-', complement[alt], v_type)))
        else:
            res.append(
                '\t'.join(('non_cds', '-', str(coord + 1), '-', str(Seq(alt).reverse_complement()), v_type)))


def process_non_gene(cds, pos, alt, v_type, ref_seq, ref_seq_compl, ref_seq_len, deletions, res):
    if cds.strand == 1:
        if cds.type == CDSType.upstream:
            if v_type == 'snp':
                res.append('\t'.join(
                    (cds.type.name, cds.name, str(pos - cds.end - 1), ref_seq[pos - 1], alt, v_type)))
            elif v_type == 'del':
                if split_dels:
                    res.append('\t'.join(
                        (cds.type.name, cds.name, str(pos - cds.end - 1), ref_seq[pos - 1], alt,
                         v_type)))
                else:
                    add_del(deletions, cds, pos)
            else:
                # ins
                res.append('\t'.join(
                    (cds.type.name, cds.name, str(pos - cds.end - 1), '-', alt, v_type)))
        else:
            if v_type == 'snp':
                res.append('\t'.join(
                    (cds.type.name, cds.name, str(pos - cds.start + 1), ref_seq[pos - 1], alt, v_type)))
            elif v_type == 'del':
                if split_dels:
                    res.append('\t'.join(
                        (cds.type.name, cds.name, str(pos - cds.start + 1), ref_seq[pos - 1], alt,
                         v_type)))
                else:
                    add_del(deletions, cds, pos)
            else:
                res.append('\t'.join(
                    (cds.type.name, cds.name, str(pos - cds.start + 1), '-', alt, v_type)))
    else:  # if strand is '-'
        if cds.type == CDSType.upstream:
            if v_type == 'snp':
                res.append('\t'.join(
                    (cds.type.name, cds.name, str(cds.start - pos - 1),
                     ref_seq_compl[ref_seq_len - pos], complement[alt], v_type)))
            elif v_type == 'del':
                if split_dels:
                    res.append('\t'.join(
                        (cds.type.name, cds.name, str(cds.start - pos - 1),
                         ref_seq_compl[ref_seq_len - pos], complement[alt], v_type)))
                else:
                    add_del(deletions, cds, pos)
            else:
                # ins
                if len(alt) == 1:
                    res.append('\t'.join(
                        (cds.type.name, cds.name, str(cds.start - pos - 1), '-', complement[alt], v_type)))
                else:
                    res.append('\t'.join(
                        (cds.type.name, cds.name, str(cds.start - pos - 1), '-',
                         str(Seq(alt).reverse_complement()), v_type)))
        else:
            if v_type == 'snp':
                res.append('\t'.join(
                    (cds.type.name, cds.name, str(cds.end - (pos - 1)), ref_seq_compl[ref_seq_len - pos],
                     complement[alt], v_type)))
            elif v_type == 'del':
                if split_dels:
                    res.append('\t'.join(
                        (cds.type.name, cds.name, str(cds.end - (pos - 1)),
                         ref_seq_compl[ref_seq_len - pos], complement[alt], v_type)))
                else:
                    add_del(deletions, cds, pos)
            else:
                if len(alt) == 1:
                    res.append('\t'.join(
                        (cds.type.name, cds.name, str(cds.end - (pos - 1)), '-', complement[alt], v_type)))
                else:
                    res.append('\t'.join(
                        (cds.type.name, cds.name, str(cds.end - (pos - 1)), '-', str(Seq(alt).reverse_complement()),
                         v_type)))


def process_gene_snp_only(gene, var_list, ref_seq, res, fout):
    snps = []
    j = 0
    is_broken = False
    if gene.strand == 1:
        while j < len(var_list):
            pos, alt, v_type = var_list[j]
            protein_pos = (pos - gene.start) // 3 + 1
            nucleotide_pos = (pos - gene.start) % 3
            # NEIGHBOURS ANALYSIS
            aa_ref, aa_alt, j = get_aminoacids_sense(ref_seq, nucleotide_pos, var_list, j)

            if aa_ref != aa_alt:
                if alt == '*':
                    # early stop codon
                    if pos - gene.start < shortening_threshold * (gene.end - gene.start + 1):
                        # new gene is too short
                        res.append('\t'.join((gene.type.name, gene.name, '-', '-', '-', 'FS')))
                        is_broken = True
                        fout.write(gene.name + '\n')
                        fout.write('reference in pos %d:\n' % protein_pos)
                        fout.write(aa_ref)
                        fout.write('\nmutated:\n')
                        fout.write(aa_alt)
                        fout.write('\n')
                        break
                snps.append('\t'.join(
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
                if alt == '*':
                    # early stop codon
                    if gene.end - pos < shortening_threshold * (gene.end - gene.start + 1):
                        # new gene is too short
                        res.append('\t'.join((gene.type.name, gene.name, '-', '-', '-', 'FS')))
                        is_broken = True
                        fout.write(gene.name + '\n')
                        fout.write('reference in pos %d:\n' % protein_pos)
                        fout.write(aa_ref)
                        fout.write('\nmutated:\n')
                        fout.write(aa_alt)
                        fout.write('\n')
                        break
                snps.append('\t'.join(
                    (gene.type.name, gene.name, str(protein_pos), aa_ref, aa_alt, v_type)))
            j += 1
    if not is_broken:
        res.extend(snps)


def process_gene_with_indels1(gene, var_list, ref_seq, res, fout):

    def align_proteins(translated_old, translated_new):
        # gap open/extension penalty combination for PAM120 is -16/-4, PAM250 -11/-1, BLOSUM50 -10/-2,
        # BLOSUM62 -7/-1, as recommended by Vingron and Waterman, Mount, and Pearson
        try:
            alignment = \
            pairwise2.align.globalds(str(translated_old)[:-1], str(translated_new)[:stop], blosum62, -7, -1)[0]
            aln_old = alignment[0]
            aln_new = alignment[1]
            insert_pos = 0
            insert_s = 0
            insert = False
            del_pos = 0
            del_s = 0
            deletion = False
            p = 0
            for i in range(len(aln_old)):
                c1 = aln_old[i]
                c2 = aln_new[i]
                if c1 != c2:
                    if c1 == '-':
                        # ins
                        if deletion and not split_dels:
                            res.append('\t'.join(
                                (gene.type.name, gene.name, str(del_pos + 1), aln_old[del_s: i], '-', 'del')))
                            deletion = False
                        if not insert:
                            insert_pos = p
                            insert_s = i
                            insert = True
                    elif c2 == '-':
                        # del
                        if insert:
                            res.append('\t'.join(
                                (gene.type.name, gene.name, str(insert_pos + 1), '-', aln_new[insert_s: i], 'ins')))
                            insert = False
                        if split_dels:
                            res.append('\t'.join(
                                (gene.type.name, gene.name, str(p + 1), c1, c2, 'del')))
                        else:
                            if not deletion:
                                del_pos = p
                                del_s = i
                                deletion = True
                        p += 1
                    else:
                        if insert:
                            res.append('\t'.join(
                                (gene.type.name, gene.name, str(insert_pos + 1), '-', aln_new[insert_s: i], 'ins')))
                            insert = False
                        res.append('\t'.join(
                            (gene.type.name, gene.name, str(p + 1), c1, c2, 'snp')))
                        p += 1
                else:
                    if insert:
                        res.append('\t'.join(
                            (gene.type.name, gene.name, str(insert_pos + 1), '-', aln_new[insert_s: i], 'ins')))
                        insert = False
                    if deletion and not split_dels:
                        res.append('\t'.join(
                            (gene.type.name, gene.name, str(del_pos + 1), aln_old[del_s: i], '-', 'del')))
                        deletion = False
                    p += 1
        except:
            print('error: stop=%d' % stop)

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
        # res.append('\t'.join(
        #     (gene.type.name, gene.name, '-', '-', '-', 'FS')))
        # continue
        # add some more letters to the end, might be new stop codon
        seq_new.append(ref_seq[gene.end: gene.end + 3 - new_len % 3 +
                                         3 * int(add_codons_to_frameshift * (gene.end - gene.start + 1))])
    else:
        seq_new.append(ref_seq[gene.end: gene.end +
                                         3 * int(add_codons_to_frameshift * (gene.end - gene.start + 1))])
    seq_old = Seq(ref_seq[gene.start - 1: gene.end])
    if len(seq_old) % 3 != 0:
        print(gene.name + ' wrong annotation')
        return
    seq_new = Seq(''.join(seq_new))
    # if new_len != len(seq_new):
    #     print('we have a problem!')
    if gene.strand == 1:
        translated_old = seq_old.translate()
        translated_new = seq_new.translate()
    else:
        translated_old = seq_old.reverse_complement().translate()
        translated_new = seq_new.reverse_complement().translate()
    # alignment = pairwise2.align.globalxx(str(translated_old), str(translated_new))[0]
    stop = -1
    for i in range(len(translated_new)):
        if translated_new[i] == '*':
            stop = i
            break
    if stop == -1 or stop == 0 or len(translated_new) < shortening_threshold * len(translated_old):
        fout.write(gene.name + '\n')
        fout.write('reference:\n')
        fout.write(str(translated_old))
        fout.write('\nmutated:\n')
        fout.write(str(translated_new))
        fout.write('\n')
        # gene is broken -> FS
        res.append('\t'.join((gene.type.name, gene.name, '-', '-', '-', 'FS')))
        return
    else:
        align_proteins(translated_old, translated_new)


def process_gene_with_indels(gene, var_list, ref_seq, res, fout):

    def align_proteins(translated_old, translated_new):
        # gap open/extension penalty combination for PAM120 is -16/-4, PAM250 -11/-1, BLOSUM50 -10/-2,
        # BLOSUM62 -7/-1, as recommended by Vingron and Waterman, Mount, and Pearson
        try:
            alignment = \
            pairwise2.align.globalds(str(translated_old)[:-1], str(translated_new), blosum62, -7, -1)[0]
            aln_old = alignment[0]
            aln_new = alignment[1]
            insert_pos = 0
            insert_s = 0
            insert = False
            del_pos = 0
            del_s = 0
            deletion = False
            p = 0
            for i in range(len(aln_old)):
                c1 = aln_old[i]
                c2 = aln_new[i]
                if c1 != c2:
                    if c1 == '-':
                        # ins
                        if deletion and not split_dels:
                            res.append('\t'.join(
                                (gene.type.name, gene.name, str(del_pos + 1), aln_old[del_s: i], '-', 'del')))
                            deletion = False
                        if not insert:
                            insert_pos = p
                            insert_s = i
                            insert = True
                    elif c2 == '-':
                        # del
                        if insert:
                            res.append('\t'.join(
                                (gene.type.name, gene.name, str(insert_pos + 1), '-', aln_new[insert_s: i], 'ins')))
                            insert = False
                        if split_dels:
                            res.append('\t'.join(
                                (gene.type.name, gene.name, str(p + 1), c1, c2, 'del')))
                        else:
                            if not deletion:
                                del_pos = p
                                del_s = i
                                deletion = True
                        p += 1
                    else:
                        if insert:
                            res.append('\t'.join(
                                (gene.type.name, gene.name, str(insert_pos + 1), '-', aln_new[insert_s: i], 'ins')))
                            insert = False
                        res.append('\t'.join(
                            (gene.type.name, gene.name, str(p + 1), c1, c2, 'snp')))
                        p += 1
                else:
                    if insert:
                        res.append('\t'.join(
                            (gene.type.name, gene.name, str(insert_pos + 1), '-', aln_new[insert_s: i], 'ins')))
                        insert = False
                    if deletion and not split_dels:
                        res.append('\t'.join(
                            (gene.type.name, gene.name, str(del_pos + 1), aln_old[del_s: i], '-', 'del')))
                        deletion = False
                    p += 1
        except:
            print('error: stop=%d' % stop)

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
        # res.append('\t'.join(
        #     (gene.type.name, gene.name, '-', '-', '-', 'FS')))
        # continue
        # add some more letters to the end, might be new stop codon
        seq_new.append(ref_seq[gene.end: gene.end + 3 - new_len % 3 +
                                         3 * int(add_codons_to_frameshift * (gene.end - gene.start + 1))])
    else:
        seq_new.append(ref_seq[gene.end: gene.end +
                                         3 * int(add_codons_to_frameshift * (gene.end - gene.start + 1))])
    seq_old = Seq(ref_seq[gene.start - 1: gene.end])
    if len(seq_old) % 3 != 0:
        print(gene.name + ' wrong annotation')
        return
    seq_new = Seq(''.join(seq_new))
    # if new_len != len(seq_new):
    #     print('we have a problem!')
    if gene.strand == 1:
        translated_old = seq_old.translate()
        translated_new = seq_new.translate()
    else:
        translated_old = seq_old.reverse_complement().translate()
        translated_new = seq_new.reverse_complement().translate()
    # alignment = pairwise2.align.globalxx(str(translated_old), str(translated_new))[0]
    start = -1
    stop = -1
    for i in range(len(translated_new)):
        if start == -1 and translated_new[i] == 'M':
            start = i
        if translated_new[i] == '*':
            stop = i
            break
    if start == -1 or stop == -1 or stop - start < shortening_threshold * len(translated_old):
        fout.write(gene.name + '\n')
        fout.write('reference:\n')
        fout.write(str(translated_old))
        fout.write('\nmutated:\n')
        fout.write(str(translated_new))
        fout.write('\n')
        # gene is broken -> FS
        res.append('\t'.join((gene.type.name, gene.name, '-', '-', '-', 'FS')))
        return
    else:
        align_proteins(translated_old, translated_new[start:stop])


def format_variants(sample_id, variants, ref_seq, ref_seq_compl, snp_to_cds):

    res = []
    ref_seq_len = len(ref_seq)
    genes_to_variants = {}
    deletions = {}
    with open(out_path + sample_id + '.variants', 'w') as f:
        for pos, alt, v_type in variants:
            try:
                cds = snp_to_cds[pos]
            except:
                cds = None

            if cds is not None:
                if cds.type == CDSType.Gene and (not cds.is_pseudogene or translate_pseudogene):
                    var_list = genes_to_variants.get(cds.name)
                    if var_list is None:
                        var_list = [(pos, alt, v_type)]
                        genes_to_variants[cds.name] = var_list
                    else:
                        var_list.append((pos, alt, v_type))
                else:
                    process_non_gene(cds, pos, alt, v_type, ref_seq, ref_seq_compl, ref_seq_len, deletions, res)
            else:
                process_non_cds(pos, alt, v_type, ref_seq_len, ref_seq_compl, deletions, res)

        if not split_dels:
            append_dels(res, deletions, ref_seq, ref_seq_compl)

        # working with protein coding genes
        with open(debug_path + sample_id + '.genes', 'w') as fout:
            for gene_name, var_list in genes_to_variants.items():
                snp_only = True
                for pos, alt, v_type in var_list:
                    if v_type != 'snp':
                        snp_only = False
                        break
                    else:
                        if len(alt) != 1:
                            print('wrong v_type!')
                gene = snp_to_cds[var_list[0][0]]
                if snp_only:
                    process_gene_snp_only(gene, var_list, ref_seq, res, fout)
                else:
                    process_gene_with_indels(gene, var_list, ref_seq, res, fout)
        if len(res) > 0:
            f.write('\n'.join(res))
            f.write('\n')
    return sample_id, res


def main():
    # mp.set_start_method('forkserver')
    if not exists(out_path):
        makedirs(out_path)
    if not exists(debug_path):
        makedirs(debug_path)
    sample_to_snps = {}
    all_snp_pos = set()
    h37rv = read_h37rv()
    h37rv_compl = str(Seq(h37rv).reverse_complement())

    sample_ids = [fname[0:-len('.variants')] for fname in listdir(path_to_variants) if fname.endswith('.variants')]

    cds_list = read_annotations(upstream_length, filter_by_gene_len=False)

    tasks = Parallel(n_jobs=thread_num)(
        delayed(read_variants)(sample_id) for sample_id in sample_ids)
    for sample_id, variants in tasks:
        sample_to_snps[sample_id] = variants
        for snp_pos, alt, v_type in variants:
            all_snp_pos.add(snp_pos)
    all_snps = list(all_snp_pos)
    print('done with variant reading')

    snp_to_cds = localize_all_variants(all_snps, cds_list)
    print('done with variant localization')

    tasks = Parallel(n_jobs=thread_num, batch_size=len(sample_ids)//thread_num + 1)(
        delayed(format_variants)(sample_id, variants, h37rv, h37rv_compl, snp_to_cds)
        for sample_id, variants in sample_to_snps.items())
    formatted_snps = {sample_id: variants for sample_id, variants in tasks}
    print('done with variant format')

    all_formatted_snps = set()
    for f_snps in formatted_snps.values():
        all_formatted_snps.update(f_snps)
    all_formatted_snps = list(all_formatted_snps)
    all_formatted_snps.sort()

    with open(out_path_snp, 'w') as f:
        f.write('\n'.join(all_formatted_snps))
        f.write('\n')
    print('printed all snps')


if __name__ == '__main__':
    main()
