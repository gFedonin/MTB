from subprocess import check_call

import pysam
from os import makedirs
from os.path import exists

from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SubsMat.MatrixInfo import blosum62
from sklearn.externals.joblib import Parallel, delayed

from core.annotations import (CDSType,
                                  read_annotations)
from core.constants import (codon_table, codon_table_compl, complement,
                                data_path, upstream_length)
from core.data_reading import read_h37rv, path_to_ref

# path_to_ids = data_path + 'all_with_pheno.txt'
path_to_ids = data_path + 'snps/combined_raw_variants_mq40_keep_complex_std_names_filtered/samples_filtered.list'
# path_to_ids = data_path + 'snps/combined_raw_variants_mq40_keep_complex_std_names_filtered/combined5.list'
# path_to_ids = data_path + 'debug1.list'
# path_to_ids = data_path + 'list_54.txt'
# path_to_variants = data_path + 'snps/pilon/raw_variants_filtered/'
# path_to_variants = data_path + 'snps/gatk_before_cortex/raw_variants_mq40_keep_complex_filtered/'
# path_to_variants = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_no_win_qual_mqm_std3_mqm30_filter_samples_first/'
path_to_variants = data_path + 'snps/combined_raw_variants_mq40_keep_complex_std_names_filtered/'
# out_path = data_path + 'snps/pilon/annotated_pg_NWds10_filtered_test/'
# out_path = data_path + 'snps/gatk_before_cortex/annotated_pg_NWds10_mq40_keep_complex_filtered_test/'
# out_path = data_path + 'snps/annotated_pg_NWds10_no_win_qual_mqm_std3_mqm30_filter_samples_first/'
out_path = data_path + 'snps/annotated_fixed_no_rep_long_del_pg_NWds10_combined_new/'
bam_path = out_path + 'bam/'
fasta_path = out_path + 'fasta/'
path_to_minimap = '/export/home/fedonin/minimap2-2.17_x64-linux/minimap2'
# path_to_snps = data_path + 'snps/skesa_mummer_raw_mum4/'
# out_path = data_path + 'snps/skesa_mum4_annotated_long_del_pg_NWds10/'
# out_path_snp = out_path + 'all_var_list.csv'
# broken_genes_path = data_path + 'snps/pilon/dead_genes_filtered_test/'
# broken_genes_path = data_path + 'snps/gatk_before_cortex/dead_genes_keep_complex_filtered_test/'
# broken_genes_path = data_path + 'snps/dead_genes_no_win_qual_mqm_std3_mqm30_filter_samples_first/'
broken_genes_path = data_path + 'snps/dead_genes_combined/'
# debug_path = data_path + 'snps/skesa_mum4_dead_genes/'

thread_num = 144


split_dels = False
translate_pseudogene = False
# add_codons_to_frameshift = 0.1
shortening_threshold = 0.9
# max_del_len = 1000
overwrite = False


def get_aminoacids_sense(ref_seq, nucleotide_pos, variants, i):
    var_num = len(variants)
    variant = variants[i]
    if nucleotide_pos == 0:
        aa0 = codon_table[ref_seq[variant.pos - 1:variant.pos + 2]]
        if i != var_num - 1:
            variant1 = variants[i + 1]
            if variant1.pos == variant.pos + 1:
                if i != var_num - 2 and variants[i + 2].pos == variant.pos + 2:
                    aa1 = codon_table[variant.alt + variant1.alt + variants[i + 2].ref]
                    return aa0, aa1, i + 2
                else:
                    aa1 = codon_table[variant.alt + variant1.alt + ref_seq[variant.pos + 1]]
                    return aa0, aa1, i + 1
            elif variant1.pos == variant.pos + 2:
                aa1 = codon_table[variant.alt + ref_seq[variant.pos] + variant1.alt]
                return aa0, aa1, i + 1
        aa1 = codon_table[variant.alt + ref_seq[variant.pos] + ref_seq[variant.pos + 1]]
        return aa0, aa1, i
    elif nucleotide_pos == 1:
        aa0 = codon_table[ref_seq[variant.pos - 2:variant.pos + 1]]
        if i != var_num - 1 and variants[i + 1].pos == variant.pos + 1:
            aa1 = codon_table[ref_seq[variant.pos - 2] + variant.alt + variants[i + 1].ref]
            return aa0, aa1, i + 1
        else:
            aa1 = codon_table[ref_seq[variant.pos - 2] + variant.alt + ref_seq[variant.pos]]
            return aa0, aa1, i
    else:
        aa0 = codon_table[ref_seq[variant.pos - 3:variant.pos]]
        aa1 = codon_table[ref_seq[variant.pos - 3:variant.pos - 1] + variant.alt]
        return aa0, aa1, i


def get_aminoacids_antisense(ref_seq, nucleotide_pos, variants, i):
    snp_num = len(variants)
    variant = variants[i]
    if nucleotide_pos == 2:
        aa0 = codon_table_compl[ref_seq[variant.pos - 1:variant.pos + 2]]
        if i != snp_num - 1:
            variant1 = variants[i + 1]
            if variant1.pos == variant.pos + 1:
                if i != snp_num - 2 and variants[i + 2].pos == variant.pos + 2:
                    aa1 = codon_table_compl[variant.alt + variant1.alt + variants[i + 2].ref]
                    return aa0, aa1, i + 2
                else:
                    aa1 = codon_table_compl[variant.alt + variant1.alt + ref_seq[variant.pos + 1]]
                    return aa0, aa1, i + 1
            elif variant1.pos == variant.pos + 2:
                aa1 = codon_table_compl[variant.alt + ref_seq[variant.pos] + variant1.alt]
                return aa0, aa1, i + 1
        aa1 = codon_table_compl[variant.alt + ref_seq[variant.pos:variant.pos + 2]]
        return aa0, aa1, i
    elif nucleotide_pos == 1:
        aa0 = codon_table_compl[ref_seq[variant.pos - 2:variant.pos + 1]]
        if i != snp_num - 1 and variants[i + 1].pos == variant.pos + 1:
            aa1 = codon_table_compl[ref_seq[variant.pos - 2] + variant.alt + variants[i + 1].ref]
            return aa0, aa1, i + 1
        aa1 = codon_table_compl[ref_seq[variant.pos - 2] + variant.alt + ref_seq[variant.pos]]
        return aa0, aa1, i
    else:
        aa0 = codon_table_compl[ref_seq[variant.pos - 3:variant.pos]]
        aa1 = codon_table_compl[ref_seq[variant.pos - 3:variant.pos - 1] + variant.alt]
        return aa0, aa1, i


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
            variant = var_list[j]
            protein_pos = (variant.pos - gene.start) // 3 + 1
            nucleotide_pos = (variant.pos - gene.start) % 3
            # NEIGHBOURS ANALYSIS
            aa_ref, aa_alt, j = get_aminoacids_sense(ref_seq, nucleotide_pos, var_list, j)

            if aa_ref != aa_alt:
                if aa_alt == '*':
                    # early stop codon
                    if variant.pos - gene.start < shortening_threshold * (gene.end - gene.start + 1):
                        # new gene is too short
                        res.append('\t'.join((gene.type.name, gene.name, '-', '-', '-', 'FS')))
                        is_broken = True
                        fout.write(gene.name + '\n')
                        # fout.write('reference in pos %d:\n' % protein_pos)
                        # fout.write(aa_ref)
                        # fout.write('\nmutated:\n')
                        # fout.write(aa_alt)
                        # fout.write('\n')
                        break
                snps.append('\t'.join(
                    (gene.type.name, gene.name, str(protein_pos), aa_ref, aa_alt, 'snp')))
            j += 1
    else:
        while j < len(var_list):
            var = var_list[j]
            protein_pos = (gene.end - var.pos) // 3 + 1
            nucleotide_pos = (gene.end - var.pos) % 3
            # NEIGHBOURS ANALYSIS
            aa_ref, aa_alt, j = get_aminoacids_antisense(ref_seq, nucleotide_pos, var_list, j)

            if aa_ref != aa_alt:
                if aa_alt == '*':
                    # early stop codon
                    if gene.end - var.pos < shortening_threshold * (gene.end - gene.start + 1):
                        # new gene is too short
                        res.append('\t'.join((gene.type.name, gene.name, '-', '-', '-', 'FS')))
                        is_broken = True
                        fout.write(gene.name + '\n')
                        # fout.write('reference in pos %d:\n' % protein_pos)
                        # fout.write(aa_ref)
                        # fout.write('\nmutated:\n')
                        # fout.write(aa_alt)
                        # fout.write('\n')
                        break
                snps.append('\t'.join(
                    (gene.type.name, gene.name, str(protein_pos), aa_ref, aa_alt, 'snp')))
            j += 1
    if not is_broken:
        res.extend(snps)


class Variant:

    def __init__(self, pos, ref, alt):
        self.pos = pos
        self.alt = alt
        self.ref = ref
        self.cds_list = []

    def __eq__(self, other):
        return self.pos == other.pos and len(self.ref) == len(self.alt)

    def __lt__(self, other):
        if self.pos == other.pos:
            if len(self.ref) < len(self.alt):
                return True
            else:
                return False
        return self.pos < other.pos


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
                                (gene.type.name, gene.name, str(p + 1), c1, c2, 'snp')))
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
    s = gene.start - 1
    # lens = []
    var_list.sort()
    for var in var_list:
        if var.pos >= gene.start:
            seq_new.append(ref_seq[s: var.pos - 1])
            seq_new.append(var.alt)
            s = var.pos + len(var.ref) - 1
        else:
            if var.pos + len(var.ref) >= gene.end and len(var.alt) - (gene.start - var.pos) < shortening_threshold * (gene.end - gene.start + 1):
                res.append('\t'.join((gene.type.name, gene.name, '-', '-', '-', 'FS')))
                return
            ref_shift = len(var.ref)%3
            alt_shift = len(var.alt)%3
            frame_shift_compensation = 3 - abs(ref_shift - alt_shift)
            seq_new.append(ref_seq[var.pos - frame_shift_compensation - 1: var.pos - 1])
            seq_new.append(var.alt)
            s = var.pos + len(var.ref) - 1
    if s < gene.end:
        seq_new.append(ref_seq[s:gene.end])
    seq_old = Seq(ref_seq[gene.start - 1: gene.end])
    seq_new_str = ''.join(seq_new)
    new_len = len(seq_new_str)
    if new_len % 3 != 0:
        if gene.strand == 1:
            seq_new.append(ref_seq[gene.end: gene.end + 3 - new_len % 3])
        else:
            seq_new.insert(0, ref_seq[gene.start - (3 - new_len % 3): gene.start])
    if len(seq_old) % 3 != 0:
        print(gene.name + ' wrong annotation')
        return
    seq_new = Seq(''.join(seq_new))
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
        if start == -1 and (translated_new[i] == translated_old[0] or translated_new[i] == 'M'):
            start = i
        if translated_new[i] == '*':
            stop = i
            break
    if start == -1 or stop == -1 or stop - start < shortening_threshold * len(translated_old):
        fout.write(gene.name + '\n')
        # fout.write('reference:\n')
        # fout.write(str(translated_old))
        # fout.write('\nmutated_extended:\n')
        # fout.write(str(translated_new))
        # fout.write('\nmutated_start-end:\n')
        # if start == -1:
        #     fout.write('no start')
        # elif stop == -1:
        #     fout.write('no stop')
        # else:
        #     fout.write(str(translated_new[start:stop]))
        # fout.write('\n')
        # gene is broken -> FS
        res.append('\t'.join((gene.type.name, gene.name, '-', '-', '-', 'FS')))
        return
    else:
        align_proteins(translated_old, translated_new[start:stop])


def gen_raw_variants(sample_id, ref_seq):
    if not exists(path_to_variants + sample_id + '.variants'):
        return sample_id, []
    if not overwrite and exists(out_path + sample_id + '.variants'):
        return sample_id, []
    mutated_genome = []
    pos = 0
    with open(path_to_variants + sample_id + '.variants', 'r') as f1:
        for line in f1.readlines():
            if line[0] == '#':
                continue
            s = line.strip().split('\t')
            v_pos = int(s[0])
            mutated_genome.append(ref_seq[pos: v_pos - 1])
            mutated_genome.append(s[2])
            pos = v_pos + len(s[1]) - 1
    mutated_genome.append(ref_seq[pos:])
    mut_str = ''.join(mutated_genome)
    with open(fasta_path + sample_id + '.fasta', 'w') as f:
        f.write('>' + sample_id + '\n')
        f.write(mut_str + '\n')
    cmd = path_to_minimap + ' -ax asm5 ' + path_to_ref + ' ' + fasta_path + sample_id + '.fasta | samtools sort - | samtools view -bS - > ' + bam_path + \
          sample_id + '.bam'
    check_call(cmd, shell=True)
    check_call('samtools index ' + bam_path + sample_id + '.bam', cwd=bam_path, shell=True)
    samfile = pysam.AlignmentFile(bam_path + sample_id + '.bam', "rb")
    raw_variants = []
    for read in samfile.fetch():
        ref_pos = read.reference_start
        read_pos = read.query_alignment_start
        # print('ref_pos = %d, contig_pos =% d' % (ref_pos, read_pos))
        # print(read.cigarstring)
        for op, l in read.cigartuples:
            if op == 0:
                #match
                # print(ref_seq[ref_pos: ref_pos + l])
                # print(mut_str[read_pos: read_pos + l])
                for i in range(l):
                    if ref_seq[ref_pos] != mut_str[read_pos]:
                        raw_variants.append(Variant(ref_pos + 1, ref_seq[ref_pos], mut_str[read_pos]))
                    ref_pos += 1
                    read_pos += 1
            elif op == 1:
                #insert
                raw_variants.append(Variant(ref_pos + 1, '', mut_str[read_pos: read_pos + l]))
                read_pos += l
            elif op == 2:
                #delete
                # raw_variants.append(Variant(ref_pos, ref_seq[ref_pos: ref_pos + l], ''))
                for i in range(l):
                    raw_variants.append(Variant(ref_pos + i + 1, ref_seq[ref_pos + i], ''))
                ref_pos += l
            elif op == 5:
                #hardclip
                read_pos += l
    samfile.close()
    # check_call('rm ' + temp_path + sample_id + '.fasta', shell=True)
    # check_call('rm ' + temp_path + sample_id + '.sam', shell=True)
    return sample_id, raw_variants


def gen_raw_variants_old(sample_id, ref_seq):
    if not exists(path_to_variants + sample_id + '.variants'):
        return sample_id, []
    if not overwrite and exists(out_path + sample_id + '.variants'):
        return sample_id, []
    mutated_genome = []
    pos = 0
    with open(path_to_variants + sample_id + '.variants', 'r') as f1:
        for line in f1.readlines():
            if line[0] == '#':
                continue
            s = line.strip().split('\t')
            v_pos = int(s[0])
            mutated_genome.append(ref_seq[pos: v_pos - 1])
            if s[-1] == 'snp':
                mutated_genome.append(s[1])
                pos = v_pos
            elif s[-1] == 'ins':
                mutated_genome.append(s[1])
                pos = v_pos - 1
            else:
                pos = v_pos
    mutated_genome.append(ref_seq[pos:])
    mut_str = ''.join(mutated_genome)
    if not exists(fasta_path + sample_id + '.fasta'):
        with open(fasta_path + sample_id + '.fasta', 'w') as f:
            f.write('>' + sample_id + '\n')
            f.write(mut_str + '\n')
    if not exists(bam_path + sample_id + '.bam'):
        cmd = path_to_minimap + ' -ax asm5 ' + path_to_ref + ' ' + fasta_path + sample_id + '.fasta | samtools sort - | samtools view -bS - > ' + bam_path + \
              sample_id + '.bam'
        check_call(cmd, shell=True)
        check_call('samtools index ' + bam_path + sample_id + '.bam', cwd=bam_path, shell=True)
    samfile = pysam.AlignmentFile(bam_path + sample_id + '.bam', "rb")
    raw_variants = []
    for read in samfile.fetch():
        ref_pos = read.reference_start
        read_pos = read.query_alignment_start
        # print('ref_pos = %d, contig_pos =% d' % (ref_pos, read_pos))
        # print(read.cigarstring)
        for op, l in read.cigartuples:
            if op == 0:
                #match
                # print(ref_seq[ref_pos: ref_pos + l])
                # print(mut_str[read_pos: read_pos + l])
                for i in range(l):
                    if ref_seq[ref_pos] != mut_str[read_pos]:
                        raw_variants.append(Variant(ref_pos + 1, ref_seq[ref_pos], mut_str[read_pos]))
                    ref_pos += 1
                    read_pos += 1
            elif op == 1:
                #insert
                raw_variants.append(Variant(ref_pos + 1, '', mut_str[read_pos: read_pos + l]))
                read_pos += l
            elif op == 2:
                #delete
                # raw_variants.append(Variant(ref_pos, ref_seq[ref_pos: ref_pos + l], ''))
                for i in range(l):
                    raw_variants.append(Variant(ref_pos + i + 1, ref_seq[ref_pos + i], ''))
                ref_pos += l
            elif op == 5:
                #hardclip
                read_pos += l
    samfile.close()
    # check_call('rm ' + temp_path + sample_id + '.fasta', shell=True)
    # check_call('rm ' + temp_path + sample_id + '.sam', shell=True)
    return sample_id, raw_variants


def format_variants(sample_id, variants, ref_seq):
    if not overwrite and exists(out_path + sample_id + '.variants'):
        return sample_id, []
    res = []
    ref_seq_len = len(ref_seq)
    gene_name_to_variants = {}
    gene_name_to_cds = {}
    deletions = {}
    ref_seq_compl = str(Seq(ref_seq).reverse_complement())
    with open(out_path + sample_id + '.variants', 'w') as f:
        for var in variants:
            if var.ref == '':
                v_type = 'ins'
            elif var.alt == '':
                v_type = 'del'
            else:
                v_type = 'snp'
            if len(var.cds_list) == 0:
                process_non_cds(var.pos, var.alt, v_type, ref_seq_len, ref_seq_compl, deletions, res)
            else:
                for cds in var.cds_list:
                    if cds.type == CDSType.Gene and (not cds.is_pseudogene or translate_pseudogene):
                        var_list = gene_name_to_variants.get(cds.name)
                        if var_list is None:
                            gene_name_to_variants[cds.name] = [var]
                            gene_name_to_cds[cds.name] = cds
                        else:
                            var_list.append(var)
                    else:
                        process_non_gene(cds, var.pos, var.alt, v_type, ref_seq, ref_seq_compl, ref_seq_len, deletions, res)

        if not split_dels:
            append_dels(res, deletions, ref_seq, ref_seq_compl)

        # working with protein coding genes
        with open(broken_genes_path + sample_id + '.genes', 'w') as fout:
            for gene_name, var_list in gene_name_to_variants.items():
                snp_only = True
                for variant in var_list:
                    if variant.ref == '' or variant.alt == '':
                        snp_only = False
                        break
                gene = gene_name_to_cds[gene_name]
                if snp_only:
                    process_gene_snp_only(gene, var_list, ref_seq, res, fout)
                else:
                    process_gene_with_indels(gene, var_list, ref_seq, res, fout)
        if len(res) > 0:
            f.write('\n'.join(res))
            f.write('\n')
    return sample_id, res


def localize_all_variants(var_list: 'list[Variant]', cds_list: 'list[CDS]', keep_genes_set: 'set[str]'=None):
    """
    Finds CDS for each gene coord from all_snps list.

    :param pos_list: list of genome coords
    :param cds_list: list of CDS
    :param keep_genes_set: a set of CDS names, if non None - the others will be filtered out
    :return: coord to CDS dictionary
    """
    var_list.sort()
    if keep_genes_set is not None:
        filtered_cds_list = [cds for cds in cds_list if cds.name in keep_genes_set]
    else:
        filtered_cds_list = cds_list
    cds_num = len(filtered_cds_list)
    var_num = len(var_list)
    if var_num == 0:
        return
    cds_i = 0
    var_i = 0
    while True:
        var = var_list[var_i]
        curr_cds = filtered_cds_list[cds_i]
        if var.pos >= curr_cds.start:
            if var.pos >= curr_cds.end:
                cds_i += 1
                if cds_i == cds_num:
                    break
                continue
            var.cds_list.append(curr_cds)
            if len(var.ref) > 1:
                # deletion or complex event
                for j in range(cds_i + 1, cds_num):
                    cds_j = filtered_cds_list[j]
                    if var.pos + len(var.ref) > cds_j.start:
                        var.cds_list.append(cds_j)
                    else:
                        break
        var_i += 1
        if var_i == var_num:
            break


def main():
    # mp.set_start_method('forkserver')
    if not exists(out_path):
        makedirs(out_path)
    if not exists(broken_genes_path):
        makedirs(broken_genes_path)
    if not exists(fasta_path):
        makedirs(fasta_path)
    if not exists(bam_path):
        makedirs(bam_path)
    sample_to_variants = {}
    h37rv = read_h37rv()

    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]

    cds_list = read_annotations(upstream_length, filter_by_gene_len=False)

    tasks = Parallel(n_jobs=thread_num)(
        delayed(gen_raw_variants)(sample_id, h37rv) for sample_id in sample_ids)
    pos_to_variants = {}
    for sample_id, variants in tasks:
        if len(variants) == 0:
            continue
        vars = []
        sample_to_variants[sample_id] = vars
        for var in variants:
            var_l = pos_to_variants.get(var.pos)
            if var_l is None:
                pos_to_variants[var.pos] = [var]
                vars.append(var)
            else:
                found = False
                for v in var_l:
                    if v.ref == var.ref and v.alt == var.alt:
                        found = True
                        vars.append(v)
                        break
                if not found:
                    var_l.append(var)
                    vars.append(var)

    var_list = []
    for pos, var_l in pos_to_variants.items():
        var_list.extend(var_l)
    var_list.sort()
    print('done with variant reading')

    localize_all_variants(var_list, cds_list)
    print('done with variant localization')

    tasks = Parallel(n_jobs=thread_num)(
        delayed(format_variants)(sample_id, variants, h37rv)
        for sample_id, variants in sample_to_variants.items())
    # formatted_snps = {sample_id: variants for sample_id, variants in tasks}
    print('done with variant format')

    # all_formatted_snps = set()
    # for f_snps in formatted_snps.values():
    #     all_formatted_snps.update(f_snps)
    # all_formatted_snps = list(all_formatted_snps)
    # all_formatted_snps.sort()
    #
    # with open(out_path_snp, 'w') as f:
    #     f.write('\n'.join(all_formatted_snps))
    #     f.write('\n')
    # print('printed all snps')


if __name__ == '__main__':
    main()
