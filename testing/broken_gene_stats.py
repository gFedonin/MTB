from os import makedirs, listdir
from os.path import exists

from Bio.Seq import Seq
from sklearn.externals.joblib import Parallel, delayed

from core.annotations import read_annotations, localize_all_variants, CDSType
from core.constants import data_path, upstream_length
from core.data_reading import read_h37rv
from data_processing.print_mutations_NW import get_aminoacids_sense, get_aminoacids_antisense

alignment_path = data_path + 'snps/gatk_bwa_mem_rev/'
path_to_variants = alignment_path + 'raw_variants_no_mq/'
out_path = alignment_path + 'broken_genes/'
out_path_stat = alignment_path + 'broken_genes.stat'

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


def process_gene_snp_only(gene, var_list, ref_seq, fout):
    j = 0
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
                        fout.write(gene.name + '\n')
                        fout.write('reference in pos %d:\n' % protein_pos)
                        fout.write(aa_ref)
                        fout.write('\nmutated:\n')
                        fout.write(aa_alt)
                        fout.write('\n')
                        break
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
                        fout.write(gene.name + '\n')
                        fout.write('reference in pos %d:\n' % protein_pos)
                        fout.write(aa_ref)
                        fout.write('\nmutated:\n')
                        fout.write(aa_alt)
                        fout.write('\n')
                        break
            j += 1


class Variant:

    def __init__(self, pos, alt, v_type):
        self.pos = pos
        self.alt = alt
        self.v_type = v_type

    def __eq__(self, other):
        return self.pos == other.pos and self.v_type == other.v_type

    def __lt__(self, other):
        if self.pos == other.pos:
            if self.v_type != 'ins':
                return True
            else:
                return False
        return self.pos < other.pos


def process_gene_with_indels(gene, var_list, ref_seq, fout):

    seq_new = []
    new_len = 0
    s = gene.start - 1
    # lens = []
    var_list.sort()
    if gene.name == 'Rv3785' and '921' in fout.name:
        a = 0
    vars = []
    for pos, alt, v_type in var_list:
        vars.append(Variant(pos, alt, v_type))
    vars = sorted(vars)
    for var in vars:
        if var.v_type == 'snp':
            seq_new.append(ref_seq[s: var.pos - 1])
            seq_new.append(var.alt)
            new_len += var.pos - s
            # lens.append(var.pos - s)
            s = var.pos
        elif var.v_type == 'ins':
            seq_new.append(ref_seq[s: var.pos])
            seq_new.append(var.alt)
            new_len += var.pos - s + len(var.alt)
            # lens.append(var.pos - s + len(var.alt))
            s = var.pos
        else:
            # del
            if s != var.pos - 1:
                seq_new.append(ref_seq[s: var.pos - 1])
                new_len += var.pos - 1 - s
                # lens.append(var.pos - 1 - s)
            s = var.pos
    seq_new.append(ref_seq[s:gene.end])
    new_len += gene.end - s
    # lens.append(gene.end - s)
    seq_len = len(Seq(''.join(seq_new)))
    if new_len != seq_len:
        print('AAAAAA!!!!')
    # additional_codon_num = 3 * (int(add_codons_to_frameshift/2 * (gene.end - gene.start + 1))//3)
    # if additional_codon_num % 3 != 0:
    #     print('ups!!!')
    if new_len % 3 != 0:
        # res.append('\t'.join(
        #     (gene.type.name, gene.name, '-', '-', '-', 'FS')))
        # continue
        # add some more letters to the end, might be new stop codon
        # seq_new.insert(0, ref_seq[gene.start - additional_codon_num: gene.start])
        # seq_new.append(ref_seq[gene.end: gene.end + 3 - new_len % 3 + additional_codon_num])
        if gene.strand == 1:
            seq_new.append(ref_seq[gene.end: gene.end + 3 - new_len % 3])
        else:
            seq_new.insert(0, ref_seq[gene.start - (3 - new_len % 3): gene.start])
    # else:
    #     seq_new.insert(0, ref_seq[gene.start - additional_codon_num: gene.start])
    #     seq_new.append(ref_seq[gene.end: gene.end + additional_codon_num])
    seq_old = Seq(ref_seq[gene.start - 1: gene.end])
    if len(seq_old) % 3 != 0:
        print(gene.name + ' wrong annotation')
        return
    seq_new = Seq(''.join(seq_new))
    if len(seq_new) % 3 != 0:
        print(gene.name + ' strange length')
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
        if start == -1 and translated_new[i] == translated_old[0]:
            start = i
        if translated_new[i] == '*':
            stop = i
            break
    if start == -1 or stop == -1 or stop - start < shortening_threshold * len(translated_old):
        fout.write(gene.name + '\n')
        fout.write('reference:\n')
        fout.write(str(translated_old))
        fout.write('\nmutated_extended:\n')
        fout.write(str(translated_new))
        fout.write('\nmutated_start-end:\n')
        if start == -1:
            fout.write('no start')
        elif stop == -1:
            fout.write('no stop')
        else:
            fout.write(str(translated_new[start:stop]))
        fout.write('\n')
        return


def format_variants(sample_id, variants, ref_seq, snp_to_cds):

    genes_to_variants = {}
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

    # working with protein coding genes
    with open(out_path + sample_id + '.genes', 'w') as fout:
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
                process_gene_snp_only(gene, var_list, ref_seq, fout)
            else:
                process_gene_with_indels(gene, var_list, ref_seq, fout)
    return sample_id


def parse_variants():
    if not exists(out_path):
        makedirs(out_path)
    sample_to_snps = {}
    all_snp_pos = set()
    h37rv = read_h37rv()

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
        delayed(format_variants)(sample_id, variants, h37rv, snp_to_cds)
        for sample_id, variants in sample_to_snps.items())
    print('done with variant parsing')


def compute_stat():
    sample_ids = [fname[0:-len('.variants')] for fname in listdir(path_to_variants) if fname.endswith('.variants')]
    cds_list = read_annotations(upstream_length, filter_by_gene_len=False)
    gene_names = {}
    for cds in cds_list:
        if cds.type == CDSType.Gene:
            gene_names[cds.name] = 0
    for sample_id in sample_ids:
        lines = open(out_path + sample_id + '.genes').readlines()
        for i in range(0, len(lines), 7):
            gene_names[lines[i].strip()] += 1
    with open(out_path_stat, 'w') as f:
        for name, c in gene_names.items():
            f.write(name + '\t' + str(c) + '\n')


if __name__ == '__main__':
    parse_variants()
    compute_stat()