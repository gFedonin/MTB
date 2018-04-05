from os import mkdir

from Bio import SeqIO
from Bio.Seq import Seq
from os.path import exists
from sklearn.externals.joblib import Parallel, delayed

from annotations import CDSType, read_annotations, localize_all_snps
from constants import codon_table, codon_table_compl, complement

path_to_ids = './data/all_with_pheno_and_snp.txt'
path_to_snps = './data/snps/raw_with_DR_with_pheno_and_snp_mc5/'
out_path = './data/snps/annotated_with_pheno_and_snp_mc5/'
out_path_snp = out_path + 'all_snp_list.csv'
path_to_annotations = './data/AL123456_rev.gff'
path_to_ref = './data/h37rv.fasta'

upstream_length = 100


def read_h37rv():
    fasta_sequences = SeqIO.parse(open(path_to_ref), 'fasta')
    sequence = ''
    for fasta in fasta_sequences:
        sequence = str(fasta.seq)
    return sequence.upper()


def read_snps(sample_id):
    snps = []
    with open(path_to_snps + sample_id + '.snp', 'r') as f1:
        lines = f1.readlines()
        for line in lines:
            s = line.strip().split('\t')
            snps.append((int(s[0]), s[1]))
    return sample_id, snps


def get_aminoacids_sense(ref_seq, nucleotide_pos, variants, i):
    var_num = len(variants)
    pos, alt, v_type = variants[i]
    if nucleotide_pos == 0:
        aa0 = codon_table[ref_seq[pos - 1:pos + 2]]
        if i != var_num - 1:
            pos1, alt1, v_type1 = variants[i + 1]
            if v_type1 == 'snp':
                if pos1 == pos + 1:
                    if i != var_num - 2 and variants[i + 2][0] == pos + 2:
                        aa1 = codon_table[alt + alt1 + variants[i + 2][1]]
                        return aa0, aa1, i + 2
                    else:
                        aa1 = codon_table[alt + alt1 + ref_seq[pos + 1]]
                        return aa0, aa1, i + 1
                elif pos1 == pos + 2:
                    aa1 = codon_table[alt + ref_seq[pos] + alt1]
                    return aa0, aa1, i + 1
            # elif v_type1 == 'ins':
            #
        aa1 = codon_table[alt + ref_seq[pos] + ref_seq[pos + 1]]
        return aa0, aa1, i
    elif nucleotide_pos == 1:
        aa0 = codon_table[ref_seq[pos - 2:pos + 1]]
        if i != var_num - 1 and variants[i + 1][0] == pos + 1:
            aa1 = codon_table[ref_seq[pos - 2] + alt + variants[i + 1][1]]
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
    pos, alt = snps[i]
    if nucleotide_pos == 2:
        aa0 = codon_table_compl[ref_seq[pos - 1:pos + 2]]
        if i != snp_num - 1:
            pos1, alt1 = snps[i + 1]
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


def format_insert(ref_seq, nucleotide_pos, pos, alt):
    ins_len = len(alt)
    if ins_len % 3 != 0:
        return 'FS'
    if nucleotide_pos == 2:
        j = 0
        ins = []
        for n in range(ins_len // 3):
            ins.append(codon_table[alt[3*j: 3 + 3*j]])
            j += 1
        return [('', ''), ''.join(ins)]
    elif nucleotide_pos == 1:
        aa0 = codon_table[ref_seq[pos - 2: pos + 1]]
        aa1 = codon_table[ref_seq[pos - 2: pos] + alt[0]]
        j = 0
        ins = []
        for n in range((ins_len - 1) // 3):
            ins.append(codon_table[alt[1 + 3*j: 4 + 3*j]])
            j += 1
        ins.append(codon_table[alt[-2:] + ref_seq[pos]])
        return [(aa0, aa1), ''.join(ins)]
    elif nucleotide_pos == 0:
        aa0 = codon_table[ref_seq[pos - 1: pos + 2]]
        aa1 = codon_table[ref_seq[pos - 1] + alt[0:1]]
        j = 0
        ins = []
        for n in range((ins_len - 2) // 3):
            ins.append(codon_table[alt[2 + 3*j: 5 + 3*j]])
            j += 1
        ins.append(codon_table[alt[-1] + ref_seq[pos: pos + 2]])
        return [(aa0, aa1), ''.join(ins)]


def format_insert_antisense(ref_seq, nucleotide_pos, pos, alt):
    ins_len = len(alt)
    if ins_len % 3 != 0:
        return 'FS'
    if nucleotide_pos == 0:
        j = 0
        ins = []
        for n in range(ins_len // 3):
            ins.append(codon_table_compl[alt[3*j: 3 + 3*j]])
            j += 1
        ins.reverse()
        return (('', ''), ''.join(ins))
    elif nucleotide_pos == 1:
        aa0 = codon_table_compl[ref_seq[pos - 2: pos + 1]]
        aa1 = codon_table_compl[ref_seq[pos - 2: pos] + alt[0]]
        j = 0
        ins = []
        for n in range((ins_len - 1) // 3):
            ins.append(codon_table_compl[alt[1 + 3*j: 4 + 3*j]])
            j += 1
        ins.append(codon_table_compl[alt[-2:] + ref_seq[pos]])
        ins.reverse()
        return ((aa0, aa1), ''.join(ins))
    elif nucleotide_pos == 2:
        aa0 = codon_table_compl[ref_seq[pos - 1: pos + 2]]
        aa1 = codon_table_compl[ref_seq[pos - 1] + alt[0:1]]
        j = 0
        ins = []
        for n in range((ins_len - 2) // 3):
            ins.append(codon_table_compl[alt[2 + 3*j: 5 + 3*j]])
            j += 1
        ins.append(codon_table_compl[alt[-1] + ref_seq[pos: pos + 2]])
        ins.reverse()
        return ((aa0, aa1), ''.join(ins))



def format_delete(ref_seq, nucleotide_pos, variants, i, cds):
    var_num = len(variants)
    pos, alt, v_type = variants[i]
    start = pos
    end = pos
    j = i + 1
    prev_pos = pos
    while j < var_num and end <= cds.end:
        pos, alt, v_type = variants[j]
        if v_type != 'del' or pos != prev_pos + 1:
            break
        j += 1
        end += 1
        prev_pos += 1
    del_len = end - start
    if del_len % 3 != 0:
        return 'FS', j - 1
    res = []
    if nucleotide_pos == 0:
        j = 0
        for n in range(del_len // 3):
            res.append((codon_table[ref_seq[pos - 1 + 3*j: pos + 2 + 3*j]], '-'))
            j += 1
    elif nucleotide_pos == 1:
        aa0 = codon_table[ref_seq[end - 2: end + 1]]
        aa1 = codon_table[ref_seq[pos - 2] + ref_seq[end:end + 2]]
        j = 0
        for n in range(del_len // 3):
            res.append((codon_table[ref_seq[pos - 2 + 3*j: pos + 1 + 3*j]], '-'))
            j += 1
        res.append((aa0, aa1))
    elif nucleotide_pos == 2:
        aa0 = codon_table[ref_seq[end - 3: end]]
        aa1 = codon_table[ref_seq[pos - 3: pos - 1] + ref_seq[end]]
        j = 0
        for n in range(del_len // 3):
            res.append((codon_table[ref_seq[pos - 3 + 3*j: pos + 3*j]], '-'))
            j += 1
        res.append((aa0, aa1))
    return res, j - 1


def format_delete_antisense(ref_seq, nucleotide_pos, variants, i, cds):
    var_num = len(variants)
    pos, alt, v_type = variants[i]
    start = pos
    end = pos
    j = i + 1
    prev_pos = pos
    while j < var_num and end <= cds.end:
        pos, alt, v_type = variants[j]
        if v_type != 'del' or pos != prev_pos + 1:
            break
        j += 1
        end += 1
        prev_pos += 1
    del_len = end - start
    if del_len % 3 != 0:
        return 'FS', j - 1
    res = []
    if nucleotide_pos == 2:
        j = 0
        for n in range(del_len // 3):
            res.append((codon_table_compl[ref_seq[pos - 1 + 3*j: pos + 2 + 3*j]], '-'))
            j += 1
    elif nucleotide_pos == 1:
        aa0 = codon_table_compl[ref_seq[end - 2: end + 1]]
        aa1 = codon_table_compl[ref_seq[pos - 2] + ref_seq[end:end + 2]]
        j = 0
        for n in range(del_len // 3):
            res.append((codon_table_compl[ref_seq[pos - 2 + 3*j: pos + 1 + 3*j]], '-'))
            j += 1
        res.append((aa0, aa1))
    elif nucleotide_pos == 0:
        aa0 = codon_table_compl[ref_seq[end - 3: end]]
        aa1 = codon_table_compl[ref_seq[pos - 3: pos - 1] + ref_seq[end]]
        j = 0
        for n in range(del_len // 3):
            res.append((codon_table_compl[ref_seq[pos - 3 + 3*j: pos + 3*j]], '-'))
            j += 1
        res.append((aa0, aa1))
    return res, j - 1


def format_variant(sample_id, variants, ref_seq, ref_seq_compl, snp_to_cds):

    res = []
    j = 0
    ref_seq_len = len(ref_seq)
    with open(out_path + sample_id + '.variants', 'w') as f:
        while j < len(variants):
            pos, alt, v_type = variants[j]
            try:
                cds = snp_to_cds[pos]
            except:
                cds = None

            if cds is not None:
                if cds.strand == 1:
                    if cds.type == CDSType.Gene:
                        protein_pos = (pos - cds.start) // 3 + 1
                        nucleotide_pos = (pos - cds.start) % 3
                        if v_type == 'snp':
                            # NEIGHBOURS ANALYSIS
                            old_aminoacid, new_aminoacid, j = get_aminoacids_sense(ref_seq, nucleotide_pos, variants, j)

                            if old_aminoacid != new_aminoacid:
                                res.append('\t'.join(
                                    (cds.type.name, cds.name, str(protein_pos), old_aminoacid, new_aminoacid)))
                        elif v_type == 'ins':
                            insertion = format_insert(ref_seq, nucleotide_pos, pos, alt)
                            if type(insertion) is tuple:
                                (aa0, aa1), ins = insertion
                                if aa0 != aa1:
                                    res.append('\t'.join(
                                        (cds.type.name, cds.name, str(protein_pos), aa0, aa1, 'snp')))
                                res.append('\t'.join(
                                    (cds.type.name, cds.name, str(protein_pos), aa0, ins, 'ins')))
                            else:
                                res.append('\t'.join(
                                    (cds.type.name, cds.name, str(protein_pos), '-', '-', 'FS')))
                        else:
                            #del
                            deletion, j = format_delete(ref_seq, nucleotide_pos, variants, j, cds)
                            if type(deletion) is list:
                                for i in range(len(deletion)):
                                    aa0, aa1 = deletion[i]
                                    if aa1 == '-':
                                        res.append('\t'.join(
                                            (cds.type.name, cds.name, str(protein_pos + i), aa0, '-', 'del')))
                                    elif aa1 != aa0:
                                        res.append('\t'.join(
                                            (cds.type.name, cds.name, str(protein_pos + i), aa0, aa1, 'snp')))
                            else:
                                res.append('\t'.join(
                                    (cds.type.name, cds.name, str(protein_pos), '-', '-', 'FS')))
                    elif cds.type == CDSType.upstream:
                        res.append('\t'.join(
                            (cds.type.name, cds.name, str(pos - cds.end - 1), ref_seq[pos - 1], alt, v_type)))
                    else:
                        res.append('\t'.join(
                            (cds.type.name, cds.name, str(pos - cds.start + 1), ref_seq[pos - 1], alt, v_type)))
                else:  # if strand is '-'
                    if cds.type == CDSType.Gene:
                        protein_pos = (cds.end - pos) // 3 + 1
                        nucleotide_pos = (cds.end - pos) % 3
                        if v_type == 'snp':
                            # NEIGHBOURS ANALYSIS
                            old_aminoacid, new_aminoacid, j = get_aminoacids_antisense(ref_seq, nucleotide_pos, variants, j)

                            if old_aminoacid != new_aminoacid:
                                res.append('\t'.join(
                                    (cds.type.name, cds.name, str(protein_pos), old_aminoacid, new_aminoacid)))
                        elif v_type == 'ins':
                            insertion = format_insert_antisense(ref_seq, nucleotide_pos, pos, alt)
                            if type(insertion) is tuple:
                                (aa0, aa1), ins = insertion
                                if aa0 != aa1:
                                    res.append('\t'.join(
                                        (cds.type.name, cds.name, str(protein_pos), aa0, aa1, 'snp')))
                                res.append('\t'.join(
                                    (cds.type.name, cds.name, str(protein_pos), aa0, ins, 'ins')))
                            else:
                                res.append('\t'.join(
                                    (cds.type.name, cds.name, str(protein_pos), '-', '-', 'FS')))
                        else:
                            #del
                            deletion, j = format_delete_antisense(ref_seq, nucleotide_pos, variants, j, cds)
                            if type(deletion) is list:
                                for i in range(len(deletion)):
                                    aa0, aa1 = deletion[len(deletion) - i - 1]
                                    if aa1 == '-':
                                        res.append('\t'.join(
                                            (cds.type.name, cds.name, str(protein_pos + i), aa0, '-', 'del')))
                                    elif aa1 != aa0:
                                        res.append('\t'.join(
                                            (cds.type.name, cds.name, str(protein_pos + i), aa0, aa1, 'snp')))
                            else:
                                res.append('\t'.join(
                                    (cds.type.name, cds.name, str(protein_pos), '-', '-', 'FS')))
                    elif cds.type == CDSType.upstream:
                        res.append('\t'.join(
                            (cds.type.name, cds.name, str(cds.start - pos - 1), ref_seq_compl[pos - 1], complement[alt], v_type)))
                    else:
                        res.append('\t'.join(
                            (cds.type.name, cds.name, str(cds.end - (pos - 1)), ref_seq_compl[ref_seq_len - pos], complement[alt])))
            else:
                coord = ref_seq_len - pos
                res.append('\t'.join(('non_cds', ' ', str(coord + 1), ref_seq_compl[coord], complement[alt])))
            j += 1
        if len(res) > 0:
            f.write('\n'.join(res))
            f.write('\n')
    return res


def main():
    if not exists(out_path):
        mkdir(out_path)
    sample_to_snps = {}
    all_snp_pos = set()
    h37rv = read_h37rv()
    h37rv_compl = str(Seq(h37rv).reverse_complement())

    sample_ids = [sample_id[:-1] for sample_id in open(path_to_ids, 'r').readlines()]

    cds = read_annotations(path_to_annotations, upstream_length)

    tasks = Parallel(n_jobs=-1)(
        delayed(read_snps)(sample_id) for sample_id in sample_ids)
    for sample_id, snps in tasks:
        sample_to_snps[sample_id] = snps
        for snp_pos, alt in snps:
            all_snp_pos.add(snp_pos)
    all_snps = list(all_snp_pos)
    all_snps.sort()
    print('done with snp')

    snp_to_cds = localize_all_snps(all_snps, cds)

    formatted_snps = Parallel(n_jobs=-1)(
        delayed(format_variant)(sample_id, snps, h37rv, h37rv_compl, snp_to_cds)
        for sample_id, snps in sample_to_snps.items()
    )
    print('done with snp format')

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
