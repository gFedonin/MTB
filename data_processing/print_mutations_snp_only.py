from os import mkdir

from Bio import SeqIO
from Bio.Seq import Seq
from os.path import exists
from sklearn.externals.joblib import Parallel, delayed

from src.core.annotations import CDSType, localize_all_variants, read_annotations
from src.core.constants import codon_table_compl, complement, codon_table, upstream_length

path_to_ids = '../../data/all_with_pheno_and_snp.txt'
path_to_snps = '../../data/snps/raw_with_DR_with_pheno_and_snp_mc5/'
out_path = '../../data/snps/annotated_with_pheno_and_snp_mc5_snp_only/'
out_path_snp = out_path + 'all_var_list.csv'
path_to_annotations = '../../data/AL123456_rev.gff'
path_to_ref = '../../data/h37rv.fasta'


def read_h37rv():
    fasta_sequences = SeqIO.parse(open(path_to_ref), 'fasta')
    sequence = ''
    for fasta in fasta_sequences:
        sequence = str(fasta.seq)
    return sequence.upper()


def read_snps(sample_id):
    snps = []
    with open(path_to_snps + sample_id + '.variants', 'r') as f1:
        lines = f1.readlines()
        for line in lines:
            s = line.strip().split('\t')
            snps.append((int(s[0]), s[1]))
    return sample_id, snps


def get_aminoacids_sense(ref_seq, nucleotide_pos, variants, i):
    var_num = len(variants)
    pos, alt = variants[i]
    if nucleotide_pos == 0:
        aa0 = codon_table[ref_seq[pos - 1:pos + 2]]
        if i != var_num - 1:
            pos1, alt1 = variants[i + 1]
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


def format_variant(sample_id, variants, ref_seq, ref_seq_compl, snp_to_cds):

    res = []
    j = 0
    ref_seq_len = len(ref_seq)
    with open(out_path + sample_id + '.variants', 'w') as f:
        while j < len(variants):
            pos, alt = variants[j]
            try:
                cds = snp_to_cds[pos]
            except:
                cds = None

            if cds is not None:
                if cds.strand == 1:
                    if cds.type == CDSType.Gene:
                        protein_pos = (pos - cds.start) // 3 + 1
                        nucleotide_pos = (pos - cds.start) % 3
                        # NEIGHBOURS ANALYSIS
                        old_aminoacid, new_aminoacid, j = get_aminoacids_sense(ref_seq, nucleotide_pos, variants, j)

                        if old_aminoacid != new_aminoacid:
                            res.append('\t'.join(
                                (cds.type.name, cds.name, str(protein_pos), old_aminoacid, new_aminoacid, 'snp')))
                    elif cds.type == CDSType.upstream:
                        res.append('\t'.join(
                            (cds.type.name, cds.name, str(pos - cds.end - 1), ref_seq[pos - 1], alt, 'snp')))
                    else:
                        res.append('\t'.join(
                            (cds.type.name, cds.name, str(pos - cds.start + 1), ref_seq[pos - 1], alt, 'snp')))
                else:  # if strand is '-'
                    if cds.type == CDSType.Gene:
                        protein_pos = (cds.end - pos) // 3 + 1
                        nucleotide_pos = (cds.end - pos) % 3
                        # NEIGHBOURS ANALYSIS
                        old_aminoacid, new_aminoacid, j = get_aminoacids_antisense(ref_seq, nucleotide_pos, variants, j)

                        if old_aminoacid != new_aminoacid:
                            res.append('\t'.join(
                                (cds.type.name, cds.name, str(protein_pos), old_aminoacid, new_aminoacid, 'snp')))
                    elif cds.type == CDSType.upstream:
                        res.append('\t'.join(
                            (cds.type.name, cds.name, str(cds.start - pos - 1), ref_seq_compl[ref_seq_len - pos], complement[alt], 'snp')))
                    else:
                        res.append('\t'.join(
                            (cds.type.name, cds.name, str(cds.end - (pos - 1)), ref_seq_compl[ref_seq_len - pos], complement[alt], 'snp')))
            else:
                coord = ref_seq_len - pos
                res.append('\t'.join(('non_cds', ' ', str(coord + 1), ref_seq_compl[coord], complement[alt], 'snp')))
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

    sample_ids = [sample_id.strip() for sample_id in open(path_to_ids, 'r').readlines()]

    cds = read_annotations(path_to_annotations, upstream_length)

    tasks = Parallel(n_jobs=-1)(
        delayed(read_snps)(sample_id) for sample_id in sample_ids)
    for sample_id, snps in tasks:
        sample_to_snps[sample_id] = snps
        for snp_pos, alt in snps:
            all_snp_pos.add(snp_pos)
    all_snps = list(all_snp_pos)
    print('done with snp')

    snp_to_cds = localize_all_variants(all_snps, cds)

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
