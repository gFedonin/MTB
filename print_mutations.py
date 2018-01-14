from bisect import bisect_left, bisect_right
from enum import Enum

from Bio import SeqIO
from sklearn.externals.joblib import Parallel, delayed

import numpy as np

from constants import complement, codon_table

path_to_ids = './data/Full_subset_filtered_snp_pheno2.txt'
path_to_snps = '/export/data/kchukreev/data/mutations_files/unformated_vcfs/realigned_vcfs/'
out_path = './data/snps/'
out_path_snp = './data/all_snp_list_filtered2.csv'
path_to_depths = './data/coverages_with_percentiles/5/'
path_to_annotations = './data/AL123456_rev.gff'
path_to_ref = './data/h37rv.fasta'


qual_thresold = 40
snp_cov_threshold = 15
sample_cov_threshold = 0.9

filter_by_upper_percentile = True
upstream_length = 100


def read_h37rv():
    fasta_sequences = SeqIO.parse(open(path_to_ref), 'fasta')
    sequence = ''
    for fasta in fasta_sequences:
        sequence = str(fasta.seq)
    return sequence.upper()


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


def read_annotations():
    genes = {}
    nc_rna = {}
    r_rna = {}
    t_rna = {}
    upstream = {}
    repeats_and_mobile_elements = []

    with open(path_to_annotations, 'r') as f:
        for line in f.readlines():
            s = line.split('\t')
            type = s[2]
            start = int(s[3])
            end = int(s[4])
            if type in ('Repeat_region', 'mobile_element'):
                repeats_and_mobile_elements.append((start, end))
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
    repeats_and_mobile_elements.sort(key=lambda tup: tup[0])
    return cds, repeats_and_mobile_elements


def read_snps(sample_id, repeats, repeats_starts, percentile):
    snps = []
    with open(path_to_snps + sample_id + '_h37rv.vcf', 'r') as f1:
        lines = f1.readlines()
        for line in lines:
            line = line[:-1]
            if line[0] != '#':
                tokens = line.split('\t')
                qual = float(tokens[5])
                if qual >= qual_thresold:
                    info = tokens[-3]
                    i = info.index('DP=')
                    j = info.index(';', i + 3)
                    cov = info[i + 3: j]
                    if int(cov) > snp_cov_threshold:
                        if filter_by_upper_percentile and cov > percentile:
                            continue
                        snp_pos = int(tokens[1])
                        in_repeat = False
                        i = bisect_left(repeats_starts, snp_pos)
                        if i == 0:
                            if snp_pos == repeats[0][0]:
                                in_repeat = True
                        else:
                            rep = repeats[i - 1]
                            if rep[1] >= snp_pos:
                                in_repeat = True
                        if not in_repeat:
                            snps.append((snp_pos, tokens[4]))
    return sample_id, snps


def read_coverage(sample_id):
    coverage = []
    with open(path_to_depths + sample_id + '.coverage') as f:
        f.readline()
        percentile = f.readline().split('\t')[-1]
        f.readline()
        for line in f.readlines():
            s = line.split()
            coverage.append((int(s[0]), int(s[1])))
    coverage.sort(key=lambda tup: tup[0])
    return sample_id, coverage, percentile


def sample_cov(sample_id, all_snps, coverage):
    # coverage = read_coverage(sample_id)
    cov = 0
    lo = 0
    # print(sample_id + ' ' + str(len(coverage)))
    for start, end in coverage:
        s = bisect_left(all_snps, start, lo=lo)
        e = bisect_right(all_snps, end, lo=s + 1)
        # print(str(start) + ' ' + str(end) + ' ' + str(s) + ' ' + str(e))
        cov += e - s
        lo = e
    # print(sample_id + ' cov ' + str(cov))
    return sample_id, cov / len(all_snps)


def translate(ref_seq, strand, start_pos, alt1, alt2, alt3):
    new_triplet = ''
    if strand == 1:


        old_triplet = ref_seq[start_pos:start_pos + 3]

        if alt1 != '-':
            new_triplet += alt1
        else:
            new_triplet += ref_seq[start_pos]

        if alt2 != '-':
            new_triplet += alt2
        else:
            new_triplet += ref_seq[start_pos + 1]

        if alt3 != '-':
            new_triplet += alt3
        else:
            new_triplet += ref_seq[start_pos + 2]

    else:

        old_triplet = ''
        for i in range(start_pos - 2, start_pos + 1):
            old_triplet += complement[ref_seq[i]]

        if alt3 != '-':
            new_triplet += complement[alt3]
        else:
            new_triplet += complement[ref_seq[start_pos]]

        if alt2 != '-':
            new_triplet += complement[alt2]
        else:
            new_triplet += complement[ref_seq[start_pos - 1]]

        if alt1 != '-':
            new_triplet += complement[alt1]
        else:
            new_triplet += complement[ref_seq[start_pos - 2]]

    return codon_table[old_triplet], codon_table[new_triplet]


def compute_snp_cov(all_snps, coverage):
    snp_cov = np.zeros(len(all_snps), dtype=int)
    lo = 0
    for start, end in coverage:
        s = bisect_left(all_snps, start, lo=lo)
        if s < len(all_snps):
            e = bisect_right(all_snps, end, lo=s + 1)
            # print(str(start) + ' ' + str(end) + ' ' + str(s) + ' ' + str(e))
            snp_cov[s:e] = 1
            lo = e
        else:
            break
    return snp_cov


def get_aminoacids_sence(ref_seq, nucleotide_pos, snps, i):
    snp_num = len(snps)
    pos, alt = snps[i]
    if nucleotide_pos == 0:
        if i != snp_num - 1:
            pos1, alt1 = snps[i + 1]
            if pos1 == pos + 1:
                pos2, alt2 = snps[i + 2]
                if i != snp_num - 2 and pos2 == pos + 2:
                    return translate(ref_seq, "+", pos - 1, alt, alt1, alt2), i + 2
                else:
                    return translate(ref_seq, "+", pos - 1, alt, alt1, '-'), i + 1
            elif pos1 == pos + 2:
                return translate(ref_seq, "+", pos - 1, alt, '-', snps[i + 2][1]), i + 1
        return translate(ref_seq, "+", pos - 1, alt, '-', '-'), i
    elif nucleotide_pos == 1:
        # if i != 0:
        #     pos_1, alt_1 = snps[i - 1]
        #     if pos_1 == pos - 1:
        #         if i != snp_num - 1 and snps[i + 1][0] == pos + 1:
        #             return translate(ref_seq, "+", pos - 2, alt_1, alt, snps[i + 1][1])
        #         else:
        #             return translate(ref_seq, "+", pos - 2, alt_1, alt, '-')
        if i != snp_num - 1 and snps[i + 1] == pos + 1:
            return translate(ref_seq, "+", pos - 2, '-', alt, snps[i + 1][1]), i + 1
        else:
            return translate(ref_seq, "+", pos - 2, '-', alt, '-'), i
    else:
        # if i != 0:
        #     pos_1, alt_1 = snps[i - 1]
        #     if pos_1 == pos - 1:
        #         if i != 1 and snps[i - 2][0] == pos - 2:
        #             return translate(ref_seq, "+", pos - 3, snps[i - 2][1], alt_1, alt)
        #         else:
        #             return translate(ref_seq, "+", pos - 3, '-', alt_1, alt)
        #     elif pos_1 == pos - 2:
        #         return translate(ref_seq, "+", pos - 3, snps[i - 2][1], '-', alt)
        return translate(ref_seq, "+", pos - 3, '-', '-', alt), i


def get_aminoacids_antisence(ref_seq, nucleotide_pos, snps, i):
    snp_num = len(snps)
    pos, alt = snps[i]
    if nucleotide_pos == 2:
        if i != snp_num - 1:
            pos1, alt1 = snps[i + 1]
            if pos1 == pos + 1:
                if i != snp_num - 2 and snps[i + 2][0] == pos + 2:
                    return translate(ref_seq, "-", pos + 1, alt, alt1, snps[i + 2][1]), i + 2
                else:
                    return translate(ref_seq, "-", pos + 1, alt, alt1, '-'), i + 1
            elif pos1 == pos + 2:
                return translate(ref_seq, "-", pos + 1, alt, '-', alt1), i + 1
        return translate(ref_seq, "-", pos + 1, alt, '-', '-'), i
    elif nucleotide_pos == 1:
        # if i != 0:
        #     pos_1, alt_1 = snps[i - 1]
        #     if pos_1 == pos - 1:
        #         if i != snp_num - 1 and snps[i + 1][0] == pos + 1:
        #             return translate(ref_seq, "-", pos, alt_1, alt, snps[i + 1][1])
        #         else:
        #             return translate(ref_seq, "-", pos, alt_1, alt, '-')
        if i != snp_num - 1 and snps[i + 1][0] == pos + 1:
            return translate(ref_seq, "-", pos, '-', alt, snps[i + 1][1]), i + 1
        return translate(ref_seq, "-", pos, '-', alt, '-'), i
    else:
        # if i != 0:
        #     pos_1, alt_1 = snps[i - 1]
        #     if pos_1 == pos - 1:
        #         if i != 1 and snps[i - 2][0] == pos - 2:
        #             return translate(ref_seq, "-", pos - 1, snps[i - 2][1], alt_1, alt)
        #         else:
        #             return translate(ref_seq, "-", pos - 1, '-', alt_1, alt)
        #     elif pos_1 == pos - 2:
        #         return translate(ref_seq, "-", pos - 1, snps[i - 2][1], '-', alt)
        return translate(ref_seq, "-", pos - 1, '-', '-', alt), i


def format_snp(sample_id, snps, ref_seq, cds_list, cds_starts, uncovered_snp_pos):

    res = []
    j = 0
    while j < len(snps):
        pos, alt = snps[j]
        if pos in uncovered_snp_pos:
            continue
        i = bisect_left(cds_starts, pos)
        cds = None
        if i == 0:
            if cds_list[0].start == pos:
                cds = cds_list[0]
        else:
            if cds_list[i - 1].end >= pos:
                cds = cds_list[i - 1]

        if cds is not None:
            if cds.strand == 1:
                if cds.type == CDSType.Gene:
                    protein_pos = (pos - cds.start) // 3 + 1
                    nucleotide_pos = (pos - cds.start) % 3

                    # NEIGHBOURS ANALYSIS
                    (old_aminoacid, new_aminoacid), j = get_aminoacids_sence(ref_seq, nucleotide_pos, snps, j)

                    if old_aminoacid != new_aminoacid:
                        res.append((cds.type.name, cds.name, str(protein_pos), old_aminoacid, new_aminoacid))
                else:
                    coord = pos - cds.start
                    res.append((cds.type.name, cds.name, coord, ref_seq[pos - 1], alt))
            else:  # if strand is '-'
                if cds.type == CDSType.Gene:
                    protein_pos = (cds.end - pos) // 3 + 1
                    nucleotide_pos = (cds.end - pos) % 3

                    # NEIGHBOURS ANALYSIS
                    (old_aminoacid, new_aminoacid), j = get_aminoacids_antisence(ref_seq, nucleotide_pos, snps, j)

                    if old_aminoacid != new_aminoacid:
                        res.append((cds.type.name, cds.name, str(protein_pos), old_aminoacid, new_aminoacid))
                else:
                    coord = cds.end - (pos - 1)
                    res.append((cds.type.name, cds.name, coord, complement[ref_seq[pos - 1]], complement[alt]))
        else:
            res.append(('non_cds', ' ', str(len(ref_seq) - pos + 1), complement[ref_seq[pos - 1]], complement[alt]))
        j += 1
    return sample_id, res


def main():
    sample_to_snps = {}
    all_snp_pos = set()
    sample_to_cov = {}
    sample_to_cov_list = {}
    percentiles = {}

    h37rv = read_h37rv()

    sample_ids = [sample_id[:-1] for sample_id in open(path_to_ids, 'r').readlines()]
    sample_num = len(sample_ids)

    cds, repeats_and_mobile_elements = read_annotations()
    repeats_starts = [x[0] for x in repeats_and_mobile_elements]
    cds_starts = [x.start for x in cds]

    coverages = Parallel(n_jobs=-1)(delayed(read_coverage)(sample_id) for sample_id in sample_ids)
    for sample_id, cov, percentile in coverages:
        sample_to_cov_list[sample_id] = cov
        percentiles[sample_id] = percentile

    tasks = Parallel(n_jobs=-1)(
        delayed(read_snps)(sample_id, repeats_and_mobile_elements, repeats_starts, percentiles[sample_id]) for sample_id in sample_ids)
    for sample_id, f_snps in tasks:
        sample_to_snps[sample_id] = f_snps
        for snp_pos, alt in f_snps:
            all_snp_pos.add(snp_pos)
    all_snps = list(all_snp_pos)
    all_snps.sort()
    snp_num = len(all_snps)
    print('done with snp')

    tasks = Parallel(n_jobs=-1)(delayed(sample_cov)(sample_id, all_snps, sample_to_cov_list[sample_id]) for sample_id in sample_ids)
    for sample_id, cov in tasks:
        sample_to_cov[sample_id] = cov
    print('done with sample coverages')

    snp_cov_list = Parallel(n_jobs=-1)(
        delayed(compute_snp_cov)(all_snps, sample_to_cov_list[sample_id]) for sample_id in sample_ids)
    snp_cov = np.zeros(snp_num, dtype=int)
    for snp_cov_arr in snp_cov_list:
        snp_cov += snp_cov_arr
    print('done with snp coverages by samples')

    # with open('./snp_cov.txt', 'w')as f:
    #     for i in range(snp_num):
    #         f.write(str(all_snps[i]) + ' ' + str(snp_cov[i]) + '\n')

    uncovered_snp_pos = set()
    for i in range(snp_num):
        if snp_cov[i]/sample_num < sample_cov_threshold:
            uncovered_snp_pos.add(all_snps[i])

    formatted_snps = Parallel(n_jobs=-1)(
        delayed(format_snp)(sample_id, snps, h37rv, cds, cds_starts, uncovered_snp_pos)
        for sample_id, snps in sample_to_snps.items() if sample_to_cov[sample_id] >= sample_cov_threshold
    )
    print('done with snp format')

    all_formatted_snps = set()
    for sample_id, f_snps in formatted_snps:
        with open(out_path + sample_id + '.snp', 'w') as f:
            for formatted_snp in f_snps:
                all_formatted_snps.add(formatted_snp)
                f.write('\t'.join(formatted_snp) + '\n')

    with open(out_path_snp, 'w') as f:
        for formatted_snp in all_formatted_snps:
            f.write(str(formatted_snp) + '\n')
    print('printed all snps')


if __name__ == '__main__':
    main()
