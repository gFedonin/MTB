from Bio import SeqIO
import numpy as np
from bisect import bisect_left, bisect_right

from sklearn.externals.joblib import Parallel, delayed

path_to_ids = './data/Full_subset_filtered_snp_pheno2.txt'
path_to_snps = '/export/data/kchukreev/data/mutations_files/unformated_vcfs/realigned_vcfs/'
out_path_aln = './data/snp_aln_no_DR_filtered2.fasta'
out_path_snp = './data/snp_pos_no_DR_filtered2.csv'
path_to_depths = './data/coverages_with_percentiles/5/'
path_to_mummer = './data/mummer1.aligns'
path_to_annotations = './data/AL123456_rev.gff'
path_to_ref = './data/h37rv.fasta'
path_to_low_coverage_ids = './data/low_covered21.txt'

qual_thresold = 40
snp_cov_threshold = 15
sample_cov_threshold = 0.9

gene_upstream_threshold = 100

filter_by_upper_percentile = True
filter_DR_genes = True
print_alignment = True

# ref_len = 4411532

list_genes = ['ahpC', 'eis', 'embA', 'embB', 'embC', 'embR', 'fabG1', 'gid', 'gyrA', 'gyrB', 'inhA', 'iniA', 'iniC',
             'katG', 'manB', 'ndh', 'pncA', 'rmlD', 'rpoB', 'rpsA', 'rpsL', 'rrs', 'tlyA', 'ethR', 'fpbC', 'iniB',
              'kasA', 'ethA', 'fabD', 'efpA', 'thyA', 'panD', 'accD6', 'fbpC', 'nat', 'folC', 'rrl', 'rpoC', 'ribD', 'rplC']

def read_h37rv():
    fasta_sequences = SeqIO.parse(open(path_to_ref), 'fasta')
    sequence = ''
    for fasta in fasta_sequences:
        sequence = str(fasta.seq)
    return sequence.upper()


def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def parse_mummer():
    aln = []
    aln1 = []
    aln2 = []
    h37rv_from = 0
    h37rv_to = 0
    with open(path_to_mummer, 'r') as f:
        for line in f:
            if line.startswith('-- BEGIN'):
                if len(aln1) > 0:
                    aln.append((h37rv_from, h37rv_to, ''.join(aln1), ''.join(aln2)))
                    aln1.clear()
                    aln2.clear()
                s = line.split(' ')
                h37rv_from = int(s[5])
                h37rv_to = int(s[7])
            else:
                if line[0] not in ('\n', ' ', '-', '=', '/'):
                    s = line.split()
                    if len(aln1) == len(aln2):
                        aln1.append(s[1])
                    else:
                        aln2.append(s[1])
    aln.append((h37rv_from, h37rv_to, ''.join(aln1), ''.join(aln2)))
    aln.sort(key=lambda tup: tup[0])
    curr_segment = aln[0]
    aln_filtered = []
    for segment in aln:
        if segment[0] <= curr_segment[1]:
            # intersection
            if curr_segment[1] - curr_segment[0] < segment[1] - segment[0]:
                curr_segment = segment
        else:
            aln_filtered.append(curr_segment)
            curr_segment = segment
    aln_filtered.append(curr_segment)
    # with open('./test.txt', 'w') as f:
    #     for segment in aln_filtered:
    #         f.write(segment[2] + '\n')
    #         f.write(segment[3] + '\n')
    return aln_filtered


def preprocess_aln(aln):
    aln_new = []
    for segment in aln:
        cannetti_seq = []
        aln1 = segment[2]
        aln2 = segment[3]
        if len(aln1) != len(aln2):
            print('wrong lens: ' + str(segment[0]) + ' ' + str(segment[1]))
        start = 0
        is_gap = False
        for i in range(len(aln1)):
            if aln1[i] == '.':
                if not is_gap:
                    cannetti_seq.append(aln2[start:i])
                    is_gap = True
            else:
                if is_gap:
                    start = i
                    is_gap = False
        if start > 0:
            cannetti_seq.append(aln2[start:len(aln2)])
            aln_new.append((segment[0], segment[1], ''.join(cannetti_seq)))
        else:
            aln_new.append((segment[0], segment[1], aln2))
    # for segment in aln_new:
    #     if segment[1] - segment[0] != len(segment[2]) - 1:
    #         print('post proc err: ' + str(segment[0]) + ' ' + str(segment[1]))
    return aln_new


def h37rv_to_canetti(aln, pos):
    for segment in aln:
        if segment[0] <= pos <= segment[1]:
            if segment[2][pos - segment[0]] == '.':
                return '-'
            return segment[2][pos - segment[0]]
    return ''


def get_genes_coords():
    coords = []
    with open(path_to_annotations, 'r') as f:
        for line in f.readlines():
            s = line.split('\t')
            if s[2] == 'Gene':
                strand = s[6]
                for gene_name in list_genes:
                    if gene_name in s[-1]:
                        if strand == '+':
                            coords.append((int(s[3]) - gene_upstream_threshold, int(s[4])))
                        else:
                            coords.append((int(s[3]), int(s[4]) + gene_upstream_threshold))
            elif s[2] in ('Repeat_region', 'mobile_element'):
                coords.append((int(s[3]), int(s[4])))
    return coords


def get_repeats_coords():
    coords = []
    with open(path_to_annotations, 'r') as f:
        for line in f.readlines():
            s = line.split('\t')
            if s[2] in ('Repeat_region', 'mobile_element'):
                coords.append((int(s[3]), int(s[4])))
    return coords


def read_coverage(sample_id):
    coverage = []
    with open(path_to_depths + sample_id + '.coverage') as f:
        f.readline()
        percentile = float(f.readline().split('\t')[-1])
        f.readline()
        for line in f.readlines():
            s = line.split()
            coverage.append((int(s[0]), int(s[1])))
    coverage.sort(key=lambda tup: tup[0])
    return sample_id, coverage, percentile


def pos_covered(coverage, pos):
    for start, end in coverage:
        if start <= pos <= end:
            return True
    return False


def read_snps(sample_id, coords, gene_starts, percentile):
    snps = {}
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
                    cov = int(info[i + 3: j])
                    if cov > snp_cov_threshold:
                        if filter_by_upper_percentile and cov > percentile:
                            continue
                        snp_pos = int(tokens[1])
                        inside_genes = False
                        if filter_DR_genes:
                            i = bisect_left(gene_starts, snp_pos)
                            if i == 0:
                                if snp_pos == gene_starts[0]:
                                    inside_genes = True
                            else:
                                if snp_pos <= coords[i - 1][1]:
                                    inside_genes = True
                            # for coord in coords:
                            #     if coord[0] <= snp_pos <= coord[1]:
                            #         inside_genes = True
                            #         break
                        if not inside_genes:
                            snps[snp_pos] = tokens[4]
    return sample_id, snps


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


def main():
    sample_to_snps = {}
    all_snp_pos = set()
    sample_to_cov = {}

    sample_ids = [sample_id[:-1] for sample_id in open(path_to_ids, 'r').readlines()]
    sample_num = len(sample_ids)

    h37rv = read_h37rv()

    aln = parse_mummer()
    aln = preprocess_aln(aln)
    print('mummer parsed')

    if filter_DR_genes:
        filter_coords = get_genes_coords()
    else:
        filter_coords = get_repeats_coords()
    filter_coords.sort(key=lambda tup: tup[0])
    gene_starts = [x[0] for x in filter_coords]

    coverages = Parallel(n_jobs=-1)(delayed(read_coverage)(sample_id) for sample_id in sample_ids)
    sample_to_cov_list = {}
    percentiles = {}
    for sample_id, cov, percentile in coverages:
        sample_to_cov_list[sample_id] = cov
        percentiles[sample_id] = percentile
    print('read sample coverages')

    tasks = Parallel(n_jobs=-1)(delayed(read_snps)(sample_id, filter_coords, gene_starts, percentiles[sample_id]) for sample_id in sample_ids)
    for sample_id, snps in tasks:
        sample_to_snps[sample_id] = snps
        for snp_pos in snps.keys():
            all_snp_pos.add(snp_pos)
    all_snps = list(all_snp_pos)
    all_snps.sort()
    snp_num = len(all_snps)
    print('snp reading done')

    tasks = Parallel(n_jobs=-1)(delayed(sample_cov)(sample_id, all_snps, sample_to_cov_list[sample_id]) for sample_id in sample_ids)
    for sample_id, cov in tasks:
        sample_to_cov[sample_id] = cov
    # for sample_id in sample_ids:
    #     sample_id, cov = sample_cov(sample_id, all_snps)
    #     sample_to_cov[sample_id] = cov
    print('done with sample coverages')

    snp_cov_list = Parallel(n_jobs=-1)(
        delayed(compute_snp_cov)(all_snps, sample_to_cov_list[sample_id]) for sample_id in sample_ids)
    snp_cov = np.zeros(snp_num, dtype=int)
    for snp_cov_arr in snp_cov_list:
        snp_cov += snp_cov_arr
    print('done with snp coverages by samples')

    with open('./snp_cov.txt', 'w')as f:
        for i in range(snp_num):
            f.write(str(all_snps[i]) + ' ' + str(snp_cov[i]) + '\n')

    uncovered_snp_pos = set()
    for i in range(snp_num):
        if snp_cov[i]/sample_num < sample_cov_threshold:
            uncovered_snp_pos.add(all_snps[i])

    with open(out_path_snp, 'w') as f:
        for snp_pos in all_snps:
            if snp_pos not in uncovered_snp_pos:
                f.write(str(snp_pos) + '\n')
    print('printed all snps')

    with open(path_to_low_coverage_ids, 'w') as f:
        for sample_id in sample_ids:
            if sample_to_cov[sample_id] < sample_cov_threshold:
                f.write(sample_id + '\n')

    if print_alignment:
        with open(out_path_aln, 'w') as f:
            f.write('>H37Rv\n')
            for snp_pos in all_snps:
                if snp_pos not in uncovered_snp_pos:
                    f.write(h37rv[snp_pos - 1])
            f.write('\n')
            f.write('>canetti\n')
            for snp_pos in all_snps:
                if snp_pos not in uncovered_snp_pos:
                    n = h37rv_to_canetti(aln, snp_pos)
                    if n != '':
                        f.write(n.upper())
                    else:
                        f.write(h37rv[snp_pos - 1])
            f.write('\n')
            for sample_id, snps in sample_to_snps.items():
                if sample_to_cov[sample_id] >= sample_cov_threshold:
                    f.write('>' + sample_id + '\n')
                    for snp_pos in all_snps:
                        if snp_pos not in uncovered_snp_pos:
                            snp_letter = snps.get(snp_pos)
                            if snp_letter is not None:
                                f.write(snp_letter)
                            else:
                                f.write(h37rv[snp_pos - 1])
                    f.write('\n')


if __name__ == '__main__':
    main()
