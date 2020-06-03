from bisect import bisect_left
from os import listdir
from os.path import exists

from sklearn.externals.joblib import Parallel, delayed

from core.annotations import path_to_annotations
from core.constants import data_path, dr_genes, upstream_length
from core.data_reading import read_h37rv

path_to_snps = data_path + 'snps/combined_raw_variants_mq40_keep_complex_std_names_filtered/'
path_to_ids = data_path + 'combined_no_dup_snp_with_pheno.list'
# path_to_ids = data_path + 'casali14.list'
out_path_aln = data_path + 'combined_mq40_keep_complex_std_names_filtered_with_DR.fasta'
# out_path_aln = data_path + 'casali.fasta'
out_path_snps = data_path + 'combined_mq40_keep_complex_std_names_filtered_with_DR.snp_list'
# out_path_snps = data_path + 'casali.snp_list'
path_to_mummer = data_path + 'mummer1.aligns'

thread_num = 141

filter_out_DR_genes = False
filter_out_PGRS = False
filter_out_recombination_hotspots = False
path_to_recombination_hotspots = data_path + 'list_of_recombining_genes.txt'


def get_intervals_to_filter_out():
    """
    Reads annotations and picks coordinates of DR genes.

    :param filter_out_DR_genes: if True adds coordinates of DR genes from list of DR genes to the results
    :param filter_short_repeats: if True adds coordinates of short tandem repeats to the result
    :return: list of tuples (begin, end), representing intervals
    """
    recombining_genes = set(l.strip() for l in open(path_to_recombination_hotspots).readlines())
    coords = []
    with open(path_to_annotations, 'r') as f:
        for line in f.readlines():
            s = line.strip().split('\t')
            if s[2] == 'Gene':
                strand = s[6]
                if filter_out_DR_genes:
                    for gene_name in dr_genes:
                        if gene_name in s[-1]:
                            if strand == '+':
                                coords.append((int(s[3]) - upstream_length, int(s[4])))
                            else:
                                coords.append((int(s[3]), int(s[4]) + upstream_length))
                if filter_out_PGRS:
                    if 'PGRS' in s[-1] or 'PE' in s[-1]:
                        if strand == '+':
                            coords.append((int(s[3]) - upstream_length, int(s[4])))
                        else:
                            coords.append((int(s[3]), int(s[4]) + upstream_length))
                if filter_out_recombination_hotspots:
                    for gene_name in recombining_genes:
                        if gene_name in s[-1]:
                            if strand == '+':
                                coords.append((int(s[3]) - upstream_length, int(s[4])))
                            else:
                                coords.append((int(s[3]), int(s[4]) + upstream_length))
    coords.sort(key=lambda tup: tup[0])
    return coords


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
                    aln.append(
                        (h37rv_from, h37rv_to, ''.join(aln1), ''.join(aln2)))
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


def read_snps(sample_id, filter_intervals):
    # sample_id = fname[0: fname.rfind('.variants')]
    variants = {}
    filter_intervals_starts = [x[0] for x in filter_intervals]
    with open(path_to_snps + sample_id + '.variants') as f1:
        for line in f1.readlines():
            if line[0] == '#':
                continue
            s = line.strip().split('\t')
            if len(s[-1]) == 1 and len(s[-2]) == 1:
                pos = int(s[0])
                inside_filtered_interval = False
                i = bisect_left(filter_intervals_starts, pos)
                if i == 0:
                    if pos == filter_intervals_starts[0]:
                        inside_filtered_interval = True
                else:
                    if i < len(filter_intervals) and filter_intervals_starts[i] == pos:
                        inside_filtered_interval = True
                    elif pos <= filter_intervals[i - 1][1]:
                        inside_filtered_interval = True
                if not inside_filtered_interval:
                    variants[pos] = s[2]
    return sample_id, variants


def read_snps_no_filter(sample_id):
    # sample_id = fname[0: fname.rfind('.variants')]
    variants = {}
    with open(path_to_snps + sample_id + '.variants') as f1:
        for line in f1.readlines():
            if line[0] == '#':
                continue
            s = line.strip().split('\t')
            if len(s[-1]) == 1 and len(s[-2]) == 1:
                pos = int(s[0])
                variants[pos] = s[2]
    return sample_id, variants


def main():
    sample_to_snps = {}
    all_snp_pos = set()

    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    # fnames = [fname for fname in listdir(path_to_snps) if fname.endswith('.variants')]

    h37rv = read_h37rv()

    aln = parse_mummer()
    aln = preprocess_aln(aln)
    print('mummer parsed')

    if filter_out_DR_genes:
        filter_intervals = get_intervals_to_filter_out()
        tasks = Parallel(n_jobs=-1)(delayed(read_snps)(sample_id, filter_intervals) for sample_id in sample_ids
                                    if exists(path_to_snps + sample_id + '.variants'))
    else:
        tasks = Parallel(n_jobs=-1)(delayed(read_snps_no_filter)(sample_id) for sample_id in sample_ids
                                    if exists(path_to_snps + sample_id + '.variants'))


    for sample_id, snps in tasks:
        sample_to_snps[sample_id] = snps
        for pos, alt in snps.items():
            all_snp_pos.add(pos)
    all_snps = list(all_snp_pos)
    all_snps.sort()

    with open(out_path_snps, 'w') as f:
        for pos in all_snps:
            f.write(str(pos))
            f.write('\n')

    with open(out_path_aln, 'w') as f:
        f.write('>H37Rv\n')
        for snp_pos in all_snps:
            f.write(h37rv[snp_pos - 1])

        f.write('\n')
        f.write('>canetti\n')
        for snp_pos in all_snps:
            n = h37rv_to_canetti(aln, snp_pos)
            if n != '':
                f.write(n.upper())
            else:
                f.write(h37rv[snp_pos - 1])
        f.write('\n')
        for sample_id, snps in sample_to_snps.items():
            f.write('>' + sample_id + '\n')
            for snp_pos in all_snps:
                snp_letter = snps.get(snp_pos)
                if snp_letter is not None:
                    f.write(snp_letter)
                else:
                    f.write(h37rv[snp_pos - 1])
            f.write('\n')


if __name__ == '__main__':
    main()
