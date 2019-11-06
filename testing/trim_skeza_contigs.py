from os import makedirs, listdir
from os.path import exists
from subprocess import check_call

from Bio import SeqIO
import pysam
from core.constants import data_path, ref_len

path_to_mummer_alns = data_path + 'mummer/dnadiff/'
path_to_skesa_contigs = data_path + 'skesa/'
out_path = data_path + 'skesa_trimmed_bam/'


class ContigAln:

    def __init__(self, cname, ref_start, ref_end, contig_start, contig_end, cigar):
        self.cname = cname
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.contig_start = contig_start
        self.contig_end = contig_end
        self.cigar = cigar


def parse_delta(sample_id, path_to_mummer_alns):
    contig_to_mapping = []
    f = open(path_to_mummer_alns + sample_id + '/' + sample_id + '.1delta')
    lines = f.readlines()
    f.close()
    cname = lines[2].split(' ')[1]
    i = 3
    while i < len(lines):
        l = lines[i]
        if l[0] == '>':
            s = l.strip().split(' ')
            cname = s[1]
            # if cname == 'Contig_188_112.007':
            #     a = 0
            i += 1
        else:
            s = l.split(' ')
            ref_start = int(s[0])
            ref_end = int(s[1])
            contig_start = int(s[2])
            contig_end = int(s[3])
            contig_aln_len = abs(contig_end - contig_start) + 1
            cigar = []
            ins_len = 0
            mandi_len = 0
            del_len = 0
            i += 1
            while(True):
                l = lines[i].strip()
                if l == '0':
                    if ins_len > 0:
                        # there was ins before
                        # add previous ins with len ins_len
                        cigar.append((1, ins_len))
                        mandi_len += ins_len
                    if del_len > 0:
                        # there was del before
                        # add previous del
                        cigar.append((2, del_len))
                    # add match with len contig_end - contig_start - indel_end
                    if contig_aln_len - mandi_len < 0:
                        print(cname + ' AAAAAAAAA!')
                        cigar_str = ''
                        for op, l in cigar:
                            if op == 0:
                                cigar_str += str(l) + 'M'
                            elif op == 1:
                                cigar_str += str(l) + 'I'
                            else:
                                cigar_str += str(l) + 'D'
                        print(cigar_str)
                    if contig_aln_len - mandi_len > 0:
                        cigar.append((0, contig_aln_len - mandi_len))
                    i += 1
                    break
                d = int(l)
                if d < 0:
                    d = -d
                    # insertion
                    if ins_len > 0:
                        # there was ins before
                        if d == 1:
                            # same ins
                            ins_len += 1
                            i += 1
                            continue
                        else:
                            # new ins
                            # add previous ins with len ins_len
                            cigar.append((1, ins_len))
                            mandi_len += ins_len
                    if del_len > 0:
                        # there was del before
                        # add previous del
                        cigar.append((2, del_len))
                        del_len = 0
                    if d > 1:
                        # add match with len d - 1
                        cigar.append((0, d - 1))
                        mandi_len += d - 1
                    # new ins
                    ins_len = 1
                else:
                    # deletion
                    if ins_len > 0:
                        # there was ins before
                        # add previous ins
                        cigar.append((1, ins_len))
                        mandi_len += ins_len
                        ins_len = 0
                    if del_len > 0:
                        # there was del before
                        if d == 1:
                            # same del
                            del_len += 1
                            # indel_end += 1
                            i += 1
                            continue
                        else:
                            # new del
                            # add previous del with len del_len
                            cigar.append((2, del_len))
                    if d > 1:
                        # add match with len d - 1
                        cigar.append((0, d - 1))
                        mandi_len += d - 1
                    # new del
                    del_len = 1
                i += 1
            contig_to_mapping.append(ContigAln(cname, ref_start, ref_end, contig_start, contig_end, cigar))
            cigar_len = 0
            for op, l in cigar:
                if op != 2:
                    cigar_len += l
            if contig_aln_len != cigar_len:
                print('%s %s %d %d' % (sample_id, cname, contig_aln_len, cigar_len))
                cigar_str = ''
                c_len = 0
                for op, l in cigar:
                    if op == 0:
                        cigar_str += str(l) + 'M'
                        c_len += l
                    elif op == 1:
                        cigar_str += str(l) + 'I'
                        c_len += l
                    else:
                        cigar_str += str(l) + 'D'
                print(cigar_str)
                print(str(c_len))
                exit(1)
    return contig_to_mapping


def trim_alns(contigs, contig_to_mapping):
    res = []
    contig_to_mapping.sort(key=lambda x: x[1])
    prev = contig_to_mapping.pop(0)
    while len(contig_to_mapping) > 0:
        contig_aln = contig_to_mapping.pop(0)
        if contig_aln.ref_start < prev.ref_end:
            if contig_aln.ref_end <= prev.ref_end:
                # this is a small contig inside the big one
                continue
            move_forward = True
            if prev.ref_end - prev.ref_start > contig_aln.ref_end - contig_aln.ref_start:
                # keep prev, trim current
                cigar = contig_aln.cigar
                if contig_aln.contig_start > contig_aln.contig_end:
                    a = 0
                else:
                    contig_pos = contig_aln.contig_start
                    ref_pos = contig_aln.ref_start
                    while ref_pos < prev.ref_end:
                        op, l = cigar.pop(0)
                        if op == 0:
                            # match
                            if ref_pos + l < prev.ref_end:
                                contig_pos += l
                                ref_pos += l
                            else:
                                contig_pos += prev.ref_end - ref_pos
                                ref_pos = prev.ref_end
                                cigar.insert(0, (0, l - prev.ref_end + ref_pos))
                        elif op == 1:
                            # insert
                            contig_pos += l
                        else:
                            # deletion
                            if ref_pos + l < prev.ref_end:
                                ref_pos += l
                            else:
                                # deletion end after prev contig end:
                                # we need to add current contig in the sorted list again
                                ref_pos = ref_pos + l
                                i = 0
                                next_aln = contig_to_mapping[i]
                                while next_aln.ref_start < ref_pos:
                                    # curr contig is located after the next one in the list
                                    i += 1
                                    next_aln = contig_to_mapping[i]
                                contig_to_mapping.insert(i, contig_aln)
                                move_forward = False
                    contig_aln.contig_start = contig_pos
                    contig_aln.ref_start = ref_pos
                if move_forward:
                    res.append(prev)
                    prev = contig_aln
            else:
                # trim prev
                a = 0


def print_bam(sample_id, contigs, contig_to_mapping, out_path):
    header = {'HD': {'VN': '1.0'}, 'SQ': [{'LN': ref_len, 'SN': 'ch1'}]}

    alns = []
    for contig_aln in contig_to_mapping:
        contig = contigs[contig_aln.cname]
        a = pysam.AlignedSegment()
        a.query_name = contig_aln.cname
        a.is_paired = False
        if contig_aln.contig_start > contig_aln.contig_end:
            a.query_sequence = str(contig.seq[contig_aln.contig_end - 1: contig_aln.contig_start].reverse_complement())
            a.is_reverse = True
        else:
            a.query_sequence = str(contig.seq[contig_aln.contig_start - 1: contig_aln.contig_end])
        a.reference_id = 0
        a.reference_start = contig_aln.ref_start - 1
        a.mapping_quality = 20
        a.cigar = contig_aln.cigar
        alns.append(a)
    alns.sort(key=lambda x: x.reference_start)
    with pysam.AlignmentFile(out_path + sample_id + '.bam', "wb", header=header) as outf:
        for a in alns:
            outf.write(a)
    # chdir(out_path)
    check_call('samtools index ' + out_path + sample_id + '.bam', shell=True)


def parse_mummer(sample_id):
    contigs = {record.id: record for record in SeqIO.parse(path_to_skesa_contigs + sample_id + '.skesa.fa', 'fasta')}
    if len(contigs) == 0:
        return 0
    contig_to_mapping = parse_delta(sample_id, path_to_mummer_alns)
    print_bam(sample_id, contigs, contig_to_mapping, out_path)
    return 1


if __name__ == '__main__':
    if not exists(out_path):
        makedirs(out_path)
    sample_ids = [fname for fname in listdir(path_to_mummer_alns)]
    # tasks = Parallel(n_jobs=-1, batch_size=len(sample_ids)//thread_num + 1)(delayed(parse_mummer)(sample_id)
    #                                                                         for sample_id in sample_ids)
    # c = 0
    # for task in tasks:
    #     c += task
    # print('%d samples processed' % c)
    for sample_id in sample_ids:
        print(sample_id)
        parse_mummer(sample_id)
    parse_mummer('SAMEA2535275')