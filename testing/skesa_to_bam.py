from os import listdir, makedirs
from os.path import exists
from subprocess import check_call

import pysam
from Bio import SeqIO
from sklearn.externals.joblib import Parallel, delayed

from core.constants import data_path, ref_len


path_to_list = data_path + 'all_with_pheno.txt'

path_to_smalt = '/export/home/fedonin/smalt/bin/'
path_to_index = data_path + 'h37rv_rev'

# path_to_mummer_alns = data_path + 'mummer/dnadiff_bbduk_q30/'
path_to_mummer_alns = data_path + 'mummer/dnadiff/'
# path_to_skesa_contigs = data_path + 'skesa_bbduk_q30/'
path_to_skesa_contigs = data_path + 'skesa/'
# out_path = data_path + 'skesa_bam_all/'
# out_path = data_path + 'skesa_bam_no_bbduk/'
out_path = data_path + 'skesa_minimap/'

thread_num = '32'


def parse_delta(sample_id, path_to_mummer_alns):
    contig_to_mapping = {}
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
            aln_list = contig_to_mapping.get(cname)
            if aln_list is None:
                aln_list = [(ref_start, ref_end, contig_start, contig_end, cigar)]
                contig_to_mapping[cname] = aln_list
            else:
                aln_list.append((ref_start, ref_end, contig_start, contig_end, cigar))
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


def print_bam(sample_id, contigs, contig_to_mapping, out_path):
    header = {'HD': {'VN': '1.0'}, 'SQ': [{'LN': ref_len, 'SN': 'ch1'}]}

    alns = []
    for contig in contigs:
        if not contig.name in contig_to_mapping:
            continue
        for ref_start, ref_end, contig_start, contig_end, cigar in contig_to_mapping[contig.name]:
            a = pysam.AlignedSegment()
            a.query_name = contig.name
            a.is_paired = False
            if contig_start > contig_end:
                a.query_sequence = str(contig.seq[contig_end - 1: contig_start].reverse_complement())
                a.is_reverse = True
            else:
                a.query_sequence = str(contig.seq[contig_start - 1: contig_end])
            a.reference_id = 0
            a.reference_start = ref_start - 1
            a.mapping_quality = 20
            a.cigar = cigar
            alns.append(a)
    alns.sort(key=lambda x: x.reference_start)
    with pysam.AlignmentFile(out_path + sample_id + '.bam', "wb", header=header) as outf:
        for a in alns:
            outf.write(a)
    # chdir(out_path)
    check_call('samtools index ' + out_path + sample_id + '.bam', shell=True)


def parse_mummer(sample_id):
    contigs = [record for record in SeqIO.parse(path_to_skesa_contigs + sample_id + '.skesa.fa', 'fasta')]
    if len(contigs) == 0:
        return 0
    contig_to_mapping = parse_delta(sample_id, path_to_mummer_alns)
    print_bam(sample_id, contigs, contig_to_mapping, out_path)
    return 1


def smalt(id):
    if not exists(out_path + id + '.bam'):
        check_call(path_to_smalt + 'smalt map -c 0.9 -x -n ' + thread_num + ' ' + path_to_index + ' ' +
                   path_to_skesa_contigs + id + '.skesa.fa | samtools sort -@ ' + thread_num +
                   ' - | samtools view -Sb -@ ' + thread_num + ' - > ' +
                   out_path + id + '.bam', shell=True)
        check_call('samtools index -@ ' + thread_num + ' ' + out_path + id + '.bam', shell=True)


def bwa_mem(id):
    if not exists(out_path + id + '.bam'):
        check_call('bwa mem -t ' + thread_num + ' ' + path_to_index + ' ' + path_to_skesa_contigs + id + '.skesa.fa ' +
                   ' | samtools sort -@ ' + thread_num + ' - | samtools view -bS -@ ' + thread_num + ' - > ' + out_path +
                   id + '.bam', shell=True)
        check_call('samtools index -@ ' + thread_num + ' ' + out_path + id + '.bam', shell=True)


def bowtie(id):
    if not exists(out_path + id + '.bam'):
        check_call('bowtie2 --very-sensitive-local -f -a -p ' + thread_num + ' -x ' + path_to_index +
                   ' -U ' + path_to_skesa_contigs + id + '.skesa.fa | samtools sort -@ ' + thread_num +
                   ' - | samtools view -bS -@ ' + thread_num + ' - > ' + out_path + id + '.bam', shell=True)
        check_call('samtools index -@ ' + thread_num + ' ' + out_path + id + '.bam', shell=True)


def minimap(id):
    if not exists(out_path + id + '.bam'):
        check_call('/export/home/fedonin/minimap2-2.17_x64-linux/minimap2 -x asm5 -a -t ' + thread_num +
                    ' ' + data_path + 'h37rv.fasta ' + path_to_skesa_contigs + id +
                   '.skesa.fa | samtools sort -@ ' + thread_num + ' - | samtools view -bS -@ ' + thread_num + ' - > ' +
                   out_path + id + '.bam', shell=True)
        check_call('samtools index -@ ' + thread_num + ' ' + out_path + id + '.bam', shell=True)


if __name__ == '__main__':
    if not exists(out_path):
        makedirs(out_path)
    # sample_ids = [fname for fname in listdir(path_to_mummer_alns)]
    sample_ids = [l.strip() for l in open(path_to_list).readlines()]
    for sample_id in sample_ids:
        if not exists(path_to_skesa_contigs + sample_id + '.skesa.fa'):
            print('no contigs for ' + sample_id)
        else:
            print(sample_id)
            minimap(sample_id)
    # parse_mummer('SAMEA2535275')