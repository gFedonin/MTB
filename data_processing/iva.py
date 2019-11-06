from os import makedirs, listdir
from os.path import exists
from subprocess import check_call

import pysam
from Bio import SeqIO
from sklearn.externals.joblib import Parallel, delayed

from core.constants import data_path, ref_len
from core.data_reading import path_to_ref

path_to_mummer = '/export/home/fedonin/mummer4/bin/'
path_to_samtools = '/export/home/fedonin/samtools-1.2/'

path_to_samples = '/export/data/kkuleshov/myc/sra/'
path_to_list = data_path + 'debug2.list'
# path_to_list = data_path + 'all_with_pheno.txt'

iva_out_path = data_path + 'iva/'
mummer_out_path = data_path + 'mummer/dnadiff_iva/'
mummer_log_path = data_path + 'mummer/log_iva/'
mummer_out_path_snps = data_path + 'mummer/snps_iva/'
final_raw_snps_path = data_path + 'snps/iva_mummer_raw_ld_mum4/'

out_path_bam = data_path + 'iva_bam/'

thread_num = '32'


def run_iva_assembly():
    if not exists(iva_out_path):
        makedirs(iva_out_path)
    sample_ids = [l.strip() for l in open(path_to_list).readlines()]
    for sample_id in sample_ids:
        if exists(iva_out_path + '/' + sample_id + '/contigs.fasta'):
            continue
        if not exists(iva_out_path + sample_id):
            makedirs(iva_out_path + sample_id)
        check_call('iva --threads ' + thread_num + ' -f ' + path_to_samples + sample_id + '/' + sample_id +
                   '_R1.fastq.gz -r ' + path_to_samples + sample_id + '/' + sample_id + '_R2.fastq.gz ' +
                   iva_out_path + sample_id + '/', shell=True)


def run_mummer(sample_id):
    if not exists(mummer_out_path + sample_id):
        makedirs(mummer_out_path + sample_id)
    if not exists(mummer_out_path + sample_id + '/' + sample_id + '.snps'):
        contigs = [l for l in open(iva_out_path + '/' + sample_id + '/assembly.fasta').readlines()]
        if len(contigs) == 0:
            print('no contigs for ' + sample_id)
            return 1
        status = check_call(path_to_mummer + 'dnadiff  --prefix=' + mummer_out_path + sample_id + '/' + sample_id + ' ' +
                            path_to_ref + ' ' + iva_out_path + '/' + sample_id + '/assembly.fasta > ' + mummer_log_path +
                            sample_id + '.dnadiff.log 2>> ' + mummer_log_path + sample_id + '.dnadiff.log', shell=True)
        if status != 0:
            print('mummer problem with ' + sample_id)
        return 1
    return 0


def mass_mummer():
    if not exists(mummer_out_path):
        makedirs(mummer_out_path)
    if not exists(mummer_log_path):
        makedirs(mummer_log_path)
    # if not exists(mummer_out_path_filtered):
    #     makedirs(mummer_out_path_filtered)
    if not exists(mummer_out_path_snps):
        makedirs(mummer_out_path_snps)
    tasks = Parallel(n_jobs=-1)(delayed(run_mummer)(l) for l in listdir(iva_out_path))
    c = 0
    for task in tasks:
        c += task
    print('%d samples processed' % c)


def parse_mummer_snps(sample_id):
    if not exists(mummer_out_path + sample_id + '/' + sample_id + '.snps'):
        print('no snps for ' + sample_id)
        return 1
    with open(final_raw_snps_path + sample_id + '.variants', 'w') as f:
        is_ins = False
        ins_str = []
        ins_pos = None
        is_del = False
        del_len = 0
        del_pos = None
        for l in open(mummer_out_path + sample_id + '/' + sample_id + '.snps').readlines()[5:]:
            s = l.strip().split()
            if s[1] == '.':
                if is_del:
                    f.write('%s\t%d\tdel\n' % (del_pos, del_len))
                    is_del = False
                    del_len = 0
                    del_pos = None
                is_ins = True
                ins_pos = s[0]
                ins_str.append(s[2])
            elif s[2] == '.':
                if is_ins:
                    f.write('%s\t%s\tins\n' % (ins_pos, ''.join(ins_str)))
                    is_ins = False
                    ins_str.clear()
                    ins_pos = None
                is_del = True
                del_pos = s[0]
                del_len += 1
            else:
                if is_del:
                    f.write('%s\t%d\tdel\n' % (del_pos, del_len))
                    is_del = False
                    del_len = 0
                    del_pos = None
                elif is_ins:
                    f.write('%s\t%s\tins\n' % (ins_pos, ''.join(ins_str)))
                    is_ins = False
                    ins_str.clear()
                    ins_pos = None
                f.write('%s\t%s\tsnp\n' % (s[0], s[2]))
    return 1


def parse_all_mummer_snps():
    if not exists(final_raw_snps_path):
        makedirs(final_raw_snps_path)
    tasks = Parallel(n_jobs=-1)(delayed(parse_mummer_snps)(l) for l in listdir(mummer_out_path))
    c = 0
    for task in tasks:
        c += task
    print('%d samples processed' % c)


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
    check_call(path_to_samtools + 'samtools index ' + out_path + sample_id + '.bam', shell=True)


def parse_mummer(sample_id):
    contigs = [record for record in SeqIO.parse(iva_out_path + '/' + sample_id + '/assembly.fasta', 'fasta')]
    if len(contigs) == 0:
        return 0
    contig_to_mapping = parse_delta(sample_id, mummer_out_path)
    print_bam(sample_id, contigs, contig_to_mapping, out_path_bam)
    return 1


def delta_to_bam():
    if not exists(out_path_bam):
        makedirs(out_path_bam)
    sample_ids = [fname for fname in listdir(mummer_out_path)]
    # tasks = Parallel(n_jobs=-1, batch_size=len(sample_ids)//thread_num + 1)(delayed(parse_mummer)(sample_id)
    #                                                                         for sample_id in sample_ids)
    # c = 0
    # for task in tasks:
    #     c += task
    # print('%d samples processed' % c)
    for sample_id in sample_ids:
        print(sample_id)
        parse_mummer(sample_id)


if __name__ == '__main__':
    run_iva_assembly()
    mass_mummer()
    parse_all_mummer_snps()
    delta_to_bam()