import gzip
from os import makedirs
from os.path import exists
from subprocess import check_call, Popen, PIPE, call

from core.constants import data_path

# path_to_ids = data_path + 'all_with_pheno.txt'
path_to_ids = data_path + 'list3.txt'
# path_to_ids = data_path + 'debug1.list'
path_to_reads = '/export/data/kkuleshov/myc/sra/'
out_path = data_path + 'storm/'
thread_num = '32'
max_iter_num = 4



def interleave(sample_id):
    f1 = gzip.open(path_to_reads + sample_id + '/' + sample_id + '_R1.fastq.gz', 'rt')
    f2 = gzip.open(path_to_reads + sample_id + '/' + sample_id + '_R2.fastq.gz', 'rt')
    proc = Popen('~/seqtk seq -a - > ' + out_path + sample_id + '/' + sample_id + '.fasta', stdin=PIPE, shell=True)
    while True:
        line = f1.readline()
        if line == "":
            break
        proc.stdin.write(line.encode())
        for i in range(3):
            proc.stdin.write(f1.readline().encode())
        for i in range(4):
            proc.stdin.write(f2.readline().encode())
    outs, errs = proc.communicate()


if __name__ == '__main__':
    sample_ids = [l.strip() for l in open(path_to_ids).readlines()]
    for sample_id in sample_ids:
        if not exists(out_path + sample_id):
            makedirs(out_path + sample_id)
        interleave(sample_id)
        call('~/RACER ' + out_path + sample_id + '/' + sample_id + '.fasta ' + out_path + sample_id + '/' + sample_id +
             '_racer.fasta 4000000', shell=True)
        query = out_path + sample_id + '/' + sample_id + '_racer.fasta'
        subject = out_path + sample_id + '/' + sample_id + '_racer.fasta'
        i = 1
        while(True):
            call('~/storm -i MergePairedEndReads --query ' + query + ' --orient 2 --subject ' + subject +
                       ' --out ' + out_path + sample_id + '/racer_iter' + str(i) + '_ -m 0 -ll 40 -rl 40 -k 29 -t ' +
                       thread_num, shell=True)
            query = out_path + sample_id + '/racer_iter' + str(i) + '_notmerged.fasta'
            subject = out_path + sample_id + '/racer_iter' + str(i) + '_merged.fasta'
            l = open(subject).readline()
            if l == '' or i == max_iter_num:
                break
            i += 1

        call('cat ' + out_path + sample_id + '/*merged.fasta > ' + out_path  + sample_id + '/' + sample_id +
                   '_merged_racer.fasta', shell=True)
        call('mv ' + out_path + sample_id + '/racer_iter' + str(i) + '_notmerged.fasta ' + out_path + sample_id + '/' +
                   sample_id + '_notmerged_racer.fasta', shell=True)
        call('rm ' + out_path + sample_id + '/racer_iter*', shell=True)
        # call('rm ' + out_path + sample_id + '/' + sample_id + '.fasta', shell=True)
        # call('~/karect -correct -threads=' + thread_num + ' -matchtype=hamming -celltype=haploid -minoverlap=90 ' +
        #      '-minoverlapper=0 -inputfile=' + sample_id + '_merged.fasta -inputfile=' + sample_id + '_notmerged.fasta', shell=True,
        #      cwd=out_path + sample_id)
        # call('rm ' + out_path + sample_id + '/res_graph*', shell=True)
        # call('rm ' + out_path + sample_id + '/temp_*', shell=True)
        # call('rm ' + out_path + sample_id + '/input_file.txt', shell=True)
        # query = out_path + sample_id + '/karect_' + sample_id + '_notmerged.fasta'
        # subject = out_path + sample_id + '/karect_' + sample_id + '_merged.fasta'
        # i = 1
        # while(True):
        #     call('~/storm -i MergePairedEndReads --query ' + query + ' --orient 2 --subject ' + subject +
        #                ' --out ' + out_path + sample_id + '/karect_racer_iter' + str(i) + '_ -m 0 -ll 40 -rl 40 -k 29 -t ' +
        #                thread_num, shell=True)
        #     query = out_path + sample_id + '/karect_racer_iter' + str(i) + '_notmerged.fasta'
        #     subject = out_path + sample_id + '/karect_racer_iter' + str(i) + '_merged.fasta'
        #     l = open(subject).readline()
        #     if l == '':
        #         break
        #     i += 1
        # call('cat ' + out_path + sample_id + '/karect_racer_iter*merged.fasta > ' + out_path  + sample_id + '/karect_' + sample_id +
        #            '_merged.fasta', shell=True)
        # call('mv ' + out_path + sample_id + '/karect_racer_iter' + str(i) + '_notmerged.fasta ' + out_path + sample_id + '/karect_' +
        #            sample_id + '_notmerged.fasta', shell=True)
        # call('rm ' + out_path + sample_id + '/karect_racer_iter*', shell=True)