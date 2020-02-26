from subprocess import call

from core.constants import data_path

path_to_big_list = data_path + 'walker18/walker18_big.list'
path_to_filtered_list = data_path + 'walker18/walker18_new.samples'
out_path1 = data_path + 'walker18/walker18_rest.list'
out_path2 = data_path + 'walker18/walker18_triplets.list'

path_to_bams = data_path + 'walker18/bwa_mem_rev_rm_dup/'
out_path = data_path + 'walker18/merged_bams/'


def filter_list():
    filtered_set = set(l.strip() for l in open(path_to_filtered_list).readlines())
    rest_f = open(out_path1, 'w')
    triples_f = open(out_path2, 'w')
    for l in open(path_to_big_list).readlines():
        if ';' in l:
            s = l.strip().split(';')
            name_new = ''
            for n in s:
                if n in filtered_set:
                    name_new = n
                    triples_f.write(l)
            if name_new != '':
                for n in s:
                    if n != name_new:
                      rest_f.write(n + '\n')
    rest_f.close()
    triples_f.close()


def merge_fastq():
    filtered_set = set(l.strip() for l in open(path_to_filtered_list).readlines())
    for l in open(out_path2).readlines():
        s = l.strip().split(';')
        name_new = ''
        names = []
        for n in s:
            if n in filtered_set:
                name_new = n
            names.append(out_path1)
        call('cat ')


def merge_bams():
    filtered_set = set(l.strip() for l in open(path_to_filtered_list).readlines())
    for l in open(out_path2).readlines():
        print(l.strip())
        s = l.strip().split(';')
        name_new = ''
        names = []
        for n in s:
            if n in filtered_set:
                name_new = n
            names.append(path_to_bams + n + '.bam')
        call('samtools merge ' + out_path + name_new + '.bam ' + ' '.join(names), shell=True)


if __name__ == '__main__':
    merge_bams()
