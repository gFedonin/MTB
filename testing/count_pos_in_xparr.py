from os import listdir

from src.core.constants import data_path

path_to_XPARR = data_path + 'xparr/mc10_mega_MP/'
filter_path = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_no_win_qual_mqm_std3_mqm30_filter_samples_first/'
path_to_pos_list = filter_path + 'filtered_xparr_variant_pos.list'

out_path_syn = filter_path + 'xparr_counts_syn.csv'
out_path_miss = filter_path + 'xparr_counts_miss.csv'


def read_xparr(path, syn_counts, nonsyn_counts):
    for line in open(path).readlines()[2:]:
        for mut in line.strip().split('\t')[-1].split(';'):
            try:
                pos = int(mut)
            except:
                mut = mut[1:-1]
            c = nonsyn_counts.get(mut)
            if c is None:
                nonsyn_counts[mut] = 1
            else:
                nonsyn_counts[mut] = c + 1
        for mut in line.strip().split('\t')[-3].split(';'):
            c = syn_counts.get(mut)
            if c is None:
                syn_counts[mut] = 1
            else:
                syn_counts[mut] = c + 1


if __name__ == '__main__':
    nonsyn_counts = {}
    syn_counts = {}
    for file in listdir(path_to_XPARR):
        if file.endswith('.xparr'):
            read_xparr(path_to_XPARR + file, syn_counts, nonsyn_counts)
    with open(out_path_syn, 'w') as f_syn:
        with open(out_path_miss, 'w') as f_miss:
            for line in open(path_to_pos_list).readlines():
                pos = line.strip()
                c = syn_counts.get(pos)
                if c is None:
                    c = 0
                cn = nonsyn_counts.get(pos)
                if cn is None:
                    cn = 0
                if c > 0 and cn == 0:
                    f_syn.write(line)
                if c + cn == 0:
                    f_miss.write(line)