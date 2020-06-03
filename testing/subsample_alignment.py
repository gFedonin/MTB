from Bio import SeqIO

from core.constants import data_path

path_to_aln = data_path + 'combined_mq40_keep_complex_std_names_filtered_no_DR.fasta'
path_to_list = data_path + 'casali14.list'
out_path = data_path + 'casali_sub.aln'


if __name__ == '__main__':
    sample_ids = set(l.strip() for l in open(path_to_list).readlines())
    seqs = [r for r in SeqIO.parse(path_to_aln, 'fasta') if r.id in sample_ids]
    SeqIO.write(seqs, out_path, 'fasta')