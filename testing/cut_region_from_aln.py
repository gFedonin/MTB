from Bio import SeqIO

path_to_aln = '../../data/codon_aln/Rifampicin_high.fasta'
coords = [8140, 8167]
sample_ids = {'H37Rv', 'SAMN03648890', 'SAMN03648892', 'SAMN03648885'}
out_path = '../../data/codon_aln_frac/Rifampicin_high.fasta'


if __name__ == '__main__':
    with open(out_path, 'w') as f:
        for r in SeqIO.parse(path_to_aln, 'fasta'):
            if r.id in sample_ids:
                f.write('>' + r.id + '\n')
                f.write(str(r.seq[coords[0]: coords[1] + 1]) + '\n')