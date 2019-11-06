from Bio import SeqIO
from sklearn.externals.joblib import Parallel, delayed

from core.annotations import CDSType, localize_all_variants, read_annotations
from core.constants import upstream_length
from core.data_reading import read_h37rv
from data_processing.print_alignment import parse_mummer, preprocess_aln, h37rv_to_canetti

path_to_ids = './data/Full_subset_filtered_snp_pheno2.txt'
path_to_snps = '/data/snps/raw_with_DR/'
out_path_aln = './data/snp_aln_genes_filtered2.fasta'
path_to_mummer = './data/mummer1.aligns'
path_to_ref = './data/h37rv.fasta'
path_to_annotations = './data/AL123456_rev.gff'


def read_snps(sample_id):
    snps = {}
    with open(path_to_snps + sample_id + '.snp', 'r') as f1:
        lines = f1.readlines()
        for line in lines:
            s = line[:-1].split('\t')
            snps[int(s[0])] = s[1]
    return sample_id, snps


def main():
    sample_to_snps = {}
    all_snp_pos = set()

    sample_ids = [sample_id[:-1] for sample_id in open(path_to_ids, 'r').readlines()]

    cds = read_annotations(path_to_annotations, upstream_length)

    h37rv = read_h37rv()

    aln = parse_mummer()
    aln = preprocess_aln(aln)
    print('mummer parsed')

    tasks = Parallel(n_jobs=-1)(delayed(read_snps)(sample_id) for sample_id in sample_ids)
    for sample_id, snps in tasks:
        sample_to_snps[sample_id] = snps
        for snp_pos in snps.keys():
            all_snp_pos.add(snp_pos)
    all_snps = list(all_snp_pos)

    snp_to_cds = localize_all_variants(all_snps, cds)

    snps_in_genes = []
    for snp in all_snps:
        cds = snp_to_cds[snp]
        if cds is not None and cds.type == CDSType.Gene:
            snps_in_genes.append(snp)

    with open(out_path_aln, 'w') as f:
        f.write('>H37Rv\n')
        for snp_pos in snps_in_genes:
            f.write(h37rv[snp_pos - 1])
        f.write('\n')
        f.write('>canetti\n')
        for snp_pos in snps_in_genes:
            n = h37rv_to_canetti(aln, snp_pos)
            if n != '':
                f.write(n.upper())
            else:
                f.write(h37rv[snp_pos - 1])
        f.write('\n')
        for sample_id, snps in sample_to_snps.items():
            f.write('>' + sample_id + '\n')
            for snp_pos in snps_in_genes:
                snp_letter = snps.get(snp_pos)
                if snp_letter is not None:
                    f.write(snp_letter)
                else:
                    f.write(h37rv[snp_pos - 1])
            f.write('\n')


if __name__ == '__main__':
    main()
