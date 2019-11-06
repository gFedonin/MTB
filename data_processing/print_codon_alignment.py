from os import listdir, makedirs
from os.path import exists

from Bio import SeqIO
from sklearn.externals.joblib import Parallel, delayed

from core.annotations import CDSType, localize_all_variants, read_annotations
from core.constants import data_path, upstream_length
from core.data_reading import read_h37rv

use_list = True
# path_to_ids = data_path + 'all_with_pheno.txt'
path_to_ids = data_path + '10drugs.sample_list'
path_to_snps = data_path + 'snps/raw_with_DR_with_indel_with_pheno_and_snp_mc10/'
# out_path_aln = data_path + 'codon_aln_with_DR_with_indel_with_pheno_and_snp_mc10.fasta'
out_path_aln = data_path + 'codon_aln_with_DR_mc10_10drugs_rep.fasta'
# out_path_snp_to_genes = data_path + \
#     'codon_aln_with_DR_with_indel_with_pheno_and_snp_mc10.coords'
out_path_snp_to_genes = data_path + \
    'codon_aln_with_DR_mc10_10drugs_rep.coords'
out_path_codon_alns = data_path + 'codon_aln_rep/'

# path_to_pheno_and_trees = data_path + 'reconstructed_mc10_mega_MP/'
path_to_high_mut_genes = '../../res/tree_was/mc10_mega_MP_RR_filter/upper.filtered_fdr0.1_pvalue0.05.genes.csv'
path_to_low_mut_genes = '../../res/tree_was/mc10_mega_MP_RR_filter/lower.filtered_fdr0.1_pvalue0.05.genes.csv'


drug_names = ('Isoniazid', 'Rifampicin', 'Ethambutol', 'Pyrazinamide', 'Streptomycin', 'Moxifloxacin', 'Ofloxacin',
                  'Amikacin', 'Capreomycin', 'Kanamycin')


def read_snps(sample_id):
    snps = {}
    with open(path_to_snps + sample_id + '.variants', 'r') as f1:
        lines = f1.readlines()
        for line in lines:
            s = line.strip().split('\t')
            if s[-1] == 'snp':
                snps[int(s[0])] = s[1]
    return sample_id, snps


def print_aln_for_all_samples():
    sample_to_snps = {}
    all_snp_pos = set()

    if use_list:
        sample_ids = [sample_id[:-1]
                    for sample_id in open(path_to_ids, 'r').readlines()]
    else:
        sample_ids = [fname[:-len('.variants')] for fname in listdir(path_to_snps) if fname.endswith('.variants')]

    cds = read_annotations(upstream_length)

    h37rv = read_h37rv()

    # aln = parse_mummer()
    # aln = preprocess_aln(aln)
    # print('mummer parsed')

    tasks = Parallel(n_jobs=-1)(delayed(read_snps)(sample_id)
                                for sample_id in sample_ids)
    for sample_id, snps in tasks:
        sample_to_snps[sample_id] = snps
        for snp_pos in snps.keys():
            all_snp_pos.add(snp_pos)
    all_snps = list(all_snp_pos)

    snp_to_cds = localize_all_variants(all_snps, cds)

    snps_in_genes = []
    pos_in_triplet = []
    last_pos = -1
    for snp_pos in all_snps:
        cds = snp_to_cds.get(snp_pos)
        if cds is not None and cds.type == CDSType.Gene:
            prot_pos = (snp_pos - cds.start)//3
            if prot_pos != last_pos:
                snps_in_genes.append(snp_pos)
                pos_in_triplet.append((snp_pos - cds.start) % 3)
                last_pos = prot_pos

    with open(out_path_snp_to_genes, 'w') as f:
        for snp_pos in snps_in_genes:
            cds = snp_to_cds[snp_pos]
            f.write('%s\t%d\t%d\t%d\n' % (cds.name, cds.start, cds.end, snp_pos))

    with open(out_path_aln, 'w') as f:
        f.write('>H37Rv\n')
        i = 0
        for snp_pos in snps_in_genes:
            triplet_pos = pos_in_triplet[i]
            f.write(h37rv[snp_pos - triplet_pos - 1: snp_pos - triplet_pos + 2])
            i += 1
        f.write('\n')
        # f.write('>canetti\n')
        # i = 0
        # for snp_pos in snps_in_genes:
        #     triplet_pos = pos_in_triplet[i]
        #     triplet = []
        #     for j in range(3):
        #         n = h37rv_to_canetti(aln, snp_pos - triplet_pos + j)
        #         if n != '':
        #             triplet.append(n.upper())
        #         else:
        #             triplet.append(h37rv[snp_pos - triplet_pos + j - 1])
        #     f.write(''.join(triplet))
        #     i += 1
        # f.write('\n')
        for sample_id, snps in sample_to_snps.items():
            f.write('>' + sample_id + '\n')
            i = 0
            for snp_pos in snps_in_genes:
                triplet_pos = pos_in_triplet[i]
                triplet = []
                for j in range(3):                
                    snp_letter = snps.get(snp_pos - triplet_pos + j)
                    if snp_letter is not None:
                        triplet.append(snp_letter)
                    else:
                        triplet.append(h37rv[snp_pos - triplet_pos + j - 1])
                f.write(''.join(triplet))
                i += 1
            f.write('\n')


def print_aln_for_drug(drug, samples, genes, sample_to_snps, snp_to_cds, h37rv, suffix):
    all_snp_pos = set()
    for sample_id in samples:
        snps = sample_to_snps[sample_id]
        for snp_pos in snps.keys():
            all_snp_pos.add(snp_pos)
    all_snps = list(all_snp_pos)
    all_snps.sort()

    snps_in_genes = []
    pos_in_triplet = []
    last_pos = -1
    for snp_pos in all_snps:
        cds = snp_to_cds.get(snp_pos)
        if cds is not None and cds.type == CDSType.Gene and cds.name in genes:
            prot_pos = (snp_pos - cds.start)//3
            if prot_pos != last_pos:
                snps_in_genes.append(snp_pos)
                pos_in_triplet.append((snp_pos - cds.start) % 3)
                last_pos = prot_pos

    with open(out_path_codon_alns + drug + suffix + '.coords', 'w') as f:
        for snp_pos in snps_in_genes:
            cds = snp_to_cds[snp_pos]
            f.write('%s\t%d\t%d\t%d\n' % (cds.name, cds.start, cds.end, snp_pos))

    with open(out_path_codon_alns + drug + suffix + '.fasta', 'w') as f:
        f.write('>H37Rv\n')
        i = 0
        for snp_pos in snps_in_genes:
            triplet_pos = pos_in_triplet[i]
            f.write(h37rv[snp_pos - triplet_pos - 1: snp_pos - triplet_pos + 2])
            i += 1
        f.write('\n')
        # f.write('>canetti\n')
        # i = 0
        # for snp_pos in snps_in_genes:
        #     triplet_pos = pos_in_triplet[i]
        #     triplet = []
        #     for j in range(3):
        #         n = h37rv_to_canetti(aln, snp_pos - triplet_pos + j)
        #         if n != '':
        #             triplet.append(n.upper())
        #         else:
        #             triplet.append(h37rv[snp_pos - triplet_pos + j - 1])
        #     f.write(''.join(triplet))
        #     i += 1
        # f.write('\n')
        for sample_id, snps in sample_to_snps.items():
            f.write('>' + sample_id + '\n')
            # if sample_id == 'SAMN03648890':
            #     a = 0
            i = 0
            for snp_pos in snps_in_genes:
                # if 1282075 < snp_pos < 1282085:
                #     a = 0
                triplet_pos = pos_in_triplet[i]
                triplet = []
                for j in range(3):
                    snp_letter = snps.get(snp_pos - triplet_pos + j)
                    if snp_letter is not None:
                        triplet.append(snp_letter)
                    else:
                        triplet.append(h37rv[snp_pos - triplet_pos + j - 1])
                f.write(''.join(triplet))
                i += 1
            f.write('\n')


def print_aln_for_all_drugs():
    if not exists(out_path_codon_alns):
        makedirs(out_path_codon_alns)
    if use_list:
        sample_ids = [sample_id[:-1]
                    for sample_id in open(path_to_ids, 'r').readlines()]
    else:
        sample_ids = [fname[:-len('.variants')] for fname in listdir(path_to_snps) if fname.endswith('.variants')]

    cds = read_annotations(upstream_length)

    h37rv = read_h37rv()

    # aln = parse_mummer()
    # aln = preprocess_aln(aln)
    # print('mummer parsed')

    sample_to_snps = {}

    tasks = Parallel(n_jobs=-1)(delayed(read_snps)(sample_id)
                                for sample_id in sample_ids)
    all_snp_pos = set()
    for sample_id, snps in tasks:
        sample_to_snps[sample_id] = snps
        for snp_pos in snps.keys():
            all_snp_pos.add(snp_pos)
    all_snps = list(all_snp_pos)

    snp_to_cds = localize_all_variants(all_snps, cds)

    drug_to_gene_set = {drug: set() for drug in drug_names}
    for l in open(path_to_high_mut_genes).readlines():
        s = l.strip().split('\t')
        drug_to_gene_set[s[0]].add(s[1])
    for drug in drug_names:
        print_aln_for_drug(drug, sample_ids, drug_to_gene_set[drug], sample_to_snps, snp_to_cds, h37rv, '_high')
    drug_to_gene_set = {drug: set() for drug in drug_names}
    for l in open(path_to_low_mut_genes).readlines():
        s = l.strip().split('\t')
        drug_to_gene_set[s[0]].add(s[1])
    for drug in drug_names:
        print_aln_for_drug(drug, sample_ids, drug_to_gene_set[drug], sample_to_snps, snp_to_cds, h37rv, '_low')


if __name__ == '__main__':
    # print_aln_for_all_samples()
    print_aln_for_all_drugs()