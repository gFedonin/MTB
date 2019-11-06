path_to_high_mut_snp_pos = '../res/tree_was/mc10_mega_MP_RR_filter/pairs.upper.filtered_fdr0.1_pvalue0.05.converted.filtered_by_mqm.csv'
path_to_low_mut_snp_pos = '../res/tree_was/mc10_mega_MP_RR_filter/pairs.lower.filtered_fdr0.1_pvalue0.05.converted.filtered_by_mqm.csv'

out_path_high_mut_genes = '../res/tree_was/mc10_mega_MP_RR_filter/upper.filtered_fdr0.1_pvalue0.05.genes.csv'
out_path_low_mut_genes = '../res/tree_was/mc10_mega_MP_RR_filter/lower.filtered_fdr0.1_pvalue0.05.genes.csv'

drug_names = ('Isoniazid', 'Rifampicin', 'Ethambutol', 'Pyrazinamide', 'Streptomycin', 'Moxifloxacin', 'Ofloxacin',
                  'Amikacin', 'Capreomycin', 'Kanamycin')


def extract_mutating_genes():
    drug_to_genes = {drug: set() for drug in drug_names}
    for l in open(path_to_high_mut_snp_pos).readlines()[1:]:
        s = l.strip().split('\t')
        if s[4] in drug_names:
            drug_to_genes[s[4]].add(s[1])
    with open(out_path_high_mut_genes, 'w') as f:
        for drug, gene_set in drug_to_genes.items():
            if len(gene_set) > 0:
                gene_list = list(gene_set)
                gene_list.sort()
                for gene in gene_list:
                    f.write('%s\t%s\n' % (drug, gene))

    drug_to_genes = {drug: set() for drug in drug_names}
    for l in open(path_to_low_mut_snp_pos).readlines()[1:]:
        s = l.strip().split('\t')
        if s[4] in drug_names:
            drug_to_genes[s[4]].add(s[1])
    with open(out_path_low_mut_genes, 'w') as f:
        for drug, gene_set in drug_to_genes.items():
            if len(gene_set) > 0:
                gene_list = list(gene_set)
                gene_list.sort()
                for gene in gene_list:
                    f.write('%s\t%s\n' % (drug, gene))
                      

if __name__ == "__main__":
    extract_mutating_genes()