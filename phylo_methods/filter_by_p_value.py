from os.path import exists

from src.core.constants import ref_len
from src.phylo_methods.convert_coordinates_na_to_aa import path_to_indel_list, localize_vars, read_Walker, print_pos

drug_list = ['PZA', 'RIF', 'AMI', 'CAP', 'CIP', 'EMB', 'ETH', 'INH', 'KAN', 'MOX', 'OFX', 'PRO', 'STR']#

fdr_threshold = 0.1


def find_max_pval():
    with open('./max_pvals_fdr0.1', 'w') as f:
        for drug in drug_list:
            if exists('./' + drug + '/' + drug.lower() + '.lower.pvalue.pairs.fdr'):
                print(drug)
                max_pval = -1
                for line in open('./' + drug + '/' + drug.lower() + '.lower.pvalue.pairs.fdr').readlines()[1:]:
                    s = line.split('\t')
                    if float(s[2])/float(s[1]) < fdr_threshold:
                        pval = float(s[0])
                        if pval > max_pval:
                            max_pval = pval
                if max_pval >= 0:
                    f.write(drug + '\t' + str(max_pval) + '\n')


def print_filtered_p_values():
    for line in open('../../res/tree_was/max_pvals_fdr0.1').readlines():
        s = line.split('\t')
        drug = s[0]
        max_pval = float(s[1])
        print(drug)
        path_to_coords = '../../res/tree_was/' + drug + '/' + drug.lower() + '.site_pairs'
        path_to_p_values = '../../res/tree_was/' + drug + '/' + drug.lower() + '.lower.pvalue'
        out_path = '../../res/tree_was/filtered_lower_fdr01/' + drug.lower() + '.lower.pvalue.pairs'
        with open(out_path, 'w') as f:
            p_val_iter = open(path_to_p_values)
            coords = open(path_to_coords)
            f.write(coords.readline())
            for line in coords.readlines():
                p_val = p_val_iter.readline().strip()
                if float(p_val) <= max_pval:
                    f.write(line)


def convert_for_all_drugs():
    indel_list = [l.strip() for l in open(path_to_indel_list, 'r').readlines()]
    for drug in drug_list:
        path_to_coords = '../../res/tree_was/filtered_lower_fdr01/' + drug.lower() + '.lower.pvalue.pairs'
        if exists(path_to_coords):
            out_path = '../../res/tree_was/filtered_lower_fdr01_converted/' + drug.lower() + '.lower.pvalue.pairs.fdr.filtered.converted'
            snps = [int(line.strip().split('\t')[0]) for line in open(path_to_coords, 'r').readlines()[1:] if line != '\n']
            pos_to_cds, indel_pos_to_cds = localize_vars(snps, indel_list)
            gene_to_Walker = read_Walker()
            with open(out_path, 'w') as out_f:
                with open(path_to_coords, 'r') as f:
                    out_f.write('type\tgene_name\tgene_pos\t' + f.readline().strip() + '\tWalker\n')
                    for line in f.readlines():
                        if line != '\n':
                            pos = int(line.strip().split('\t')[0])
                            if pos > ref_len:
                                s = indel_list[pos - ref_len - 1].split('_')
                                del_pos = int(s[1])
                                cds = indel_pos_to_cds.get(del_pos)
                                print_pos(out_f, cds, del_pos, gene_to_Walker, line.strip(), s[0])
                            else:
                                cds = pos_to_cds.get(pos)
                                print_pos(out_f, cds, pos, gene_to_Walker, line.strip(), 'snp')


if __name__ == '__main__':
    # find_max_pval()
    print_filtered_p_values()
    convert_for_all_drugs()
