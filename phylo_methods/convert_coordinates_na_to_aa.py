from bisect import bisect_left

from src.core.annotations import read_annotations, localize_all_snps, CDSType, CDS
from src.core.constants import upstream_length, data_path, ref_len
from src.core.data_reading import read_h37rv

drug_list = ['PZA', 'RIF', 'AMI', 'CAP', 'CIP', 'EMB', 'ETH', 'INH', 'KAN', 'MOX', 'OFX', 'PRO', 'STR']#
path_to_dict = data_path + 'dictionaries_indels/Walker_dictionary.txt'#data_path + 'dictionaries_indels/Walker_dictionary.txt'
path_to_indel_list = data_path + 'xparr/mc10_mega_MP/indel_list.txt'#data_path + 'xparr/mc10_mega_MP/indel_list.txt'


def print_pos(out_f, cds, pos, gene_to_Walker, line, type, cds_list=None):
    if cds is None:
        if cds_list is None:
            out_f.write(type + '\tnon_cds\t' + str(ref_len - pos + 1) + '\t' + line + '\tmiss\n')
        else:
            left_gene, right_gene = find_closest_genes(cds_list, pos)
            if left_gene is None:
                out_f.write(type + '\tnon_cds\t' + str(ref_len - pos + 1) + '\t' + line + '\tmiss\t-\t' +
                            right_gene.name + '\n')
            elif right_gene is None:
                out_f.write(type + '\tnon_cds\t' + str(ref_len - pos + 1) + '\t' + line + '\tmiss\t' + left_gene.name +
                            '\t-\n')
            else:
                out_f.write(type + '\tnon_cds\t' + str(ref_len - pos + 1) + '\t' + line + '\tmiss\t' + left_gene.name +
                            '\t' + right_gene.name + '\n')
        return
    if cds.strand == 1:
        if cds.type == CDSType.Gene:
            gene_pos = str((pos - cds.start) // 3 + 1)
            if type in ('del', 'ins'):
                walker = gene_to_Walker.get(cds.name + '\t-')
            else:
                walker = gene_to_Walker.get(cds.name + '\t' + gene_pos)
            if walker is not None:
                if cds_list is None:
                    out_f.write(type + '\t' + cds.name + '\t' + gene_pos + '\t' + line + '\t' + walker + '\n')
                else:
                    out_f.write(type + '\t' + cds.name + '\t' + gene_pos + '\t' + line + '\t' + walker + '\t\t\n')
            else:
                if cds_list is None:
                    out_f.write(type + '\t' + cds.name + '\t' + gene_pos + '\t' + line + '\tmiss\n')
                else:
                    out_f.write(type + '\t' + cds.name + '\t' + gene_pos + '\t' + line + '\tmiss\t\t\n')
        elif cds.type == CDSType.upstream:
            if type in ('del', 'ins'):
                walker = gene_to_Walker.get(cds.name + '\t-')
            else:
                walker = gene_to_Walker.get(cds.name + '\t' + str(pos - cds.end - 1))
            if walker is not None:
                if cds_list is None:
                    out_f.write(type + '\t' + cds.name + '\t' + str(pos - cds.end - 1) + '\t' + line + '\t' + walker + '\n')
                else:
                    out_f.write(type + '\t' + cds.name + '\t' + str(pos - cds.end - 1) + '\t' + line + '\t' + walker + '\t\t\n')
            else:
                if cds_list is None:
                    out_f.write(type + '\t' + cds.name + '\t' + str(pos - cds.end - 1) + '\t' + line + '\tmiss\n')
                else:
                    out_f.write(type + '\t' + cds.name + '\t' + str(pos - cds.end - 1) + '\t' + line + '\tmiss\t\t\n')
        else:
            if type in ('del', 'ins'):
                walker = gene_to_Walker.get(cds.name + '\t-')
            else:
                walker = gene_to_Walker.get(cds.name + '\t' + str(pos - cds.start + 1))
            if walker is not None:
                if cds_list is None:
                    out_f.write(type + '\t' + cds.name + '\t' + str(pos - cds.start + 1) + '\t' + line + '\t' + walker + '\n')
                else:
                    out_f.write(type + '\t' + cds.name + '\t' + str(pos - cds.start + 1) + '\t' + line + '\t' + walker + '\t\t\n')
            else:
                if cds_list is None:
                    out_f.write(type + '\t' + cds.name + '\t' + str(pos - cds.start + 1) + '\t' + line + '\tmiss\n')
                else:
                    out_f.write(type + '\t' + cds.name + '\t' + str(pos - cds.start + 1) + '\t' + line + '\tmiss\t\t\n')
    else:
        if cds.type == CDSType.Gene:
            gene_pos = str((cds.end - pos) // 3 + 1)
            if type in ('del', 'ins'):
                walker = gene_to_Walker.get(cds.name + '\t-')
            else:
                walker = gene_to_Walker.get(cds.name + '\t' + gene_pos)
            if walker is not None:
                if cds_list is None:
                    out_f.write(type + '\t' + cds.name + '\t' + gene_pos + '\t' + line + '\t' + walker + '\n')
                else:
                    out_f.write(type + '\t' + cds.name + '\t' + gene_pos + '\t' + line + '\t' + walker + '\t\t\n')
            else:
                if cds_list is None:
                    out_f.write(type + '\t' + cds.name + '\t' + gene_pos + '\t' + line + '\tmiss\n')
                else:
                    out_f.write(type + '\t' + cds.name + '\t' + gene_pos + '\t' + line + '\tmiss\t\t\n')
        elif cds.type == CDSType.upstream:
            if type in ('del', 'ins'):
                walker = gene_to_Walker.get(cds.name + '\t-')
            else:
                walker = gene_to_Walker.get(cds.name + '\t' + str(cds.start - pos - 1))
            if walker is not None:
                if cds_list is None:
                    out_f.write(type + '\t' + cds.name + '\t' + str(cds.start - pos - 1) + '\t' + line + '\t' + walker + '\n')
                else:
                    out_f.write(type + '\t' + cds.name + '\t' + str(cds.start - pos - 1) + '\t' + line + '\t' + walker + '\t\t\n')
            else:
                if cds_list is None:
                    out_f.write(type + '\t' + cds.name + '\t' + str(cds.start - pos - 1) + '\t' + line + '\tmiss\n')
                else:
                    out_f.write(type + '\t' + cds.name + '\t' + str(cds.start - pos - 1) + '\t' + line + '\tmiss\t\t\n')
        else:
            if type in ('del', 'ins'):
                walker = gene_to_Walker.get(cds.name + '\t-')
            else:
                walker = gene_to_Walker.get(cds.name + '\t' + str(cds.end - (pos - 1)))
            if walker is not None:
                if cds_list is None:
                    out_f.write(type + '\t' + cds.name + '\t' + str(cds.end - (pos - 1)) + '\t' + line + '\t' + walker + '\n')
                else:
                    out_f.write(type + '\t' + cds.name + '\t' + str(cds.end - (pos - 1)) + '\t' + line + '\t' + walker + '\t\t\n')
            else:
                if cds_list is None:
                    out_f.write(type + '\t' + cds.name + '\t' + str(cds.end - (pos - 1)) + '\t' + line + '\tmiss\n')
                else:
                    out_f.write(type + '\t' + cds.name + '\t' + str(cds.end - (pos - 1)) + '\t' + line + '\tmiss\t\t\n')


def localize_vars(snps: 'list[int]', indel_list: 'list[str]'=None)->'tuple[dict[int, list[CDS]], dict[int, list[CDS]]]':
    """
    Reads annotations and maps snps and indels to CDS by positions

    :param snps: list[int] list of snp pos
    :param indel_list: list[str] list of indels
    :return: tuple[dict[int, list[CDS]], dict[int, list[CDS]]] map from snp pos to CDS and map from indel pos to CDS
    """
    cds_list = read_annotations(upstream_length)
    pos_to_cds = localize_all_snps(snps, cds_list)
    if indel_list is None:
        return pos_to_cds, None
    indel_pos_set = set()
    for m in indel_list:
        if m != '':
            indel_pos_set.add(int(m.split('_')[1]))
    indel_pos_list = list(indel_pos_set)
    indel_pos_to_cds = localize_all_snps(indel_pos_list, cds_list)
    return pos_to_cds, indel_pos_to_cds


def read_Walker():
    gene_to_Walker = {}
    with open(path_to_dict, 'r') as f:
        for line in f.readlines():
            s = line.strip().split('\t')
            gene_to_Walker[s[1] + '\t' + s[2]] = s[-1]
    return gene_to_Walker


def convert_for_all_drugs():
    indel_list = [l.strip() for l in open(path_to_indel_list, 'r').readlines()]
    for drug in drug_list:
        path_to_coords = '../../res/tree_was/' + drug + '/' + drug.lower() + '.site_pairs'#'../../res/tree_was/filtered/' + drug.lower() + '.upper.pvalue.pairs'
        out_path = '../../res/tree_was/site_pairs_converted/' + drug.lower() + '.site_pairs.converted'#'../../res/tree_was/converted/' + drug.lower() + '.upper.pvalue.pairs.converted'
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


def print_pos_simple(out_f, cds, pos, gene_to_Walker, type):
    if cds is None:
        out_f.write(type + '\tnon_cds\t' + str(ref_len - pos + 1) + '\tmiss\n')
        return
    if cds.strand == 1:
        if cds.type == CDSType.Gene:
            gene_pos = str((pos - cds.start) // 3 + 1)
            if type in ('del', 'ins'):
                walker = gene_to_Walker.get(cds.name + '\t-')
            else:
                walker = gene_to_Walker.get(cds.name + '\t' + gene_pos)
            if walker is not None:
                out_f.write(type + '\t' + cds.name + '\t' + gene_pos + '\t' + walker)
            else:
                out_f.write(type + '\t' + cds.name + '\t' + gene_pos + '\tmiss\n')
        elif cds.type == CDSType.upstream:
            if type in ('del', 'ins'):
                walker = gene_to_Walker.get(cds.name + '\t-')
            else:
                walker = gene_to_Walker.get(cds.name + '\t' + str(pos - cds.end - 1))
            if walker is not None:
                out_f.write(type + '\t' + cds.name + '\t' + str(pos - cds.end - 1) + '\t' + walker)
            else:
                out_f.write(type + '\t' + cds.name + '\t' + str(pos - cds.end - 1) + '\tmiss\n')
        else:
            if type in ('del', 'ins'):
                walker = gene_to_Walker.get(cds.name + '\t-')
            else:
                walker = gene_to_Walker.get(cds.name + '\t' + str(pos - cds.start + 1))
            if walker is not None:
                out_f.write(type + '\t' + cds.name + '\t' + str(pos - cds.start + 1) + '\t' + walker)
            else:
                out_f.write(type + '\t' + cds.name + '\t' + str(pos - cds.start + 1) + '\tmiss\n')
    else:
        if cds.type == CDSType.Gene:
            gene_pos = str((cds.end - pos - 1) // 3 + 1)
            if type in ('del', 'ins'):
                walker = gene_to_Walker.get(cds.name + '\t-')
            else:
                walker = gene_to_Walker.get(cds.name + '\t' + gene_pos)
            if walker is not None:
                out_f.write(type + '\t' + cds.name + '\t' + gene_pos + '\t' + walker)
            else:
                out_f.write(type + '\t' + cds.name + '\t' + gene_pos + '\tmiss\n')
        elif cds.type == CDSType.upstream:
            if type in ('del', 'ins'):
                walker = gene_to_Walker.get(cds.name + '\t-')
            else:
                walker = gene_to_Walker.get(cds.name + '\t' + str(cds.start - pos - 1))
            if walker is not None:
                out_f.write(type + '\t' + cds.name + '\t' + str(cds.start - pos - 1) + '\t' + walker)
            else:
                out_f.write(type + '\t' + cds.name + '\t' + str(cds.start - pos - 1) + '\tmiss\n')
        else:
            if type in ('del', 'ins'):
                walker = gene_to_Walker.get(cds.name + '\t-')
            else:
                walker = gene_to_Walker.get(cds.name + '\t' + str(cds.end - (pos - 1)))
            if walker is not None:
                out_f.write(type + '\t' + cds.name + '\t' + str(cds.end - (pos - 1)) + '\t' + walker)
            else:
                out_f.write(type + '\t' + cds.name + '\t' + str(cds.end - (pos - 1)) + '\tmiss\n')


def convert_full_indel_list():
    indel_list = [l.strip() for l in open(path_to_indel_list, 'r').readlines()]
    out_path = path_to_indel_list + '.converted'
    pos_to_cds, indel_pos_to_cds = localize_vars([], indel_list)
    gene_to_Walker = read_Walker()
    with open(out_path, 'w') as out_f:
        out_f.write('type\tgene_name\tgene_pos\tWalker\n')
        for line in indel_list:
            s = line.strip().split('_')
            pos = int(s[1])
            cds = indel_pos_to_cds.get(pos)
            print_pos_simple(out_f, cds, pos, gene_to_Walker, s[0])


def convert_vars_to_drugs():
    indel_list = [l.strip() for l in open(path_to_indel_list, 'r').readlines()]
    path_to_coords = '../../res/tree_was/tr(vars2drugs).cor2pcor.R.out'
    out_path = '../../res/tree_was/tr(vars2drugs).cor2pcor.R.out.converted'
    snps = [int(line.strip().split('\t')[1]) for line in open(path_to_coords, 'r').readlines()[1:] if line != '\n']
    pos_to_cds, indel_pos_to_cds = localize_vars(snps, indel_list)
    gene_to_Walker = read_Walker()
    with open(out_path, 'w') as out_f:
        with open(path_to_coords, 'r') as f:
            out_f.write('type\tgene_name\tgene_pos\t' + f.readline().strip() + '\tWalker\n')
            for line in f.readlines():
                if line != '\n':
                    pos = int(line.strip().split('\t')[1])
                    if pos > ref_len:
                        s = indel_list[pos - ref_len - 1].split('_')
                        del_pos = int(s[1])
                        cds = indel_pos_to_cds.get(del_pos)
                        print_pos(out_f, cds, del_pos, gene_to_Walker, line.strip(), s[0])
                    else:
                        cds = pos_to_cds.get(pos)
                        print_pos(out_f, cds, pos, gene_to_Walker, line.strip(), 'snp')


def convert_single_file(header=True):
    indel_list = [l.strip() for l in open(path_to_indel_list, 'r').readlines()]
    path_to_coords = '../../res/tree_was/vars.set.txt'
    out_path = '../../res/tree_was/vars.set.converted'
    if header:
        snps = [int(line.strip().split('\t')[0]) for line in open(path_to_coords, 'r').readlines()[1:] if line != '\n']
    else:
        snps = [int(line.strip().split('\t')[0]) for line in open(path_to_coords, 'r').readlines() if line != '\n']
    pos_to_cds, indel_pos_to_cds = localize_vars(snps, indel_list)
    gene_to_Walker = read_Walker()
    with open(out_path, 'w') as out_f:
        with open(path_to_coords, 'r') as f:
            if header:
                out_f.write('type\tgene_name\tgene_pos\t' + f.readline().strip() + '\tWalker\n')
            else:
                out_f.write('type\tgene_name\tgene_pos\txparr_pos\tWalker\n')
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


def find_closest_genes(cds_list, pos):
    cds_starts = [cds.start for cds in cds_list]
    i = bisect_left(cds_starts, pos)
    if i == 0:
        return (None, cds_list[0])
    elif i == len(cds_list):
        return (cds_list[i - 1], None)
    else:
        return (cds_list[i - 1], cds_list[i])


def convert_single_file_with_closest_genes(header=True):
    cds_list = read_annotations(upstream_length)
    indel_list = [l.strip() for l in open(path_to_indel_list, 'r').readlines()]
    path_to_coords = '../../res/tree_was/vars.set.txt'
    out_path = '../../res/tree_was/vars.set.converted_with_genes'
    if header:
        snps = [int(line.strip().split('\t')[0]) for line in open(path_to_coords, 'r').readlines()[1:] if line != '\n']
    else:
        snps = [int(line.strip().split('\t')[0]) for line in open(path_to_coords, 'r').readlines() if line != '\n']
    pos_to_cds, indel_pos_to_cds = localize_vars(snps, indel_list)
    gene_to_Walker = read_Walker()
    with open(out_path, 'w') as out_f:
        with open(path_to_coords, 'r') as f:
            if header:
                out_f.write('type\tgene_name\tgene_pos\t' + f.readline().strip() + '\tWalker\tleft_gene\tright_gene\n')
            else:
                out_f.write('type\tgene_name\tgene_pos\txparr_pos\tWalker\tleft_gene\tright_gene\n')
            for line in f.readlines():
                if line != '\n':
                    pos = int(line.strip().split('\t')[0])
                    if pos > ref_len:
                        s = indel_list[pos - ref_len - 1].split('_')
                        del_pos = int(s[1])
                        cds = indel_pos_to_cds.get(del_pos)
                        print_pos(out_f, cds, del_pos, gene_to_Walker, line.strip(), s[0], cds_list)
                    else:
                        cds = pos_to_cds.get(pos)
                        print_pos(out_f, cds, pos, gene_to_Walker, line.strip(), 'snp', cds_list)


if __name__ == '__main__':
    convert_single_file(False)
    # convert_for_all_drugs()