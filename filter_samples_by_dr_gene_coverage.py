from bisect import bisect_left, bisect_right

from sklearn.externals.joblib import Parallel, delayed

from annotations import read_annotations
from constants import dr_genes

path_to_ids = './data/all_mtb_with_pheno_and_snp.txt'
out_path_gene_cov = './res/sample_gen_cov.txt'
path_to_depths = './data/coverages_with_percentiles/t5p5/'
path_to_annotations = './data/AL123456_rev.gff'
out_path = './data/dr_covered_mtb_with_pheno_and_snp.txt'

gene_cov_threshold = 0.7
upstream_length = 100


def get_gene_annotations():
    cds_list = read_annotations(path_to_annotations, upstream_length)
    res = {}
    gene_names = set(dr_genes)
    for cds in cds_list:
        if cds.name in gene_names or cds.synonym in gene_names:
            if cds.name not in res.keys():
                res[cds.name] = cds
            else:
                cds_prev = res[cds.name]
                cds_prev.start = min(cds_prev.start, cds.start)
                cds_prev.end = max(cds_prev.end, cds.end)
    return list(res.values())


def test_sample_cov(sample_id, genes):
    coverage = []
    cov_starts = []
    with open(path_to_depths + sample_id + '.coverage') as f:
        f.readline()
        f.readline()
        f.readline()
        for line in f.readlines():
            s = line.split()
            st = int(s[0])
            coverage.append((st, int(s[1])))
            cov_starts.append(st)
    # coverage.sort(key=lambda tup: tup[0])
    gene_to_cov = []
    covered_genes_num = 0
    for cds in genes:
        i = bisect_left(cov_starts, cds.start)
        j = bisect_right(cov_starts, cds.end, lo=i)
        gen_len = cds.end - cds.start + 1
        cov_len = 0
        if i > 0 and cds.start < coverage[i - 1][1]:
            cov_len += min(coverage[i - 1][1], cds.end) - cds.start + 1
        for k in range(i, j - 1):
            cov_len += coverage[k][1] - coverage[k][0] + 1
        if j - 1 < len(coverage) and i != j:
            cov_len += min(coverage[j - 1][1], cds.end) - coverage[j - 1][0] + 1
        cov = cov_len / gen_len
        gene_to_cov.append((cds.name, str(cov)))
        if cov >= gene_cov_threshold:
            covered_genes_num += 1
    return sample_id, covered_genes_num, gene_to_cov


def main():
    genes = get_gene_annotations()
    with open(path_to_ids, 'r') as f:
        with open(out_path_gene_cov, 'w') as out_cov:
            with open(out_path, 'w') as out:
                tasks = Parallel(n_jobs=-1)(delayed(test_sample_cov)(name.strip(), genes) for name in f.readlines())
                for task in tasks:
                    sample_id, covered_genes_num, gene_to_cov = task
                    out_cov.write(sample_id + '\n')
                    for name, cov in gene_to_cov:
                        out_cov.write(name + ' ' + cov + '\n')
                    if covered_genes_num == len(genes):
                        out.write(sample_id + '\n')


if __name__ == '__main__':
    main()
