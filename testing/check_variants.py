from sklearn.externals.joblib import Parallel, delayed

path_to_snps = './data/snps/raw_with_DR_with_indel_with_pheno_and_snp_mc10/'
path_to_ids = './data/all_with_pheno_and_snp.txt'

def read_variants(sample_id):
    snps = []
    with open(path_to_snps + sample_id + '.variant', 'r') as f1:
        lines = f1.readlines()
        for line in lines:
            s = line.strip().split('\t')
            snps.append((int(s[0]), s[1], s[2]))
    return sample_id, snps


sample_ids = [sample_id[:-1] for sample_id in open(path_to_ids, 'r').readlines()]
tasks = Parallel(n_jobs=-1)(
    delayed(read_variants)(sample_id) for sample_id in sample_ids)
for sample_id, variants in tasks:
    for pos, alt, v_type in variants:
        if v_type == 'snp' and len(alt) != 1:
            print(sample_id + ' ' + str(pos))
        elif v_type == 'ins' and ('-' in alt or len(alt) == 0):
            print(sample_id + ' ' + str(pos))