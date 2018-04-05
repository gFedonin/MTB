path_to_full_set = './data/dr_covered_with_pheno_and_snp.txt'#'./data/all_with_pheno_ids.txt'
filter_paths = ['./data/subsets/Walker_subset.txt']#['./data/No_pheno.txt', './data/No_snp.txt', './data/trash.txt']['./data/No_snp.txt', './data/trash.txt']
out_path = './data/dr_covered_with_pheno_and_snp_Walker.txt'#''./data/dr_covered_with_pheno_and_snp_no_Walker.txt'./data/all_mtb_with_pheno_and_snp.txt'

intersect = True

full_set = []
with open(path_to_full_set, 'r') as f:
    for line in f.readlines():
        full_set.append(line.strip())

filter_set = set()
for path in filter_paths:
    with open(path, 'r') as f:
        for line in f.readlines():
            filter_set.add(line.strip())
print(str(len(filter_set)) + ' in filter set')

c = 0
with open(out_path, 'w') as f:
    for line in full_set:
        if line in filter_set:
            if intersect:
                f.write(line + '\n')
            c += 1
        else:
            if not intersect:
                f.write(line + '\n')

if intersect:
    print(str(c) + ' selected')
else:
    print(str(c) + ' filtered out')
