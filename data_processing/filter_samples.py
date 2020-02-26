from core.constants import data_path

list_to_filter = data_path + 'combined_no_dup.list'
no_pheno = data_path  + 'no_pheno.txt'
list_with_mut = data_path + 'combined_filtered.list'

out_path = data_path + 'combined_no_dup_with_pheno.list'


if __name__ == '__main__':
    drop_list = set(l for l in open(no_pheno).readlines())
    keep_list = set(l for l in open(list_with_mut).readlines())
    with open(out_path, 'w') as f:
        for l in open(list_to_filter).readlines():
            if l in drop_list:
                continue
            if l in keep_list:
                f.write(l)