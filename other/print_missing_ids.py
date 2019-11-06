from os import listdir

from core.constants import data_path

path_to_res_folder = data_path + 'snps/gatk_before_cortex/gvcf_force_active/'
path_to_full_list = data_path + 'all_with_pheno.txt'
suffix = '_h37rv.g.vcf'
out_path = data_path + 'missing.list'


if __name__ == '__main__':
    sample_ids = [l.strip() for l in open(path_to_full_list).readlines()]
    processed = set()
    for fname in listdir(path_to_res_folder):
        i = fname.find(suffix)
        if i != -1:
            processed.add(fname[:i])
    with open(out_path, 'w') as f:
        for sample_id in sample_ids:
            if sample_id not in processed:
                f.write(sample_id + '\n')
