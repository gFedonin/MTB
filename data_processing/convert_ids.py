from os import makedirs, listdir
from os.path import exists
from subprocess import check_call, call

from core.constants import data_path

data_set = 'patric'
# data_set = 'missing'
# data_set = 'combined_nocoll'
path_to_list_in = data_path + data_set + '/' + data_set + '.list'
# path_to_list_in = data_path + data_set + '/Genoscreen.list'
# path_to_list_in = data_path + data_set + '/missing.list'
path_to_list_out = data_path + data_set + '/' + data_set + '_converted.samples'
metadata_path = data_path + data_set + '/metadata/'
# metadata_path = data_path + 'snps/combined_nocoll/combined_metadata/'
# metadata_path = data_path + data_set + '/Genoscreen_metadata/'


batch_size = 100


def get_metadata():
    if not exists(metadata_path):
        makedirs(metadata_path)
    samples = [l.strip() for l in open(path_to_list_in).readlines()]
    for i in range(0, len(samples), batch_size):
        if not exists(metadata_path + str(i) + '.csv'):
            terms = '+OR+'.join(samples[i: min(i + batch_size, len(samples))])
            check_call('wget -O ' + metadata_path + str(i) + '.csv ' +
                       '\'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=' +
                       terms + '\'', shell=True)


def get_biosample_from_sra():
    samples = [l.strip() for l in open(path_to_list_in).readlines()]
    with open(path_to_list_out, 'w') as f:
        for i in range(0, len(samples), batch_size):
            for l in open(metadata_path + str(i) + '.csv').readlines()[1:-1]:
                f.write(l.strip().split(',')[25] + '\n')


path_to_cryptic = data_path + data_set + '/CRyPTIC.csv'
# path_to_cryptic_mapping = data_path + data_set + '/CRyPTIC_map.csv'
path_to_cryptic_mapping = data_path + data_set + '/Genoscreen.csv'
path_to_cryptic_not_found = data_path + data_set + '/Genoscreen.missing'


def get_err_and_biosample():
    samples = [l.strip() for l in open(path_to_list_in).readlines()]
    lib_name_to_ids = []
    for l in open(path_to_cryptic).readlines()[1:]:
        s = l.split(',')
        if len(s) < 26:
            print(s)
            exit(1)
        lib_name_to_ids.append((s[11], s[0], s[25]))
    with open(path_to_cryptic_mapping, 'w') as f:
        with open(path_to_cryptic_not_found, 'w') as fm:
            for sample in samples:
                l_name = sample + '.'
                found = False
                for lib_name, err, samea in lib_name_to_ids:
                    if l_name in lib_name:
                        found = True
                        f.write(sample + '\t' + err + '\t' + samea + '\n')
                        break
                if not found:
                    fm.write(sample + '\n')


# old_samples = data_path + 'all_with_pheno.txt'
old_samples = data_path + 'base_coll18_walker18_johnsen19_farhat19.samples'
# new_samples = data_path + data_set + '/' + data_set + '_converted.samples'
new_samples = data_path + data_set + '/' + 'cryptic_converted.samples'
# really_new_samples = data_path + data_set + '/' + data_set + '_converted_new.samples'
really_new_samples = data_path + data_set + '/' + 'cryptic_converted_new.samples'


def pick_new():
    old_set = {l.strip() for l in open(old_samples).readlines()}
    with open(really_new_samples, 'w') as f:
        for l in open(new_samples).readlines():
            if l.strip() not in old_set:
                f.write(l)


really_new_samples_err = data_path + data_set + '/' + data_set + '_new.samples'


def id_mapping():
    samples = [l.strip() for l in open(path_to_list_in).readlines()]
    samn_to_err = {}
    for i in range(0, len(samples), batch_size):
        for l in open(metadata_path + str(i) + '.csv').readlines()[1:-1]:
            s = l.strip().split(',')
            samn_to_err[s[25]] = s[0]
    with open(really_new_samples_err, 'w') as f:
        for l in open(really_new_samples).readlines():
            f.write(samn_to_err[l.strip()] + '\n')


path_to_downloaded = data_path + data_set + '/downloaded.samples'
path_to_canonical = data_path + data_set + '/downloaded_canonical.samples'


def id_to_run_id():
    metadata = []
    samples = [l.strip() for l in open(path_to_list_in).readlines()]
    for i in range(0, len(samples), batch_size):
        for l in open(metadata_path + str(i) + '.csv').readlines()[1:-1]:
            metadata.append(l.strip())
    with open(path_to_canonical, 'w') as f:
        for l in open(path_to_downloaded).readlines():
            err = l.strip()
            for meta in metadata:
                if err in meta:
                    f.write(meta.split(',')[0] + '\n')
                    break


# path_to_err = data_path + data_set + '/' + 'covered.list'
path_to_err = data_path + data_set + '/' + data_set + '_trimmed.list'
path_to_samea = data_path + data_set + '/' + data_set + '_samea.samples'
path_to_mapping = data_path + data_set + '/' + data_set + '_id_mapping.csv'


def get_samea():
    can_to_samea = {}
    for fname in listdir(metadata_path):
        for l in open(metadata_path + fname).readlines()[1:-1]:
            s = l.split(',')
            if len(s) >= 1:
                can_to_samea[s[0]] = s[25]
    with open(path_to_samea, 'w') as f1:
        with open(path_to_mapping, 'w') as f2:
            for l in open(path_to_err).readlines():
                sid = l.strip()
                if sid in can_to_samea:
                    f1.write(can_to_samea[sid] + '\n')
            for l in open(path_to_err).readlines():
                sid = l.strip()
                if sid in can_to_samea:
                    f2.write(l.strip() + '\t' + can_to_samea[sid] + '\n')


path_to_files = data_path + 'snps/combined_raw_variants_mq40_keep_complex/'
out_path = data_path + 'snps/combined_raw_variants_mq40_keep_complex_std_names/'
# path_to_merged_metadata = data_path + 'snps/combined_nocoll/combined_metadata/merged'
path_to_merged_metadata = data_path + 'snps/combined_nocoll/combined_metadata/merged.metadata'
exclude_list = data_path + 'coll18/coll18.list'
no_meta_list = data_path + 'snps/combined_nocoll/combined_metadata/nometa.list'
repeats_list = data_path + 'snps/combined_nocoll/combined_metadata/repeats.list'


def rename_files():
    if not exists(out_path):
        makedirs(out_path)
    meta_to_samea = []
    for l in open(path_to_merged_metadata).readlines()[1:-1]:
        s = l.split(',')
        if len(s) > 1:
            meta_to_samea.append((l, s[25]))
    exclude_set = {l.strip() for l in open(exclude_list).readlines()}
    c = 0
    c1 = 0
    c2 = 0
    c3 = 0
    c4 = 0
    repeats_f = open(repeats_list, 'w')
    nometa_f = open(no_meta_list, 'w')
    for fname in listdir(path_to_files):
        if fname.endswith('.variants'):
            c += 1
            sample = fname[0: fname.find('.')]
            if 'SAMN' in sample or 'SAMEA' in sample:
                if not exists(out_path + fname):
                    call('cp ' + path_to_files + fname + ' ' + out_path, shell=True)
                else:
                    c3 += 1
                c1 += 1
            else:
                if not sample in exclude_set:
                    samea = None
                    for meta, samea_id in meta_to_samea:
                        if sample in meta:
                            samea = samea_id
                            break
                    if not samea is None:
                        if not exists(out_path + samea + '.variants'):
                            call('cp ' + path_to_files + fname + ' ' + out_path, shell=True)
                            call('mv ' + out_path + fname + ' ' + out_path + samea + '.variants', shell=True)
                        else:
                            repeats_f.write(sample + '\t' + samea + '\n')
                            c3 += 1
                        c2 += 1
                        # call('mv ' + out_path + fname + ' ' + out_path + samea + '_h37rv.g.vcf', shell=True)
                        # call('cp ' + path_to_files + sample + '_h37rv.g.vcf.idx ' + out_path, shell=True)
                        # call('mv ' + out_path + sample + '_h37rv.g.vcf.idx ' + out_path + samea + '_h37rv.g.vcf.idx', shell=True)
                    else:
                        nometa_f.write(sample + '\n')
                        print(sample)
                else:
                    c4 += 1
                    if not exists(out_path + sample + '.variants'):
                        call('cp ' + path_to_files + fname + ' ' + out_path, shell=True)
    repeats_f.close()
    nometa_f.close()
    print('total ' + str(c))
    print('samea ' + str(c1))
    print('converted ' + str(c2))
    print('existed ' + str(c3))
    print('excluded ' + str(c4))


def find_non_tuberculosis():
    for l in open(path_to_merged_metadata).readlines()[1:-1]:
        if 'Run' not in l and l.strip() != '':
            if ',Mycobacterium tuberculosis,' not in l:
                print(l.strip())


nontub_path = data_path + 'nontub.list'
filtered_list_path = data_path + 'big_filtered.list'


def get_full_list():
    nontub_err = set()
    nontub_samea = set()
    for l in open(nontub_path).readlines():
        s = l.strip().split('\t')
        nontub_err.add(s[0])
        if len(s) > 1:
            nontub_samea.add(s[1])
    with open(filtered_list_path, 'w') as f:
        for fname in listdir(out_path):
            sample = fname[:fname.index('.')]
            if sample not in nontub_err and sample not in nontub_samea:
                f.write(sample + '\n')


if __name__ == '__main__':
    # get_err_and_biosample()
    get_metadata()
    # get_biosample_from_sra()
    # pick_new()
    # id_mapping()
    # id_to_canonical()
    # get_samea()
    # rename_files()
    # find_non_tuberculosis()
    # get_full_list()