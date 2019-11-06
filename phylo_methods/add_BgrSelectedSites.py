from os.path import exists

drug_list = ['PZA', 'RIF', 'AMI', 'CAP', 'CIP', 'EMB', 'ETH', 'INH', 'KAN', 'MOX', 'OFL', 'PRO', 'STR']
path_to_filter = '../filtered_xparr_variant_pos_bam_filtered.list'
upper_or_lower = 'lower'

if __name__ == '__main__':
    for drug in drug_list:
        if exists('./' + drug):
            path = './' + drug + '/' + drug.lower() + '.' + upper_or_lower + '.pvalue.sites.fdr.prm'
            prm = None
            if exists(path):
                prm = open(path).readlines()
            else:
                path = './' + drug + '/' + drug.lower() + '.' + upper_or_lower + '.pvalue.pairs.fdr.prm'
                if exists(path):
                    prm = open(path).readlines()
            if prm is None:
                print('No prm found for ' + drug + ' ' + upper_or_lower)
                continue
            with open(path, 'w') as f:
                found = False
                for l in prm:
                    if l.startswith('BgrSelectedSites='):
                        found = True
                        f.write('BgrSelectedSites="' + path_to_filter + '"\n')
                    else:
                        f.write(l)
                        if not l.endswith('\n'):
                            f.write('\n')
                if not found:
                    f.write('BgrSelectedSites="' + path_to_filter + '"\n')