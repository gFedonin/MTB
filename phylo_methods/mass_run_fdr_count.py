import os

from os.path import exists
from sklearn.externals.joblib import Parallel, delayed

drug_list = ['PZA', 'RIF', 'AMI', 'CAP', 'CIP', 'EMB', 'ETH', 'INH', 'KAN', 'MOX', 'OFL', 'PRO', 'STR']
upper_or_lower = 'upper'

epistat = '~/devel/epistat.6/estimate_fdr.pl '
in_suffix1 = '.' + upper_or_lower + '.pvalue.pairs.fdr.prm'
out_suffix1 = '.' + upper_or_lower + '.pvalue.pairs_filter.fdr'
in_suffix2 = '.' + upper_or_lower + '.pvalue.sites.fdr.prm'
out_suffix2 = '.' + upper_or_lower + '.pvalue.sites_filter.fdr'


def run(drug):
    if exists('./' + drug):
        os.chdir('./' + drug)
        in_path = drug.lower() + in_suffix1
        if exists(in_path):
            os.system(epistat + in_path + ' > ' + drug.lower() + out_suffix1 + ' 2> ' + drug.lower() + '.' +
                      upper_or_lower + '.log')
        else:
            in_path = drug.lower() + in_suffix2
            if exists(in_path):
                os.system(epistat + in_path + ' > ' + drug.lower() + out_suffix2 + ' 2> ' + drug.lower() + '.' +
                          upper_or_lower + '.log')
        return 1
    return 0


if __name__ == '__main__':
    tasks = Parallel(n_jobs=len(drug_list))(delayed(run)(drug) for drug in drug_list)
    c = 0
    for task in tasks:
        c += task
    print('computed for %d drugs' % c)
