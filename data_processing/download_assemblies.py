import time
from urllib.error import HTTPError

from Bio import Entrez

from core.constants import data_path

path_to_list = data_path + 'ref_strains.list1'
path_for_saving = data_path + 'ref_strains/'


def fetch_given_seq(seq_id):
    Entrez.email = "gennady.fedonin@gmail.com"
    res_file = open(path_for_saving + seq_id + '.fasta', 'w')
    try:
        handle = Entrez.efetch(db="nucleotide", rettype='fasta_cds_na', retmode="txt", id=seq_id)
    except HTTPError:
        print('### get an exception ###')
        time.sleep(10)
        try:
            handle = Entrez.efetch(db="nucleotide", rettype='fasta_cds_na', retmode="txt", id=seq_id)
        except HTTPError:
            print('### fetch failed ###')
            failed_refseq = open(path_for_saving + 'failed_refseq.txt', 'a')
            failed_refseq.write(seq_id + '\n')
            failed_refseq.close()
            return 'fail'

    res_file.write(handle.read())
    handle.close()
    res_file.close()

    return 'success'


if __name__ == '__main__':
    # for l in open(path_to_list).readlines():
    #     fetch_given_seq(l.strip())
    fetch_given_seq('GCA_000270345.1')