from Bio.Seq import Seq
from sklearn.metrics import confusion_matrix, make_scorer, f1_score, roc_auc_score

# data_path = '/export/data/fedonin/MTB/data/'
data_path = '../../data/'

complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', '-': '-'}


aminoacids = ['C', 'D', 'S', 'Q', 'K', 'I', 'P', 'T', 'F', 'N', 'G', 'H', 'L', 'R', 'W', 'A', 'V', 'E', 'Y', 'M']


codon_table = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
            'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
            'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
            'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
            'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
            'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
            'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
            'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
            'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
            'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
            'GGG': 'G', 'TAA': '*', 'TAG': '*', 'TGA': '*'}

codon_table_compl = {}
for k,v in codon_table.items():
    codon_table_compl[str(Seq(k).reverse_complement())] = v

stop_codons = ['TAA', 'TAG', 'TGA']

start_codons = ['TTG', 'CTG', 'ATG']

drug_names = {'STR': 'Streptomycin', 'INH': 'Isoniazid', 'EMB': 'Ethambutol', 'ETH': 'Ethionamide',
              'CIP': 'Ciprofloxacin', 'OFX': 'Ofloxacin', 'PZA': 'Pyrazinamide', 'RIF': 'Rifampicin',
              'AMK': 'Amikacin', 'CAP': 'Capreomycin', 'KAN': 'Kanamycin', 'SM': 'Streptomycin',
              'MOX': 'Moxifloxacin', 'PTH': 'Prothionamide', 'RPM': 'Rifampicin','AK': 'Amikacin',
              'FLQ':['Ciprofloxacin', 'Moxifloxacin', 'Ofloxacin'], 'AMI':'Amikacin'}
              # 'AMI': ['Amikacin', 'Capreomycin', 'Kanamycin']}


dr_genes = ['ahpC', 'eis', 'embA', 'embB', 'embC', 'embR', 'fabG1', 'gid', 'gyrA', 'gyrB', 'inhA', 'iniA', 'iniC',
             'katG', 'manB', 'ndh', 'pncA', 'rmlD', 'rpoB', 'rpsA', 'rpsL', 'rrs', 'tlyA', 'ethR', 'fpbC', 'iniB',
              'kasA', 'ethA', 'fabD', 'efpA', 'thyA', 'panD', 'accD6', 'fbpC', 'nat', 'folC', 'rrl', 'rpoC', 'ribD', 'rplC']

upstream_length = 100
ref_len = 4411532

def tp(y_true, y_pred): return confusion_matrix(y_true, y_pred)[1, 1]
def tn(y_true, y_pred): return confusion_matrix(y_true, y_pred)[0, 0]
def fp(y_true, y_pred): return confusion_matrix(y_true, y_pred)[0, 1]
def fn(y_true, y_pred): return confusion_matrix(y_true, y_pred)[1, 0]

custom_scoring = {'tp' : make_scorer(tp), 'tn' : make_scorer(tn), 'fp' : make_scorer(fp), 'fn' : make_scorer(fn), 'f1' : f1_score, 'auc': roc_auc_score}
