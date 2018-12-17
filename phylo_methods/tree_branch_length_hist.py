import matplotlib.pyplot as plt

path_to_branch_lengths = '../../data/tree_with_pheno_and_snp_mc5_mega.branch_lengths'
out_path = '../../data/tree_with_pheno_and_snp_mc5_hist_mega.png'

if __name__ == '__main__':
    branch_lengths = [float(l.strip()) for l in open(path_to_branch_lengths, 'r').readlines()]
    plt.title('Histogram of branch lengths')
    plt.xlabel('Branch length')
    plt.ylabel('Number of nodes (incl. tips)')
    n, bins, patches = plt.hist(branch_lengths, 50, density=True, facecolor='g', alpha=0.75)
    plt.axis([0, 0.002, 0, 6000])
    plt.savefig(out_path)
