from matplotlib import pyplot as plt
import math
import numpy as np

path_to_mut_counts = '../res/testing/gatk_before_vs_freebayes_before_std3.freebayes_unique'
out_path = '../res/testing/freebayes_before.png'
title = 'freebayes_unique'


if __name__ == "__main__":
    data = []
    for l in open(path_to_mut_counts).readlines():
        s = l.strip().split('\t')
        data.append(int(s[1]))
    weights = np.ones_like(data)/float(len(data))
    bins = np.linspace(math.ceil(min(data)), math.floor(max(data)), 30)    
    # plt.xlim([min(data)-5, max(data)+5])
    plt.hist(data, bins=bins, alpha=0.5, weights=weights)
    plt.title(title)
    plt.xlabel('sample_num')
    plt.ylabel('count')

    plt.savefig(out_path)