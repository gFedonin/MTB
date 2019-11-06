from bisect import bisect_left

from sklearn.externals.joblib import Parallel, delayed

import pysam as ps
from core.annotations import read_annotations
from core.constants import data_path, upstream_length

path_to_ids = data_path + 'all_with_pheno_and_snp.txt'#'test.list'
path_to_bam = '/export/data/kchukreev/data/bam_files/gatk/realigned_bams_with_known/'
path_to_short_tandem_repeats = data_path + 'h37rv.fasta.2.7.7.80.10.20.10.dat'
out_path = data_path + 'filtered_bams/'

threshold = 5


def read_repeats():
    coords = []
    with open(path_to_short_tandem_repeats) as f:
        for line in f.readlines()[15:]:
            s = line.split()
            coords.append((int(s[0]), int(s[1]) + 1))
    coords.sort(key=lambda tup: tup[0])
    return coords


def filter_bam(sample_id, coords):
    coord_starts = [s for s,e in coords]
    samfile = ps.AlignmentFile(path_to_bam + sample_id + '_h37rv.bam', "rb")
    with ps.AlignmentFile(out_path + sample_id + '_h37rv.bam', "wb", header=samfile.header) as outf:
        for read in samfile.fetch():
            i = bisect_left(coord_starts, read.reference_start)
            keep_this_read = True
            if i != 0:
                # check the left one
                if read.reference_start < coords[i - 1][1] + threshold:
                    # this read intersects the interval
                    # interval is not inside the read alignment
                    keep_this_read = False
            if i != len(coords):
                # check the right ones
                for coord in coords[i:]:
                    if read.reference_end > coord[0] - threshold:
                        # this read intersects the interval
                        if read.reference_start > coord[0] - threshold:
                            # interval is not inside the read alignment
                            keep_this_read = False
                            break
                        if read.reference_end < coord[1] + threshold:
                            # interval is not inside the read alignment
                            keep_this_read = False
                            break
                    else:
                        break
            if keep_this_read:
                outf.write(read)
    return 1


if __name__ == '__main__':
    coords = read_repeats()
    tasks = Parallel(n_jobs=-1)(delayed(filter_bam)(sample_id.strip(), coords) for sample_id in open(path_to_ids).readlines())
    c = 0
    for task in tasks:
        c += task
    print('%d bam files filtered\n' % c)
