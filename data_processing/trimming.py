from os import system, makedirs, listdir
from os.path import exists
from subprocess import check_call, call

from sklearn.externals.joblib import Parallel, delayed

from core.constants import data_path

data_set = 'patric'
dataFolder = data_path + data_set + '/' + data_set + '_raw/'
# dataFolder = data_path + 'coll18/coll18_raw/'
# dataFolder = '/export/data/kkuleshov/myc/sra/'
# outFolder = data_path + 'walker18/walker18_trimmed/'
outFolder = data_path + data_set + '/' + data_set + '_trimmed/'
# outFolder = data_path + data_set + '/' + 'test/'
# outFolder = data_path + 'coll18/coll18_trimmed_bbduk_q30/'
# outFolder = data_path + 'kkuleshov_bbduk_q30/'
# listPath = data_path + 'walker18/downloaded.list'
listPath = data_path + data_set + '/' + data_set+ '_new.samples'
# listPath = data_path + data_set + '/Genoscreen_err.list'
# listPath = data_path + 'coll18/coll18_supp.samples'
# listPath = data_path + 'all_with_pheno.txt'
pathToAdapters = "/export/home/fedonin/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa"
pathToAdaptersS = "/export/home/fedonin/Trimmomatic-0.38/adapters/TruSeq3-SE-2.fa"
# pathToAdapters = "/export/home/fedonin/bbmap/resources/adapters.fa"
pathToTrimmomatic = "/export/home/fedonin/Trimmomatic-0.38/trimmomatic-0.38.jar"
path_to_bbtools = '/export/home/fedonin/bbmap/'
threadNum = '4'
sample_in_parallel = 8
suffix1 = "_1.fastq.gz"
suffix2 = "_2.fastq.gz"
suffixS = '.fastq.gz'
# suffix1 = "_R1.fastq.gz"
# suffix2 = "_R2.fastq.gz"

clipQ30 = True

overwrite = False


def trim_sample(name):
    if exists(outFolder + name):
        call('rm ' + name + '*', shell=True, cwd=outFolder + name)
    else:
        makedirs(outFolder + name)
    if exists(dataFolder + name + suffix1) and exists(dataFolder + name + suffix2):
        if clipQ30:
            call("java -jar " + pathToTrimmomatic +
                   " PE -threads " + threadNum + " -phred33" +
                   " -trimlog " + outFolder + name + '/' + name + ".log " + dataFolder + name + suffix1 + " " +
                   dataFolder + name + suffix2 + " " + outFolder + name + '/' + name + "_p1.fastq.gz " +
                   outFolder + name + '/' + name + "_u1.fastq.gz " + outFolder + name + '/' + name + "_p2.fastq.gz " +
                   outFolder + name + '/' + name + "_u2.fastq.gz " + " ILLUMINACLIP:" + pathToAdapters +
                   ":1:30:10 SLIDINGWINDOW:10:30", shell=True)
        else:
            call("java -jar " + pathToTrimmomatic +
                   " PE -threads " + threadNum + " -phred33" +
                   " -trimlog " + outFolder + name + '/' + name + ".log " + dataFolder + name + suffix1 + " " +
                   dataFolder + name + suffix2 + " " + outFolder + name + '/' + name + "_p1.fastq.gz " +
                   outFolder + name + '/' + name + "_u1.fastq.gz " + outFolder + name + '/' + name + "_p2.fastq.gz " +
                   outFolder + name + '/' + name + "_u2.fastq.gz " + " ILLUMINACLIP:" + pathToAdapters + ":1:30:10",
                 shell=True)
        return 1
    elif exists(dataFolder + name + suffixS):
        if clipQ30:
            system("java -jar " + pathToTrimmomatic + " SE -threads " + threadNum + " -phred33" +
                   " -trimlog " + outFolder + name + '/' + name + ".log " + dataFolder + name + suffixS + " " +
                   outFolder + name + '/' + name + "_S.fastq.gz " + " ILLUMINACLIP:" + pathToAdapters +
                   ":1:30:10 SLIDINGWINDOW:10:30")
        else:
            system("java -jar " + pathToTrimmomatic + " SE -threads " + threadNum + " -phred33" +
                   " -trimlog " + outFolder + name + '/' + name + ".log " + dataFolder + name + suffixS + " " +
                   outFolder + name + '/' + name + "_S.fastq.gz " + " ILLUMINACLIP:" + pathToAdapters + ":1:30:10")
        return 1
    else:
        return 0


def trim_all():
    if not exists(outFolder):
        makedirs(outFolder)
    trimmed = set(l for l in listdir(outFolder))
    raw = set(l for l in listdir(dataFolder))
    samples_to_trim = []
    for l in open(listPath).readlines():
        sample_id = l.strip()
        if sample_id + suffix1 in raw:
            if sample_id not in trimmed:
                samples_to_trim.append(sample_id)
            else:
                if not exists(outFolder + sample_id + '/' + sample_id + "_p1.fastq.gz") or not exists(outFolder +
                                                                    sample_id + '/' + sample_id + "_p2.fastq.gz"):
                    samples_to_trim.append(sample_id)
    tasks = Parallel(n_jobs=sample_in_parallel)(delayed(trim_sample)(sample_id) for sample_id in samples_to_trim)
    c = 0
    for task in tasks:
        c += task
    print('%d samples processed' % c)


def trim_bbduk():
    for l in open(listPath).readlines():
        name = l.strip()
        r1 = dataFolder + name + suffix1
        r2 = dataFolder + name + suffix2
        # r1 = dataFolder + name + '/' + name + suffix1
        # r2 = dataFolder + name + '/' + name + suffix2
        if not exists(r1) or not exists(r2):
            continue
        if not exists(outFolder + name):
            makedirs(outFolder + name)
        r1_save = outFolder + name + '/' + name + "_p1.fastq.gz"
        r2_save = outFolder + name + '/' + name + "_p2.fastq.gz"
        if not overwrite and exists(r1_save) and exists(r2_save):
            continue
        stat_path = outFolder + name + '/' + name + '.stat'
        log_path = outFolder + name + '/' + name + '.log'
        check_call("{}bbduk.sh -Xmx1g tossbrokenreads=t in1={} in2={} out1={} out2={} ref={} stats={} ktrim=r k=23 mink=11 hdist=1 " \
        "tpe tbo trimq=30 qtrim=r > {}".format(
            path_to_bbtools, r1, r2, r1_save, r2_save, pathToAdapters, stat_path, log_path), shell=True)


if __name__ == '__main__':
    if not exists(outFolder):
        makedirs(outFolder)
    # trim_bbduk()
    trim_all()
