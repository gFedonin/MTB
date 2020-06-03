import gzip
from os import listdir, makedirs
from os.path import exists
from subprocess import check_call, call
import matplotlib.pyplot as plt
import pysam
from Bio import SeqIO
from core.constants import data_path, ref_len
from core.data_reading import path_to_ref, read_h37rv
from data_processing.translate_mutations import Variant
from joblib import Parallel, delayed

from data_processing.parse_gatk_variants_keep_complex import in_filter_interval, check_var, no_filters, min_dp, \
    min_alt_frac, get_intervals_to_filter_out

path_to_minimap = '/export/home/fedonin/minimap2-2.17_x64-linux/minimap2'
path_to_h37rv = data_path + 'h37rv.fasta'
path_to_index = data_path + 'h37rv_rev'
path_to_piccard = '/export/home/fedonin/picard.jar'
path_to_GATK = '/export/home/fedonin/gatk-4.1.4.0/gatk'

path_to_asm_list = data_path + 'ref_seqs_no_gaps.list'
path_to_assemblies = data_path + 'ref_seqs/'
path_to_assemblies_rev = data_path + 'ref_seqs_rev/'
path_to_mappings = data_path + 'ref_seqs_mapping/'
path_to_art = '/export/home/fedonin/art_bin_MountRainier/'
path_to_simulated_reads = data_path + 'simulated/'

dataset = 'coll18'
path_to_ids = data_path + 'snps/combined_raw_variants_mq40_keep_complex_std_names_filtered/samples_filtered.list'
path_to_meta = '../../res/new_dataset_stats/merged.metadata'
# path_to_ids = data_path + 'all_with_pheno.txt'
# path_to_raw = data_path + dataset + '/' + dataset + '_trimmed/'
path_to_raw = '/export/data/kkuleshov/myc/sra/'
# stat_path = data_path + 'sample_stats/' + dataset + '_sample_stats.csv'
stat_path = data_path + 'sample_stats/' + 'base_sample_stats.csv'

read_len = ['50', '75', '100', '150', '200', '250']
# read_len = ['50']
coverage = ['20', '30', '45', '75', '100', '125', '150', '200']
read_len_to_ss = {'50': 'NS50', '75': 'NS50', '100': 'HS20', '150': 'HS25', '200': 'MSv3', '250': 'MSv3'}

thread_num = 32
# read_num = 10000


def parse_meta():
    f = open(path_to_meta)
    s = f.readline().strip().split(',')
    bio_sample_index = 0
    for bio_sample_index in range(len(s)):
        if s[bio_sample_index] == 'BioSample':
            break
    err_to_samea = {}
    for l in f.readlines():
        s = l.strip().split(',')
        if len(s) > 1:
            err_to_samea[s[0]] = s[bio_sample_index]
    f.close()
    return err_to_samea


def rev_compl():
    for fname in listdir(path_to_assemblies):
        if fname.endswith('.fasta'):
            asm = [r for r in SeqIO.parse(path_to_assemblies + fname, 'fasta')][0]
            asm_rev = asm.reverse_complement()
            with open(path_to_assemblies + fname[:-6] + '_rev.fasta', 'w') as f:
                f.write('>' + asm.id + '\n')
                f.write(str(asm_rev.seq) + '\n')


def map_assembly(asm_name):
    cmd = path_to_minimap + ' -ax asm5 ' + path_to_ref + ' ' + path_to_assemblies_rev + '\'' + asm_name + '\'' + \
          '_rev.fasta | samtools sort - -o ' + path_to_mappings + \
          asm_name.replace(' ', '_') + '.bam'
    check_call(cmd, shell=True)
    check_call('samtools index ' + path_to_mappings + asm_name.replace(' ', '_') + '.bam', cwd=path_to_mappings,
               shell=True)


def map_all():
    # rev_compl()
    for fname in listdir(path_to_assemblies):
        if fname.endswith('.fasta'):
            map_assembly(fname[:-6])


def parse_bam(asm_name, ref_seq):
    assembly = [r for r in SeqIO.parse(path_to_assemblies_rev + asm_name + '_rev.fasta', 'fasta')][0].seq
    samfile = pysam.AlignmentFile((path_to_mappings + asm_name + ".bam").replace(' ', '_'), "rb")
    raw_variants = []
    for read in samfile.fetch():
        ref_pos = read.reference_start
        read_pos = read.query_alignment_start
        # print('ref_pos = %d, contig_pos =% d' % (ref_pos, read_pos))
        # print(read.cigarstring)
        for op, l in read.cigartuples:
            if op == 0:
                # match
                # print(ref_seq[ref_pos: ref_pos + l])
                # print(mut_str[read_pos: read_pos + l])
                for i in range(l):
                    if ref_seq[ref_pos] != assembly[read_pos]:
                        raw_variants.append(Variant(ref_pos + 1, ref_seq[ref_pos], assembly[read_pos]))
                    ref_pos += 1
                    read_pos += 1
            elif op == 1:
                # insert
                raw_variants.append(Variant(ref_pos + 1, '', assembly[read_pos: read_pos + l]))
                read_pos += l
            elif op == 2:
                # delete
                # raw_variants.append(Variant(ref_pos, ref_seq[ref_pos: ref_pos + l], ''))
                for i in range(l):
                    raw_variants.append(Variant(ref_pos + i + 1, ref_seq[ref_pos + i], ''))
                ref_pos += l
            elif op == 5:
                # hardclip
                read_pos += l
    samfile.close()
    raw_variants.sort()
    with open(path_to_mappings + asm_name + '.variants', 'w') as f:
        for var in raw_variants:
            f.write('%d\t%s\t%s\n' % (var.pos, var.ref, var.alt))


def true_variants():
    ref_seq = read_h37rv()
    for fname in listdir(path_to_assemblies):
        if fname.endswith('.fasta'):
            parse_bam(fname[:-6], ref_seq)


def simulate_reads(asm_name, rl, cov):
    if not exists(path_to_simulated_reads + rl + '/' + cov + '/\"' + asm_name + '1.fq.gz\"'):
        insert_size = 2 * int(rl) + 50
        insert_std = '10'
        cmd = path_to_art + 'art_illumina -na -p -q -ss ' + read_len_to_ss[rl] + ' -i \"' + \
              path_to_assemblies + asm_name + '.fasta\" -l ' + rl + ' -f ' + cov + \
              ' -m ' + str(insert_size) + ' -s ' + insert_std + \
              ' -o \"' + path_to_simulated_reads + rl + '/' + cov + '/' + asm_name + '_\"'
        check_call(cmd, shell=True)
        cmd = 'gzip \"' + path_to_simulated_reads + rl + '/' + cov + '/' + asm_name + '_1.fq\"'
        check_call(cmd, shell=True)
        cmd = 'gzip \"' + path_to_simulated_reads + rl + '/' + cov + '/' + asm_name + '_2.fq\"'
        check_call(cmd, shell=True)
        return 1
    else:
        return 0


def simulate_all():
    for rl in read_len:
        if not exists(path_to_simulated_reads + rl):
            makedirs(path_to_simulated_reads + rl)
        for cov in coverage:
            if not exists(path_to_simulated_reads + rl + '/' + cov):
                makedirs(path_to_simulated_reads + rl + '/' + cov)
            tasks = Parallel(n_jobs=-1)(delayed(simulate_reads)(fname[:-6], rl, cov)
                                        for fname in listdir(path_to_assemblies) if fname.endswith('.fasta'))
            c = 0
            for task in tasks:
                c += task
            print('rl=%s cov=%s n_files=%d' % (rl, cov, c))
            # for fname in listdir(path_to_assemblies):
            #     if fname.endswith('.fasta'):
            #         simulate_reads(fname[:-6], rl, cov)


def sample_stat():
    sample_ids = {l.strip() for l in open(path_to_ids).readlines()}
    with open(stat_path, 'w') as fout:
        for sample_id in listdir(path_to_raw):
            if sample_id not in sample_ids:
                continue
            av_read_len = 0
            read_num = 0
            # if exists(path_to_raw + sample_id + '/' + sample_id + '_p1.fastq.gz'):
            if exists(path_to_raw + sample_id + '/' + sample_id + '_R1.fastq.gz'):
                try:
                    # f = gzip.open(path_to_raw + sample_id + '/' + sample_id + '_p1.fastq.gz')
                    f = gzip.open(path_to_raw + sample_id + '/' + sample_id + '_R1.fastq.gz')
                    lines = f.readlines()
                    for i in range(1, len(lines), 4):
                        av_read_len += len(lines[i])
                        read_num += 1
                    f.close()
                    fout.write(sample_id + '\t' + str(av_read_len / read_num) + '\t' + str(2 * read_num) + '\n')
                except:
                    fout.write(sample_id + '\t0\t0\n')


path_to_sample_stats = '../../res/new_dataset_stats/sample_stats/'
path_to_hist_rl = '../../res/new_dataset_stats/rl.png'


def read_len_hist():
    read_lens = []
    for fname in listdir(path_to_sample_stats):
        for l in open(path_to_sample_stats + fname).readlines():
            s = l.strip().split()
            read_lens.append(float(s[1]))
    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
    fig.suptitle('Av read len', fontsize=20)
    plt.xlabel('read len', fontsize=18)
    plt.ylabel('sample num', fontsize=16)
    axs.hist(read_lens, bins=20)
    plt.savefig(path_to_hist_rl)


path_to_hist_rn = '../../res/new_dataset_stats/rn.png'
path_to_depths = '../../res/new_dataset_stats/rn.stat'
path_to_filtered_out_err = '../../res/new_dataset_stats/no_such_err.list'
path_to_filtered_out_samea = '../../res/new_dataset_stats/no_such_samea.list'


def read_num_hist():
    err_to_samea = parse_meta()
    filtered_out_samea = []
    filtered_out_err = []
    read_nums = []
    depths = []
    sample_ids = {l.strip() for l in open(path_to_ids).readlines()}
    for fname in listdir(path_to_sample_stats):
        for l in open(path_to_sample_stats + fname).readlines():
            s = l.strip().split()
            biosample = err_to_samea.get(s[0])
            if biosample is not None and biosample not in sample_ids:
                filtered_out_err.append(s[0])
                filtered_out_samea.append(biosample)
                # print(s[0])
                # print(biosample)
                continue
            read_nums.append(float(s[2]) * float(s[1]) / ref_len)
            depths.append((s[0], float(s[2]) * float(s[1]) / ref_len))
    fig, axs = plt.subplots(1, 1, tight_layout=True)
    fig.suptitle('Sequencing depth', fontsize=20)
    plt.xlabel('depth', fontsize=18)
    plt.ylabel('sample num', fontsize=16)
    read_nums.sort()
    axs.hist(read_nums[20:-20], bins=40)
    plt.savefig(path_to_hist_rn)
    depths.sort(key=lambda x: x[1])
    with open(path_to_depths, 'w') as f:
        for sid, d in depths:
            f.write(sid + '\t' + str(d) + '\n')
    with open(path_to_filtered_out_err, 'w') as f:
        f.write('\n'.join(filtered_out_err))
    with open(path_to_filtered_out_samea, 'w') as f:
        f.write('\n'.join(filtered_out_samea))


def bwa_mem(id, rl, cov):
    thread_num = 10
    if not exists(path_to_simulated_reads + rl + '/' + cov + '/' + id + '.bam'):
        check_call('bwa mem -t ' + str(thread_num) + ' -R \'@RG\\tID:1\\tSM:' + id
                   + '\' ' + path_to_index + ' \"' + path_to_simulated_reads + rl + '/' + cov + '/' + id + '_1.fq.gz\" \"' +
                   path_to_simulated_reads + rl + '/' + cov + '/' + id + '_2.fq.gz\" | samtools sort -@ ' + str(thread_num) +
                   ' - -o \"' + path_to_simulated_reads + rl + '/' + cov + '/' + id + '.bam\"',
                   shell=True)
        return 1
    else:
        return 0


def map_all_reads():
    tasks = Parallel(n_jobs=14)(delayed(bwa_mem)(fname[:-6], rl, cov) for rl in read_len for cov in coverage
                                        for fname in listdir(path_to_assemblies) if fname.endswith('.fasta'))
    c = 0
    for task in tasks:
        c += task
    print('n_files=%d' % c)


def mark_dup(id, rl, cov):
    c = 'java -Xmx4g -jar ' + path_to_piccard + ' MarkDuplicates INPUT=\"' + path_to_simulated_reads + rl + '/' + cov + \
        '/' + id + '.bam\" OUTPUT=\"' + path_to_simulated_reads + rl + '/' + cov + '/' + id + '_mark_dup.bam\"' + \
        ' METRICS_FILE=\"' + path_to_simulated_reads + rl + '/' + cov + '/' + id + \
        '.metrix\" REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT' + \
        ' > \"' + path_to_simulated_reads + rl + '/' + cov + '/' + id + '.log\" 2>> \"' + path_to_simulated_reads + rl + \
        '/' + cov + '/' + id + '.log\"'
    call(c, shell=True)
    check_call('samtools index \"' + path_to_simulated_reads + rl + '/' + cov + '/' + id + '_mark_dup.bam\"', shell=True)
    return 1


def mark_all_dup():
    for rl in read_len:
        for cov in coverage:
            tasks = Parallel(n_jobs=10)(delayed(mark_dup)(fname[:-6], rl, cov)
                                        for fname in listdir(path_to_assemblies) if fname.endswith('.fasta'))
            c = 0
            for task in tasks:
                c += task
            print('n_files=%d' % c)


def snp_call_gvcf(id, rl, cov):
    if exists(path_to_simulated_reads + rl + '/' + cov + '/' + id + '_mark_dup.bam'):
        if not exists(path_to_simulated_reads + rl + '/' + cov + '/' + id + '_h37rv.g.vcf'):
            cln = path_to_GATK + ' --java-options "-Xmx1g" HaplotypeCaller -VS SILENT --native-pair-hmm-threads 2 ' \
                                 '-R ' + path_to_h37rv + ' -I \"' + \
                           path_to_simulated_reads + rl + '/' + cov + '/' + id + '_mark_dup.bam\" -O \"' + \
                  path_to_simulated_reads + rl + '/' + cov + '/' + id + \
                           '_h37rv.g.vcf\" -ERC GVCF --smith-waterman FASTEST_AVAILABLE -ploidy 1 > \"' + \
                  path_to_simulated_reads + rl + '/' + cov + '/' + id + '_hc.log\" 2>> \"' + path_to_simulated_reads + \
                  rl + '/' + cov + '/' + id + '_hc.log\"'
            call(cln, shell=True)  # --native-pair-hmm-threads 4
            return 1
        else:
            return 0
    else:
        print(path_to_simulated_reads + rl + '/' + cov + '/' + id + '_mark_dup.bam not found')
        return 0


def snp_call_all():
    tasks = Parallel(n_jobs=thread_num//2)(delayed(snp_call_gvcf)(fname[:-6], rl, cov)
                                           for rl in read_len for cov in coverage
                                for fname in listdir(path_to_assemblies) if fname.endswith('.fasta'))
    c = 0
    for task in tasks:
        c += task
    print('n_files=%d' % c)


def gen_map_file():
    with open(path_to_simulated_reads + 'simulated.map', 'w') as f:
        for rl in read_len:
            for cov in coverage:
                for fname in open(path_to_asm_list).readlines():
                    fname = fname.strip()
                    f.write(fname + '_' + rl + '_' + cov + '\t' + path_to_simulated_reads + rl +
                            '/' + cov + '/' + fname + '_h37rv.g.vcf\n')


def run_GenomicsDBImport():
    temp_path = path_to_simulated_reads + 'tmp/'
    log_path_db = path_to_simulated_reads + 'db.log'
    path_to_map_file = path_to_simulated_reads + 'simulated.map'
    path_to_db = path_to_simulated_reads + 'simulated.db'
    if not exists(temp_path):
        makedirs(temp_path)
    call(path_to_GATK + ' --java-options "-Xmx4g -Xms4g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenomicsDBImport ' + \
        '--genomicsdb-workspace-path ' + path_to_db + ' --batch-size 1000 -L ch1:1-' + \
         str(ref_len) + ' --sample-name-map ' + path_to_map_file + ' --tmp-dir=' + temp_path + ' --reader-threads ' + \
         str(thread_num) + ' > ' + log_path_db + ' 2>> ' + log_path_db, shell=True)#


def run_GenotypeGVCFs():
    path_to_db = path_to_simulated_reads + 'simulated.db'
    merged_vcf_path = path_to_simulated_reads + 'merged.vcf'
    log_path_genotype = path_to_simulated_reads + 'genotype.log'
    temp_path = path_to_simulated_reads + 'tmp/'
    call(path_to_GATK + ' --java-options "-Xmx4g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenotypeGVCFs ' +
        '--sample-ploidy 1 -R ' + path_to_h37rv + ' -V gendb://' + path_to_db + ' -O ' + merged_vcf_path + ' ' +
         '--tmp-dir=' + temp_path + ' > ' +
    log_path_genotype + ' 2>> ' + log_path_genotype, shell=True)


def extract_sample_gatk(sample_id, out_path_vcf, big_vcf, rl, cov):
    if not exists(out_path_vcf + sample_id + '.vcf'):
        call(path_to_GATK + ' --java-options "-Xmx1g -Xms1g" SelectVariants --exclude-non-variants true ' +
                              '--remove-unused-alternates true -OVI false -R ' + path_to_h37rv + ' -sn \"' +
               sample_id + '_' + rl + '_' + cov + '\" -V ' + big_vcf + ' -O \"' + out_path_vcf + sample_id + '.vcf\"', shell=True)
        return 1
    return 0


def split_gatk():
    tasks = Parallel(n_jobs=32)(delayed(extract_sample_gatk)(fname.strip(),
                                                             path_to_simulated_reads + rl + '/' + cov + '/',
                                                             path_to_simulated_reads + 'merged.vcf', rl, cov)
                                for rl in read_len for cov in coverage
                                for fname in open(path_to_asm_list).readlines())
    c = 0
    for task in tasks:
        c += task
    print('%d samples extracted' % c)


def parse(id, rl, cov, filter_intervals):
    filter_intervals_starts = [x[0] for x in filter_intervals]
    with open(path_to_simulated_reads + rl + '/' + cov + '/' + id + '.variants', 'w') as fout:
        var_list = []
        ann_dict = {}
        f = open(path_to_simulated_reads + rl + '/' + cov + '/' + id + '.vcf')
        lines = f.readlines()
        fout.write('# source ' + path_to_simulated_reads + rl + '/' + cov + '/' + id + '\n')
        for l in lines:
            if l.startswith('#'):
                continue
            s = l.strip().split('\t')
            annotations = s[7].split(';')
            variants = s[4].split(',')
            stats = s[-1].split(':')
            # stat_names = s[-2].split(':')
            read_counts = [int(c) for c in stats[1].split(',')]
            total_depth = sum(read_counts)
            if total_depth < min_dp:
                continue
            if len(variants) == 1:
                if read_counts[-1]/total_depth >= min_alt_frac and read_counts[-1] >= min_dp:
                    var_list.append((int(s[1]), s[3], s[4], annotations))
            else:
                scores = stats[-1].split(',')
                for i in range(1, len(scores)):
                    if scores[i] == '0':
                        if read_counts[i]/total_depth >= min_alt_frac and read_counts[i] >= min_dp:
                            var_list.append((int(s[1]), s[3], variants[i - 1], annotations))
                        break
            for a in annotations:
                if ',' in a:
                    continue
                i = a.find('=')
                vals = ann_dict.get(a[0:i])
                if vals is None:
                    vals = []
                    ann_dict[a[0:i]] = vals
                vals.append(float(a[i + 1:]))
        f.close()
        res = []
        for pos, ref, alt, ann in var_list:
            if (no_filters or check_var(ann)) and not in_filter_interval(pos, filter_intervals_starts, filter_intervals):
                res.append((pos, ref, alt))
        var_list = list(res)
        var_list.sort(key=lambda x: x[0])
        for pos, alt, t in var_list:
            fout.write('%d\t%s\t%s\n' % (pos, alt, t))
    return 1


def parse_all():
    filter_intervals = get_intervals_to_filter_out()
    tasks = Parallel(n_jobs=10)(delayed(parse)(fname.strip(), rl, cov, filter_intervals)
                                                                for rl in read_len for cov in coverage
                                                                    for fname in open(path_to_asm_list).readlines())
    c = 0
    for task in tasks:
        c += task
    print('%d samples processed' % c)


def rename_all():
    for rl in read_len:
        for cov in coverage:
            for fname in open(path_to_asm_list).readlines():
                fname = fname.strip()
                no_spaces = fname.replace(' ', '_')
                if no_spaces != fname:
                    # call('mv \"' + path_to_simulated_reads + rl + '/' + cov + '/' + fname + '_h37rv.g.vcf\" ' +
                    #      path_to_simulated_reads + rl + '/' + cov + '/' + no_spaces + '_h37rv.g.vcf', shell=True)
                    # call('mv \"' + path_to_simulated_reads + rl + '/' + cov + '/' + fname + '_h37rv.g.vcf.idx\" ' +
                    #      path_to_simulated_reads + rl + '/' + cov + '/' + no_spaces + '_h37rv.g.vcf.idx', shell=True)
                    call('mv \"' + path_to_simulated_reads + rl + '/' + cov + '/' + fname + '_h37rv.g.vcf.idx\" ' +
                         path_to_simulated_reads + rl + '/' + cov + '/' + no_spaces + '_h37rv.g.vcf.idx', shell=True)


if __name__ == '__main__':
    # map_all()
    # simulate_all()
    # sample_stat()
    # read_len_hist()
    # read_num_hist()
    # map_all_reads()
    # mark_all_dup()
    # rename_all()
    # gen_map_file()
    # snp_call_all()
    # true_variants()
    # run_GenomicsDBImport()
    # run_GenotypeGVCFs()
    split_gatk()