from core.constants import data_path

path_to_lists = [data_path + 'base_coll18_walker18_johnsen19_farhat19.samples',
                 data_path + 'walker18/walker18_samea.samples',
                 data_path + 'walker18/samea.list', data_path + 'yang17/yang17_converted.samples']
merged_list_path = data_path + 'combined_samea.list'
list_to_filter_path = data_path + 'patric/patric.list'
list_to_filter_mtb_path = data_path + 'patric/patric_mtb_samea.list'
filtered_list_path = data_path + 'patric/patric_new.samples'
metadata_path = data_path + 'patric/metadata.merged'
filtered_run_ids_path = data_path + 'patric/patric_new.list'


def filter_non_mtb():
    run_to_meta = {}
    sample_to_meta = {}
    samea_to_meta = {}
    for l in open(metadata_path).readlines():
        s = l.strip().split(',')
        if len(s) > 25 and s[0] != 'Run':
            run_to_meta[s[0]] = s
            samea_to_meta[s[25]] = s
            sample_to_meta[s[24]] = s
    with open(list_to_filter_mtb_path, 'w') as f:
        for l in open(list_to_filter_path).readlines():
            sid = l.strip()
            if sid.startswith('SAM'):
                meta = samea_to_meta.get(sid)
                if meta is not None:
                    if 'Mycobacterium tuberculosis' in meta[28]:
                        f.write(sid + '\n')
            elif sid.startswith('SRS'):
                meta = sample_to_meta[sid]
                if 'Mycobacterium tuberculosis' in meta[28]:
                    f.write(meta[25] + '\n')
            else:
                print(sid)


def merge_lists():
    sample_ids = set()
    for file in path_to_lists:
        for l in open(file).readlines():
            sample_ids.add(l.strip())
    with open(merged_list_path, 'w') as f:
        f.write('\n'.join(sorted(list(sample_ids))))


def filter_by_list():
    old_samples = {l.strip() for l in open(merged_list_path).readlines()}
    with open(filtered_list_path, 'w') as f:
        for l in open(list_to_filter_mtb_path).readlines():
            if l.strip() not in old_samples:
                f.write(l.strip() + '\n')


def convert_samea_to_run():
    samea_to_run = {}
    for l in open(metadata_path).readlines():
        s = l.strip().split(',')
        if len(s) > 25 and s[0] != 'Run':
            if s[25] in samea_to_run:
                samea_to_run[s[25]].append(s[0])
            else:
                samea_to_run[s[25]] = [s[0]]
    with open(filtered_run_ids_path, 'w') as f:
        for l in open(filtered_list_path).readlines():
            f.write('\n'.join(samea_to_run[l.strip()]) + '\n')


if __name__ == '__main__':
    # merge_lists()
    # filter_non_mtb()
    # filter_by_list()
    convert_samea_to_run()