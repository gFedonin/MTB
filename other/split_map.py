from core.constants import data_path

path_to_big_map = data_path + 'combined_ids_nocoll.map'
batch_size = 2000
prefix = 'combined_no_coll'


if __name__ == '__main__':
    sample_ids = []
    lines = []
    for l in open(path_to_big_map).readlines():
        lines.append(l)
        sample_ids.append(l.split('\t')[0])
    for i in range(len(sample_ids)//batch_size):
        with open(data_path + prefix + str(i + 1) + '.list', 'w') as f:
            if i != len(sample_ids)//batch_size - 1:
                f.write('\n'.join(sample_ids[i*batch_size: (i + 1)*batch_size]))
            else:
                f.write('\n'.join(sample_ids[i * batch_size: len(sample_ids)]))
        with open(data_path + prefix + str(i + 1) + '.map', 'w') as f:
            if i != len(sample_ids)//batch_size - 1:
                f.write(''.join(lines[i*batch_size: (i + 1)*batch_size]))
            else:
                f.write(''.join(lines[i * batch_size: len(sample_ids)]))