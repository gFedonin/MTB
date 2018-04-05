table_to_convert = './res/simple_full.csv'
path_to_model = './res/model_full.csv'
out_path = './res/by_drug.csv'

drug_to_model = {}
with open(path_to_model, 'r') as f:
    f.readline()
    for line in f.readlines():
        s = line.strip().split('\t')
        drug_to_model[s[0]] = s[-1]

with open(table_to_convert, 'r') as f:
    with open(out_path, 'w') as out:
        dicts = f.readline().strip().split('\t')[2:]
        for line in f.readlines():
            s = line.strip().split('\t')
            out.write(s[0] + '\n')
            out.write("Охарактеризовано образцов: " + s[1] + "\n")
            out.write("Словарь\tЧисло мутаций в словаре\tЧисло мутаций из словаря в датасете\tPPV\tNPV\tЧувствительность\tСпецифичность\n")
            for i in range(len(dicts)):
                out.write(dicts[i])
                vals = s[i + 2].split(';')
                for v in vals:
                    out.write('\t' + v)
                out.write('\n')
            out.write('Модель')
            vals = drug_to_model[s[0]].split(';')
            for v in vals:
                out.write('\t' + v)
            out.write('\n')
            out.write('\n\n')

