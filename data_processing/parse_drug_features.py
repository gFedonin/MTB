from pathlib import Path

from sklearn.externals.joblib import Parallel, delayed

path_to_data = './data/drugs/txt/'
path_to_drug_list = './data/drugs.list'
out_path = './data/drug_features3_txt.csv'


def read_drug_features(drug):
    res = set()
    file = Path(path_to_data + drug + '.xml')
    if file.exists():
        with open(path_to_data + drug + '.xml') as f:
            for line in f.readlines():
                i = line.find('<Name>')
                if i != -1:
                    j = line.find('</Name>', i + 1)
                    res.add(line[i + 6: j])
    else:
        with open(path_to_data + drug + '.txt') as f:
            for line in f.readlines():
                s = line[:-1].split(' ')
                res.add(' '.join(s[2:]))
    return drug, res


if __name__ == '__main__':
    with open(path_to_drug_list, 'r') as f:
        drug_list = f.readlines()
    drug_list = [x[:-1] for x in drug_list]
    drug_features = Parallel(n_jobs=-1)(delayed(read_drug_features)(drug) for drug in drug_list)
    all_features = set()
    for drug, fset in drug_features:
        for f in fset:
            all_features.add(f)
    filtered_features = []
    for feature in all_features:
        count = 0
        for drug, fset in drug_features:
            if feature in fset:
                count += 1
        if 1 < count < len(drug_list) - 1:
            filtered_features.append(feature)

    features_to_drugs = {}
    for feature in filtered_features:
        features_to_drugs[feature] = []
    for feature in filtered_features:
        for drug, fset in drug_features:
            if feature in fset:
                features_to_drugs[feature].append(drug)

    # filtered_features2 = []
    drug_str_to_features = {}
    for feature in filtered_features:
        drugs = features_to_drugs[feature]
        drug_str = ''.join(drugs)
        # if drug_str not in drug_str_to_features:
            # filtered_features2.append(feature)
        if drug_str in drug_str_to_features.keys():
            drug_str_to_features[drug_str].append(feature)
        else:
            f_set = [feature]
            drug_str_to_features[drug_str] = f_set

    # with open(out_path, 'w') as f:
    #     drugs = ''
    #     for drug, fset in drug_features:
    #         drugs += '\t' + drug
    #     drugs += '\n'
    #     f.write(drugs)
    #     for feature in filtered_features2:
    #         f.write(feature)
    #         for drug, fset in drug_features:
    #             if feature in fset:
    #                 f.write('\t1')
    #             else:
    #                 f.write('\t0')
    #         f.write('\n')
    with open(out_path, 'w') as f:
        drugs = ''
        for drug in drug_list:
            drugs += '\t' + drug
        drugs += '\n'
        f.write(drugs)
        for drug_str, f_set in drug_str_to_features.items():
            for feature in f_set:
                f.write(feature)
                for drug in drug_list:
                    if drug in drug_str:
                        f.write('\t1')
                    else:
                        f.write('\t0')
                f.write('\n')
            f.write('\n')
