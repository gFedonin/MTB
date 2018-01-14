import os

path_to_data = './data/pheno/'


if __name__ == '__main__':
    for (dirpath, dirnames, filenames) in os.walk(path_to_data):
        for filename in filenames:
            name = filename[0:-6]
            os.system('nohup ./BayesTraitsV3 ' + name + '.nexus ' + name + '.pheno < config.txt > ' + name + '.out 2>> ' + name + '.out &')
