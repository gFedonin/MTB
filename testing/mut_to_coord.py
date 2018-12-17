from src.core.annotations import read_annotations, CDSType
from src.core.constants import upstream_length

gene = 'katG'
pos = '315'


if __name__ == '__main__':
    cds_list = read_annotations(upstream_length)
    for cds in cds_list:
        if cds.name == gene and cds.type == CDSType.Gene:
            if cds.strand == 1:
                print(str(cds.start + 3*(int(pos) - 1)))
            else:
                print(str(cds.end - 3*(int(pos) - 1)))
            break