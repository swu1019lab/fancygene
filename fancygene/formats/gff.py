# -*- coding: utf-8 -*-
# @Time    : 2024/8/21 19:24
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : gff.py

import re
import gzip
from collections import defaultdict

def extract_coord(file_path, file_format='gff3', feature='gene', key='ID'):
    """
    Extract the coordinates of the specified feature from the gff file
    :param file_path: str, the path of the gff file
    :param file_format: str, the format of the gff file, default is 'gff3'
    :param feature: str, the feature to be extracted, default is 'gene'
    :param key: str, the key of the feature, default is 'ID'
    :return: dict, the coordinates of the specified feature
    """
    coord = defaultdict(list)
    patterns = {
        'gff3': re.compile(r'(?P<key>\w+)=(?P<value>[^;]+)'),
        'gtf2': re.compile(r'(?P<key>\w+) "(?P<value>[^"]+)"')
    }
    pattern = patterns[file_format]
    # whether the gff/gtf file is compressed or not
    if file_path.endswith('.gz'):
        open_func = gzip.open
    else:
        open_func = open
    with open_func(file_path, 'rt') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            line = line.strip().split('\t')
            if line[2] == feature:
                attributes = line[8].split(';')
                attributes = {pattern.match(attr).group('key'): pattern.match(attr).group('value') for attr in attributes}
                if key in attributes:
                    coord[attributes[key]].append([line[0], int(line[3]), int(line[4]), line[6]])
                else:
                    print(attributes)
                    raise ValueError('The key of the feature is not found in the attributes')
    return coord
