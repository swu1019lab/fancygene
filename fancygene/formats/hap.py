# -*- coding: utf-8 -*-
# @Time    : 2024/8/21 19:41
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : hap.py

import pandas as pd

def read_hap(file_path, keep_hap: list = None, rename_hap: bool = False, hap_map: dict = None):
    """
    Read the haplotype file

    :param file_path: str, the path of the haplotype file
    :param keep_hap: list, the haplotypes to be kept, default is None
    :param rename_hap: bool, whether to rename the haplotypes, default is False
    :param hap_map: dict, the mapping of the haplotypes, need to be specified when rename_hap is True, default is None
    :return: dict, the haplotype data including the haplotype group, snp information and genotype data
    """
    hap_data = dict()
    df = pd.read_csv(file_path, header=[0, 1, 2, 3], index_col=0)
    # only keep the specified haplotypes
    if keep_hap:
        df = df[df.loc[:, ('haplotypes', slice(None), slice(None), slice(None))].iloc[:, 0].isin(keep_hap)]
    df = df.sort_values(by=df.columns[-1], axis=0)
    # get the haplotype group
    hap_data['group'] = df.iloc[:, -1].reset_index().droplevel([1, 2, 3], axis=1)
    # rename the haplotype names
    if rename_hap and hap_map:
        hap_data['group'] = hap_data['group'].replace(hap_map)
    # get the snp information
    hap_data['chrom'] = df.columns.get_level_values(0).tolist()
    hap_data['pos'] = df.columns.get_level_values(1).tolist()
    hap_data['ref'] = df.columns.get_level_values(2).tolist()
    hap_data['alt'] = df.columns.get_level_values(3).tolist()
    # get the genotype data
    hap_data['geno'] = df.iloc[:, :-1].values
    return hap_data
