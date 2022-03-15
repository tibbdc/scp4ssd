#!/home/zhangjq/anaconda3/envs/september/bin/python
# -*-coding:utf-8 -*-
'''
@File    :   predict.py
@Time    :   2022/02/21 17:11:24
@Author  :   zhangjq
@Version :   1.0
@Contact :   zhangjq@tib.cas.cn
@License :   (C)Copyright 2022-2023, zhangjq
@Desc    :   Enjoy your dinner
'''

import os
import re
import numpy as np
import pandas as pd
from joblib import load

from features.api import calc_feat

import argparse

def create_logger(name, silent=False, to_disk=True, log_file=None):
    """Create a new logger"""
    import logging
    import time
    from time import strftime, gmtime
    import random
    import sys
    from datetime import datetime

    # setup logger
    log = logging.getLogger(name)
    log.setLevel(logging.DEBUG)
    log.propagate = False
    formatter = logging.Formatter(fmt='%(message)s', datefmt='%Y/%m/%d %I:%M:%S')
    if not silent:
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.DEBUG)
        ch.setFormatter(formatter)
        log.addHandler(ch)
    if to_disk:
        log_file = log_file if log_file is not None else strftime("log/log_%m%d_%H%M.txt", gmtime())
        if type(log_file) == list:
            for filename in log_file:
                fh = logging.FileHandler(filename, mode='w')
                fh.setLevel(logging.INFO)
                fh.setFormatter(formatter)
                log.addHandler(fh)
        if type(log_file) == str:
            fh = logging.FileHandler(log_file, mode='w')
            fh.setLevel(logging.INFO)
            fh.setFormatter(formatter)
            log.addHandler(fh)
    return log

######################################################################################
#      python predict.py --fasta ./examples/example.fna --out example_out.csv        #
######################################################################################

def predict(data: pd.DataFrame, return_prob: bool = False):
    clf = load('./models/scp4ssd.joblib')

    prediction = np.expand_dims(clf.predict(data), axis=1)
    
    if return_prob:
        prob = clf.predict_proba(data)
        return np.concatenate((prediction, prob), axis=1)
    return prediction
    
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--fasta', '-f', help='fasta file of atcg sequences', required=True)
    parser.add_argument(
        '--out', '-o', help='output csv file name', required=True)
    parser.add_argument(
        '--verb', '-v', help='(opt) shows program progress', default=False, required=False)
    # parser.add_argument("--help", '-h', description="")

    return parser.parse_args()

def main():
    # read args
    args = parse_args()

    # calculate nucleotide sequence features
    df = calc_feat(args.fasta)

    df.to_csv(f'{args.out}_feature.tsv', sep='\t', index=0)

    # predict & save result
    savedir = args.out + '_prediction.txt'
    with open(savedir, 'w') as f:
        f.write(str(predict(df))) # load model && predict

if __name__ == "__main__":
    # test_Ecoli_GCF_ASM584v2()
    main()