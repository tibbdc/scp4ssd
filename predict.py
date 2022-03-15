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
#          python predict.py --fasta example.fna --out example_out.csv                        #
######################################################################################

def predict_fasta(data: pd.DataFrame, return_prob: bool = False):
    # datadir = os.path.join('./result', re.split('/.', fasta)[-2], 'csv')
    # if not os.path.exists(datadir):
    #     calc_feat(fasta)
    # data = pd.read_csv(datadir)
    # print("Save feature successfully... start predicting...")
    # assert len(data) > 
    clf = load('./september/model/september.joblib')

    prediction = np.expand_dims(clf.predict(data), axis=1)
    prob = clf.predict_proba(data)
    # print('-'*50)
    # print(prediction.shape)
    # print(prob.shape)
    # print('-'*50)
    if return_prob:    
        return np.concatenate((prediction, prob), axis=1)
    return prediction

def test_Ecoli_GCF_ASM584v2():
    # fasta = './data/GCF_000005845.2_ASM584v2_cds_from_genomic.fna'
    # fasta = './data/lt3000_4246.fna'
    # fasta = './data/test_3.fna'
    outdir = f'dachang4138_prediction.txt'
    # if not os.path.exists(outdir):
    #     os.makedirs(outdir)
    fasta = './examples/example.fna'

    datadir = './september/features/dachang4318_sep_feat-with-3000bp.csv'
    
    pd.DataFrame(predict_fasta()).to_csv(outdir)
    print("Successfully saved.")
    
    
    
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

    # predict & save result
    savedir = args.out
    with open(savedir, 'w') as f:
        f.write(str(predict_fasta(df))) # load model && predict

if __name__ == "__main__":
    test_Ecoli_GCF_ASM584v2()
    # main()