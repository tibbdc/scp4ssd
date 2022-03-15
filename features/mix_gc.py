#!/home/zhangjq/anaconda3/envs/september/bin/python
#-*- coding: UTF-8 -*-
'''
@File    :   mix_gc.py
@Time    :   2022/03/15 13:08:21
@Author  :   zhangjq
@Version :   1.0
@Contact :   zhangjq@tib.cas.cn
@License :   (C)Copyright 2022-2023, zhangjq
@Desc    :   Enjoy dinner :D
'''

from collections import Counter
from Bio.SeqUtils import MeltingTemp as TM
import re
# import collections


__all__ = [
    "calc_tm",
    "count_gc",
    "calc_tm_high",
    "calc_tm_low",
    "calc_seq_gc_high",
    "calc_seq_gc_low",
    "calc_terminal_high",
    "calc_terminal_low",
    "calc_delta_gc",
    "calc_delta_tm",
    "mers",
    "calc_poly_run",
    "calc_g_quad_motifs",
    "calc_i_motifs",
    "calc_patternruns",
]

def calc_tm(seq):
    return TM.Tm_NN(seq, dnac1=250.0, saltcorr=7)

def count_gc(seq):
    counts = dict(Counter(seq))  # 用于追踪value出现次数
    '''
    >>> obj = Counter('aabbccc')
    >>> print obj
    Counter({'c': 3, 'a': 2, 'b': 2})
    '''
    if 'G' in counts.keys() and 'C' in counts.keys():
        val = float(counts['G'] + counts['C']) / len(seq)
    elif 'G' in counts.keys():
        val = float(counts['G']) / len(seq)
    elif 'C' in counts.keys():
        val = float(counts['C']) / len(seq)
    else:
        val = 0
    return val

def calc_tm_high(seq, window_length):
    count_tm_high = 0
    for i in range(len(seq)-window_length):
        subseq = seq[i:i+window_length]
        if calc_tm(subseq) >= 70:
            count_tm_high += 1
    return count_tm_high

def calc_tm_low(seq, window_length):
    count_tm_low = 0
    for i in range(len(seq)-window_length):
        subseq = seq[i:i+window_length]
        if calc_tm(subseq) <= 40:
            count_tm_low += 1
    return count_tm_low

def calc_seq_gc_high(seq, window_length=20):
    count_high = 0
    for i in range(len(seq)-window_length):
        subseq = seq[i:i+window_length]
        if window_length == 20:
            if count_gc(subseq) >= 0.8:#0.8
                count_high += 1
        if window_length == 100:
            if count_gc(subseq) >= 0.7:
                count_high += 1
    return count_high

def calc_seq_gc_low(seq, window_length=20):
    count_low = 0
    for i in range(len(seq)-window_length):
        subseq = seq[i:i+window_length]
        if window_length == 20:
            if count_gc(subseq) <= 0.2:# 0.2
                count_low += 1
        if window_length == 100:
            if count_gc(subseq) <= 0.3:
                count_low += 1
    return  count_low


def calc_terminal_high(seq, window_length):
    count_terminal_high = 0
    for i in range(window_length - 20):  # 滑动窗口
        subseq_5 = seq[i:i + 20]  # 取两个端
        s5_gc = count_gc(subseq_5)
        subseq_3 = seq[len(seq) - window_length + i:len(seq) - window_length + i + 20]  # 取两个端
        s3_gc = count_gc(subseq_3)

        if s5_gc > 0.7:
            count_terminal_high += 1
        if s3_gc > 0.7:
            count_terminal_high += 1
    return count_terminal_high

def calc_terminal_low(seq, window_length):
    count_terminal_low = 0
    for i in range(window_length - 20):  # 滑动窗口
        subseq_5 = seq[i:i + 20]  # 取两个端
        s5_gc = count_gc(subseq_5)
        subseq_3 = seq[len(seq) - window_length + i:len(seq) - window_length + i + 20]  # 取两个端
        s3_gc = count_gc(subseq_3)

        if s5_gc < 0.3:
            count_terminal_low += 1
        if s3_gc < 0.3:
            count_terminal_low += 1
    return count_terminal_low

def calc_delta_gc(seq, window_length):
    count = 0
    threshold = 0.5
    size = len(seq)
    gc_list = []
    for i in range(size - window_length):
        subseq = seq[i:i + window_length]
        gc_list.append(count_gc(subseq))
    seq_checker = [None, None]  # eliminate repeat count
    for i in range(len(gc_list) - 80):
        window = gc_list[i:i + 80]
        max_val = max(window)
        min_val = min(window)
        if max_val - min_val > threshold:
            max_val_seq = seq[i + window.index(max_val):i + window.index(max_val) + 20]
            min_val_seq = seq[i + window.index(min_val):i + window.index(min_val) + 20]
            if max_val_seq != seq_checker[1] or min_val_seq != seq_checker[0]:  # if this subseq is a new subseq, count++
                seq_checker = [min_val_seq, max_val_seq]
                count += 1
    return count


def calc_delta_tm(seq, window_length):
    count = 0
    threshold = 30
    size = len(seq)
    tm_list = []
    for i in range(size - window_length):
        subseq = seq[i:i + window_length]
        tm_list.append(calc_tm(subseq))
    seq_checker = [None, None]  # eliminate repeat count
    for i in range(len(tm_list) - 80):  # silde window
        window = tm_list[i:i + 80]
        max_val = max(window)
        min_val = min(window)
        if max_val - min_val > threshold:
            max_val_seq = seq[i + window.index(max_val):i + window.index(max_val) + 20]
            min_val_seq = seq[i + window.index(min_val):i + window.index(min_val) + 20]
            if max_val_seq != seq_checker[1] or min_val_seq != seq_checker[0]:  # if this subseq is a new subseq, count++
                seq_checker = [min_val_seq, max_val_seq]
                count += 1
    return count

# def calc_size(seq):
#     return len(seq)

def mers(length):
    """Generates multimers for sorting through list of 10mers based on user
    specification. Multimers generated act as the keys for generating a
    hashtable to eliminate undesired sequence patterns from those 10mers not
    found in the genome.

    Usage: mers(N) = 4^(N) unique Nmers
    """
    seq_list = ['']
    counter = 0
    while counter < length:
        for seq in seq_list:
            if len(seq) == counter:
                for x in ['A', 'T', 'C', 'G']:
                    seq_list.append(seq + x)
        counter += 1
    last_N_Mers = 4 ** length
    return seq_list[len(seq_list) - last_N_Mers:]

def calc_poly_run(seq):
    counter = 0
    cons = {'A': 13, 'C': 9, 'G': 9, 'T': 13}
    for key in cons:
        rule = "%s{%i,1000}" % (key, cons[key])
        polyruns = re.finditer(r"%s" % (rule), seq)
#         count_of_polyrun = reg_ex(rule, sequence, count_of_polyrun, "poly N run")
        one_poly = [[n.start(), n.end()] for n in polyruns]
        if len(one_poly) > 0:
            for i in one_poly:
                counter += 1
    return counter


def calc_g_quad_motifs(seq):
    rule = "G{3,5}[ATGC]{2,7}" * 3
    counter = 0

    g_quad_motifs = re.finditer(r"%s" % (rule), seq)

    motif = [[n.start(), n.end()] for n in g_quad_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter


def calc_i_motifs(seq):
    rule = "C{2,5}[ATGC]{2,7}" * 3
    counter = 0

    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter


def calc_patternruns(seq):
    counter = 0
    mer_size = [3, 2]
    mer_number = [6, 10]

    #     print(list(range(len(mer_size)))

    for i in range(len(mer_size)):
        for j in mers(mer_size[i]):
            rule = j * mer_number[i]
            pattern_runs = re.finditer(r"%s" % (rule), seq)
            motif = [[n.start(), n.end()] for n in pattern_runs]
            if len(motif) > 0:
                for xxi in motif:
                    counter += 1
    return counter

# def my_func(seq):
#     calc_patternruns(seq)
#     calc_i_motifs(seq)
#     calc_g_quad_motifs(seq)
#     calc_poly_run(seq)
#     calc_size(seq)
#     calc_delta_tm(seq,20)
#     calc_delta_gc(seq, 20)
#     calc_terminal_high(seq, 50)
#     calc_terminal_low(seq,50)
#     calc_seq_gc_high(seq, window_length=20)
#     calc_seq_gc_high(seq, window_length=100)
#     calc_seq_gc_low(seq, window_length=20)
#     calc_seq_gc_low(seq, window_length=100)
#     calc_tm_high(seq, 20)
#     calc_tm_low(seq, 20)
#     count_gc(seq)

# if __name__ == "__main__":
#     seq = 'TGCATGATCTACGTGCGTCACATGCACGCGTACCAACCAGAAAAGCCAACCTGCGGGTTGGCTTTTTTATGCAATCACTTCTCTGTTGGCACGAAAAGGGCAATAAGATTTACGGATTACTATCTTGACATGAAGTGTTAGACGTCATATAATCGTGGTGTTCCGCCAGCCTGCCAACCGCTTCAGATCCAGAAATGGAAAGTTGAAGTGAGGCAGGTCCGGTAGCAACTCGAAAGAGTGAGAAAAGAGGGGAGCGGGAAACCGCTCCCCTTTTTTCGTTTAGACTACTATTGTCGTTTATTATCGCAACAGAGGGAAGTTCACTGACCTATTGACATAGGCAAGCCAGTATAGTATAATCACATACGCAGGCGAGTGAGCATAATCTTACCGAACTAGGAATAGTAAGTGGTAAGAAGGCCTGACCGTAATAAGCCTGAAAAGGCGACCAAAAAGGGGGGATTTTATCTCCCCTTTAATTTTTCCTGAATAAGCACTGTTGATAATCGCAATCTGTCTCTTCGTGAAAAGTAGCTTGACAATCGCTGTCTACGTGAATATAATGAATTTTCAGCCCCCGGGTCGCCATCCATTTTGGCGTCGAAAGACGAAGTAAAATGAAGGCGAGACCGATATCAACTGGAAGCAGTGTCTGGTAGTCCTGGTAAGACGCGAACAGCGTCGCATCAGGCATATTGCCAACTAGAGACTCCCTGTATCGTTGAAAAGTTGGAACACTGTGAATCCTATTACTGATTGACATGACTCTCCAGCTGTGCTATAATTGTACTATCCATCGCAAGACAGATGTGTTTAAGAGCTAGAAATAGCACGTTTAAATAAGGCTAGTCCGTTTTCAACTTGAAAAAGTGATAACAAAGCCGGGTAATTCCCGGCTTTGTTGTATCGTGAACGACACTACTATTTCTTACGAGATACTTATTCTGGAAGCAACGGTAGAAATCATCCTTAGCGAAAGCTAAGGATTTTTTTTATCTGTTACACTGCGCCATATGCTCAGATTCAGTAGACCGCTGTTGTAGTAATGCAGACACTTGCGGTCCATCTCG'
#     print my_func(seq)


