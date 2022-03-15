#!/home/zhangjq/anaconda3/envs/september/bin/python

'''
@File    :   calc_repeat1.py
@Author  :   renshuai
@Version :   1.0
@Contact :   rens@tib.cas.cn
@License :   (C)Copyright 2022-2023, renshuai
'''

import sys
import os
# import pandas as pd
sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
import numpy as np
from collections import defaultdict
import sys 
sys.path.append("..") 
# import utils.FastFinder
from utils.RepeatFinder import call_repeatFinder

class seq_assess_repeat(object):
    def __init__(self, sequence=None, begin=0):
        # defaults
        self.rc = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        self.sequence = sequence.upper()
        self.begin = begin
        self.size = len(sequence)
        if sequence:
            self.load(sequence)

    def load(self, sequence):
        # input validation
        try:
            if 'U' in sequence.seq or 'u' in sequence.seq:
                # print(sequence.seq
                sequence = sequence.back_transcribe()
            sequence1 = str(sequence.seq).upper()
            sequence2 = str(sequence.seq.reverse_complement()).upper()

        except:
            try:
                if 'U' in sequence or 'u' in sequence:
                    # print(sequence
                    sequence = sequence.back_transcribe()
                sequence1 = str(sequence).upper()
                sequence2 = str(sequence.reverse_complement()).upper()
            except:
                try:
                    sequence1 = str(sequence)
                    sequence1 = sequence1.upper()
                    # print(sequence1
                    bases = list(sequence1)
                    bases = reversed([self.rc.get(base, base) for base in bases])
                    sequence2 = ''.join(bases)
                    # print(sequence2
                except:
                    # print("The sequence provided was neither a Biopython record, Biopython Seq class, nor a valid string with your DNA sequence. Check what you are passing to seq_assess."
                    self.bad_seqs.append(
                        "The sequence provided was neither a Biopython record, Biopython Seq class, nor a valid string with your DNA sequence. Check what you are passing to seq_assess.")

        self.sequence = sequence1.upper()
        self.sequence1 = sequence1.upper()
        self.sequence2 = sequence2.upper()

    def run(self):
        self.find_repeats(RL=8, Rmax=11)


    def find_repeats(self, RL=8,
                     Rmax=11,
                     Rthreshold=[0.69, 0.40, 0.90, 0.60],
                     Rwindow=[70, 500, 60],
                     Rtandem=5):
        ### code for testing the following six rules:
        # A) Any sequence that contains more that 69 % of  bases that are part of repeats 8 bases or longer *   "X= 69 Y=8" 10
        # B) Any single repeated sequence contains more that 40 % of  bases that are part of repeats 8 bases or longer *    "X=40 Y=8"  6
        # C) Any sequence that contains more that 90 % of  bases that are part of repeats 8 bases or longer within a window of 70 length *  "X=90 Y=8 Z=70" 10
        # D) Any sequence that contains more that 60 % of  bases that are part of repeats 8 bases or longer within a window of 500 length for sequences longer than 1000 *  "X=60 Y=8 Z=500 L=1000" 10
        # E) a tandem repeat of length 5 or greater beginning or ending at either termini of the sequence   X=5 10
        # F) repeats more than 12 bp long risk impede synthesis

        ## outputs for Operon Calculator

        # metrics for tracking repeats for ML analysis
        self.tandem = 0  # number of tandem repeats


        self.high_local_density_1 = 0  # number of windows that exceed local density screen (0.90 of sequence is in some way repetitive for window of 70bp)
        self.high_local_density_2 = 0  # number of windows that exceed local density screen (0.60 of sequence is in some way repetitive for window of 500bp)
        self.high_specific_density = 0  # number of windows that exceed specific density screen (0.40 of sequence is repetitive to a given repeat)
        self.high_total_density = 0  # whether sequence exceeds repeat screen (0.69 of sequence is repetitive in some way)

        self.terminal_repeat = 0
        self.repeat9 = 0

        # ============================================
        # calc_repeats_10_15_20_25_40_;arge
        # ============================================
        # RP = FastFinder()
        # self.RepeatDict = RP.get_repeat_dict([self.sequence], RL, verbose=False)

        self.sequence = [self.sequence]
        self.RepeatDict = call_repeatFinder(self.sequence, k_low=9, k_high=100, verb=False)

        total = np.zeros(self.size, dtype=bool)
        total2 = np.zeros(self.size)

        repeat9_numerator = 0

        # ==================================================
        # RepeatDitc --> location_list --> 1 repeat_9 2 terminal_repeat 3 val=[]
        # ==================================================
        for i in self.RepeatDict.keys():
            val = np.zeros(self.size, dtype=bool)
            # ==================================================
            # calculate location_list
            # ==================================================
            for j in self.RepeatDict[i].keys():
                location_list = list(self.RepeatDict[i][j][0])
                # ==================================================
                # calculate repeat 9 mer
                # ==================================================
                for k in location_list:
                    if len(i) > 8:
                        for l in range(8, len(i)):
                            repeat9_numerator += 9
                    # ==================================================
                    # calculate terminal repeat
                    # ==================================================
                    if (self.begin + k) < 60 or (self.begin + k) >= self.size - 60:
                        self.terminal_repeat += 1

                    # ==================================================
                    # calculate val = []
                    # ==================================================
                    for v in range(k, k + len(i)):
                        val[v] = True
            # ============================================
            # calculate ind_fract using val
            # ============================================
            ind_fract = float(sum(val.astype(int))) / self.size
            # ==================================================
            # calculate high_specific_density using ind_fract
            # ==================================================
            if ind_fract > Rthreshold[1]:  # B)
                self.high_specific_density += 1
            # ==================================================
            # calculate total, total2 using val
            # ==================================================
            total += val
            total2 += val.astype(int)

        # ============================================
        # calculate total_fract using total, size
        # calculate high_total_density using total_fract
        # ============================================
        total_fract = float(sum(total.astype(int))) / self.size
        if total_fract > Rthreshold[0]:  # A)
            self.high_total_density += 1

        # ==================================================
        # calculate repeat9 using  repeat9_numerator, size
        # ==================================================
        self.repeat9 = repeat9_numerator / float(self.size)    # gen_9_metric???
        # ============================================
        # calc_local_density
        # ============================================
        for i in range(len(Rwindow) - 1):  # C) & D)
            for j in range(self.size - Rwindow[i]):
                peek = float(sum(total[j:j + Rwindow[i]].astype(int)))
                if peek / Rwindow[i] > Rthreshold[2 + i]:
                    # self.bad_seqs.append("%i bp window starting at %i is too repeat rich. Redesign to remove "%(Rwindow[i],j))
                    if i == 0:
                        self.high_local_density_1 += 1
                    if i == 1:
                        self.high_local_density_2 += 1

        # E)
        # ============================================
        # calculate tandem
        # ============================================
        # split the terminal part and call tandem_subroutine to calculate tandem
        # using tandem_subroutine()
        # using get_tandem_dict()
        # ============================================
        term1 = self.sequence1[:Rwindow[2]]
        self.tandem_subroutine(term1, Rtandem, "five")
        term2 = self.sequence1[len(self.sequence1) - Rwindow[2]:]
        self.tandem_subroutine(term2, Rtandem, "three")

    def tandem_subroutine(self, term, Rtandem, end):
        # subroutine that finds tandem repeats using regex

        run = self.get_tandem_dict(term, Rtandem, 5, 0)
        if len(run) > 0:
            for i in run.keys():
                for j in range(len(run[i])):
                    self.tandem += 1

    def is_tandem(self, seq, start, gap, rulen, max_x):
        x = 0
        for i in range(rulen):
            if seq[start + i] != seq[start + i + gap + rulen]:
                x += 1
            if x > max_x:
                return False
        return True

    def hdist(self, seq, start, gap, rulen, max_x):
        x = 0
        for i in range(rulen):
            # print(len(seq)
            # print(start+i
            # print(start+i+gap+rulen
            if seq[start + i] != seq[start + i + gap + rulen]:
                x += 1
            if x > max_x:
                return False
        self.tandem += 1
        return True

    def get_tandem_dict(self, seq, rulen, max_g, max_x):
        tr_chains = defaultdict(list)
        for i in range(len(seq) - rulen - rulen - 1):
            for g in range(min(len(seq) - i - rulen - max_g, max_g + 1)):
                if self.hdist(seq, i, g, rulen, max_x):
                    tr_chains[i].append(i + rulen + g)
        return tr_chains

if __name__ == "__main__":
    seq = 'AACGCAGTCAGGCACCGTGTATGAAATCTAACAATGCGCTCATCGTCATCCTCGGCACCGTCACCCTGGATGCTGTAGGCATAGGCTTGGTTATGCCGGTGCTCTTCAACGCGGGAATACTGTTAAGATTCATATATTACCTCTCAAATAATTAGCGGAAGGTTGCGATGAAACCGGCAAGGCTCTCTCAAACTGTCGTTGCGCCCGGATGTTGGGGTGAGTTGCCCTGGGGCAATTACTACCGTGAGGCGCTGGAACAGCAGCTAAATCCGTGGTTTGCAAAAATGTATGGTTTCCATTTGCTTAAAATCGGTAATTTAAGCGCAGAAATCAATTCCGAAGCGTGCGCGGTCTCCCATCAGGTGAATGTTTCATCGCAGGGGTCGCCGATGCAGGTTCTGGCCGATCCGCTACATCTTCCTTTTGCAGATAAATCCGTCGATGTTTGTCTGCTGGCGCATACTTTGCCGTGGTGTACCGACCCGCACCGTTTATTGCGGGAAGCCGACCGCGTATTGATTGATGACGGTTGGCTGGTCATTAGTGGATTTAACCCGCTGAGTTTGATGGGGTTACGTAAACTGGTACCCGTTTTACGTAAAACACCGCCCTATAATAGTCGGATGTTTACCCTTATGCGGCAACTGGACTGGCTGTCTTTACTCAATTTCGAAGTGCTACATTATAGCCGTTTTCATGTCTTACCCTGGAAAAAGCAGGGGGGGCGGCTTTTAAATACGCATATCCCGGCGCTGGGCTGTTTACAGCTTATTGTGGCCCGTAAGCGGACCATCCCGCTTACGCTTAATCCGCTGCGACATAATAAAAGTAAAACCCCTATCCGCCAGACCGTTGGCGCCACCCGGCAATATCGCAAACCGGATGGCTAAGCTTCCGCCTGGTAGCCACTATCTTCCTGCGTGGGATTCATTGCCGCGGCGCGCGCCAGTTCATCACAACGCTCGTTTTCGGGATGGCCTGCATGGCCTTTGACCCATACCCATTTGATCTGATGCTGACCTAACGCAGCATCGAGACGTTTCCAGAGATCGACATTTTTTACGGGTTTCTTTTCCGCTGTTTTCCAGCCGCGTTTCTTCCAGTTATGAATCCATTGGGTAATTCCTTGCCGCACATATTGGCTGTCGGTGCTCAATGTTACTTCGCAATGTTCTTTTAACGCTTCAAGCGCGACGATCGCCGCCATCAGTTCCATACGGTTATTGGTGGTCAGCGTGTAACCTTCACTAAACGTTTTTTCATGACCGCGATAGCGTAGGATAGCGCCATAACCACCAGGCCCTGGATTCCCCAGGCAAGAGCCATCGGTGAAAATTTCTACCTGTTTAAGCATCTCTGGTAGACTTCCTGTAATTGAAATCGATAACAAAACGCA'
    assess = seq_assess_repeat(seq, begin=0)
    assess.run()
    print("Processing....")
    print(assess.repeat9)                 # 2.25716104392
    print(assess.high_local_density_1)       # 474
    print(assess.high_local_density_2)       # 711
    print(assess.high_specific_density)      # 0
    print(assess.high_total_density)         # 0
    print(assess.tandem)                     # 2
    print(assess.terminal_repeat)            # 3

    # df = pd.read_csv('TIB_repeat2.csv').iloc[:, 0:8]
    # # print(df)
    # for index, row in df.iterrows():
    #     print(index)
    #     #     print(row)
    #     # print(row[0]) # 'AT..CG'
    #     assess = seq_assess(row[0], begin=0)
    #     assess.run()
    #     # result = calc_repeat(row[0])  # data = {'repeat_25': 0, 'freq_repeat': 4, 'repeat_40': 30, 'repeat_large': 29, 'longest_repeat': 85, 'repeat_20': 11, 'repeat_10': 17, 'repeat_15': 0}
    #     #     print(df.iloc[index,1])
    #     df.iloc[index, 1] = assess.repeat9
    #     df.iloc[index, 2] = assess.high_local_density_1
    #     df.iloc[index, 3] = assess.high_local_density_2
    #     df.iloc[index, 4] = assess.high_specific_density
    #     df.iloc[index, 5] = assess.high_total_density
    #     df.iloc[index, 6] = assess.tandem
    #     df.iloc[index, 7] = assess.terminal_repeat
    # #
    # #
    # df.to_csv("TIB_repeat_feature22.csv")