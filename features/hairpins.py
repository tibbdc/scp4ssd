#!/home/zhangjq/anaconda3/envs/september/bin/python
# coding=utf-8

import sys 
sys.path.append("..") 
from utils.PyVRNA import PyVRNA
from bisect import bisect
import collections
# from string import maketrans
from collections import OrderedDict


class seq_assess_hairpin(object):

    # 初始化函数，定义变量

    def __init__(self, sequence=None, begin=0, verbose=False):
        self.rc = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        self.comp_table = "".maketrans('ATGC', 'TACG')
        self.verbose = verbose
        self.begin = begin
        if sequence:
            self.load(sequence)

    # 传入序列，将序列预处理（包括小写变大写，得到反向互补序列）

    def load(self, sequence):
        try:
            if 'U' in sequence.seq or 'u' in sequence.seq:
                print(sequence.seq)
                sequence = sequence.back_transcribe()
            sequence1 = str(sequence.seq).upper()
            sequence2 = str(sequence.seq.reverse_complement()).upper()

        except:
            try:
                if 'U' in sequence or 'u' in sequence:
                    print(sequence)
                    sequence = sequence.back_transcribe()
                sequence1 = str(sequence).upper()
                sequence2 = str(sequence.reverse_complement()).upper()
            except:
                try:
                    sequence1 = str(sequence)
                    sequence1 = sequence1.upper()

                    bases = list(sequence1)
                    bases = reversed([self.rc.get(base, base) for base in bases])
                    sequence2 = ''.join(bases)

                except:
                    print("The sequence provided was neither a Biopython record, Biopython Seq class, nor a valid string with your DNA sequence. Check what you are passing to seq_assess.")
        self.size = len(sequence)
        self.sequence = sequence1.upper()
        self.sequence1 = sequence1.upper()
        self.sequence2 = sequence2.upper()

    # 计算gc含量:
    # 输入：字符串的seq序列，如：“AAATTTCCC”
    # 输出：浮点型的数字GC含量,如：0.12

    def count_gc(self, seq, window):
        counts = dict(collections.Counter(seq))
        if 'G' in counts.keys() and 'C' in counts.keys():
            val = float(counts['G'] + counts['C']) / window
        elif 'G' in counts.keys():
            val = float(counts['G']) / window
        elif 'C' in counts.keys():
            val = float(counts['C']) / window
        else:
            val = 0
        return val

    #
    def run(self):
        self.find_hairpins_windowed()
        self.h_metrics = OrderedDict([("long_hairpins", self.long_hairpins),
                                      ("strong_hairpins", self.strong_hairpins),
                                      ("wide_hairpins", self.wide_hairpins),
                                      ("longest_hairpin", self.longest_hairpin),
                                      ("richest_hairpin", self.richest_hairpin),
                                      ("palindromes", self.palindromes),
                                      ("terminal_hairpins", self.terminal_hairpins)])

    # 计算预测出来的hairpin的feature：

    # 传入lengths：list of length of subhairpin,
    # details：list of tuples (subsequence（子序列的字符串）,substructure（子序列的结构）,
    #                        start_pos of str（子序列的起始点），mfe（最小自由能数值）)
    # gquads:list of gquads（预测出来的一种特殊hairpins）

    # 函数运算过程中对于各个feature值进行累加！

    def hairpin_processing(self, index, lengths, details, gquads, l_threshold, gc_threshold, mfe, rc=False):
        for i, l in enumerate(lengths):
            if l > l_threshold and details[i][3] < mfe and l_threshold < 16:
                hpinx = details[i][0]  # 'AATCG'
                structure = details[i][1]  # (((...)))
                gc = self.count_gc(hpinx, len(hpinx))  # 0.333
                if rc:
                    location = self.size - (index + details[i][2] + len(hpinx))
                else:
                    location = index + details[i][2]
                location_tup = (self.begin + location, self.begin + location + len(hpinx))
                if l > self.longest_hairpin:
                    self.longest_hairpin = l

                if gc > gc_threshold:
                    if hpinx not in self.h_info.keys():
                        if any(i < 60 for i in location_tup) or any(i > (self.size - 60) for i in location_tup):
                            self.terminal_hairpins += 1

                        self.h_info[hpinx] = {
                            "locations": [(self.begin + location, self.begin + location + len(hpinx))],
                            "types": ["strong hairpin"]}  # Adding offset (begin) if sequence is a subsequence
                        self.strong_hairpins += 1

                        self.hairpins[hpinx] = (self.begin + location, structure)
                        if gc > self.richest_hairpin:
                            self.richest_hairpin = gc

                    else:
                        continue
            elif l > l_threshold and details[i][3] < mfe and l_threshold > 16 and l_threshold < 20:
                hpinx = details[i][0]
                structure = details[i][1]
                gc = self.count_gc(hpinx, len(hpinx))
                if rc:
                    location = self.size - (index + details[i][2] + len(hpinx))
                else:
                    location = index + details[i][2]
                location_tup = (self.begin + location, self.begin + location + len(hpinx))
                if hpinx not in self.h_info.keys():
                    if any(i < 60 for i in location_tup) or any(i > (self.size - 60) for i in location_tup):
                        self.terminal_hairpins += 1
                    self.h_info[hpinx] = {"locations": [(self.begin + location, self.begin + location + len(hpinx))],
                                          "types": [
                                              "long hairpin"]}  # Adding offset (begin) if sequence is a subsequence
                    self.wide_hairpins += 1
                    self.hairpins[hpinx] = (self.begin + location, structure)
                else:
                    continue
            elif l > l_threshold and details[i][3] < mfe and l_threshold > 20:
                hpinx = details[i][0]
                structure = details[i][1]
                gc = self.count_gc(hpinx, len(hpinx))
                if rc:
                    location = self.size - (index + details[i][2] + len(hpinx))
                else:
                    location = index + details[i][2]
                location_tup = (self.begin + location, self.begin + location + len(hpinx))
                if hpinx not in self.h_info.keys():
                    if any(i < 60 for i in location_tup) or any(i > (self.size - 60) for i in location_tup):
                        self.terminal_hairpins += 1
                    self.h_info[hpinx] = {"locations": [(self.begin + location, self.begin + location + len(hpinx))],
                                          "types": [
                                              "long hairpin"]}  # Adding offset (begin) if sequence is a subsequence
                    self.long_hairpins += 1
                    self.hairpins[hpinx] = (self.begin + location, structure)
            else:
                continue
        for gquad in gquads:
            hpinx = gquad[0]
            structure = gquad[1]
            if rc:
                location = self.size - (index + gquad[2] + len(hpinx))
            else:
                location = index + gquad[2]
            location_tup = (self.begin + location, self.begin + location + len(hpinx))
            if hpinx not in self.h_info.keys():
                if any(i < 60 for i in location_tup) or any(i > (self.size - 60) for i in location_tup):
                    self.terminal_hairpins += 1
                self.h_info[hpinx] = {"locations": [(self.begin + location, self.begin + location + len(hpinx))],
                                      "types": ["gquadruplex"]}  # Adding offset (begin) if sequence is a subsequence
                self.hairpins[hpinx] = (self.begin + location, structure)
        return

    #  计算出各个feature的最终结果
    # 最终结果分由两个部分组成：1. 通过 hairpin_processing 先计算出预测hairpin的各个feature值；
    #                           2. 通过fast_airpin  计算出特定的hairpin的各个feature值，并在上面结果基础上累加

    def find_hairpins_windowed(self, Hwindow=[50, 100, 50],
                               HL=[11, 17, 21],
                               HGC=0.8,
                               HMFE=[-10, -15, -20],
                               pal=11):
        self.h_info = {}
        self.hairpins = {}

        self.long_hairpins = 0
        self.strong_hairpins = 0
        self.wide_hairpins = 0
        self.longest_hairpin = 0
        self.richest_hairpin = 0
        self.palindromes = 0
        self.terminal_hairpins = 0
        self.g_quad_predictions = 0

        self.DNAmodel = PyVRNA(parameter_file="dna_mathews2004.par", dangles=0, noGU=True,
                               gquad=True)

        for j in range(self.size - 100):
            seq1 = self.sequence1[j:j + 100]
            seq2 = self.sequence2[self.size - (j + 100):self.size - j]
            for i, h in enumerate(Hwindow):
                lengths1, details1, gquads1 = self.get_hairpin_lengths(seq1, h)
                lengths2, details2, gquads2 = self.get_hairpin_lengths(seq2, h)
                self.hairpin_processing(j, lengths1, details1, gquads1, HL[i], HGC, HMFE[i])
                self.hairpin_processing(j, lengths2, details2, gquads2, HL[i], HGC, HMFE[i], rc=True)

        self.fast_hairpins(fast=False)

    # 计算预测的hairpins的信息
    # 传入字符串类型的序列，以及整型的maxspan（最大跨距）
    # 传出lengths：list of length of subhairpin,
    # details：list of tuples (subsequence,substructure,start_pos of str，mfe),
    # gquads:list of gquads

    def get_hairpin_lengths(self, seq, maxspan):
        # hairpin length subroutine
        # details = list of tuples (subsequence,substructure,start_pos of str)
        self.DNAmodel.settings.max_bp_span = maxspan
        fold = self.DNAmodel.RNAfold(seq)
        parsed_fold = self.DNAmodel.vienna2bp(fold.structure)
        bp_x = parsed_fold.bpx
        bp_y = parsed_fold.bpy
        bp_gquad = parsed_fold.gquad
        lengths = []
        details = []
        gquads = []
        while len(bp_x) > 0:
            indx = bisect(bp_x, bp_y[0])
            lengths.append(indx)  # save number of bp in each hairpin (height)
            details.append((seq[bp_x[0] - 1:bp_y[0]], fold.structure[bp_x[0] - 1:bp_y[0]], bp_x[0] - 1,
                            fold.energy))  # save string and location of each hairpin
            bp_x = bp_x[indx:]
            bp_y = bp_y[indx:]
        for i in range(len(bp_gquad) // 3):
            first = bp_gquad[3 * i][0]
            last = bp_gquad[3 * i][3]
            gquads.append((seq[first - 1:last], fold.structure[first - 1:last], first - 1))
        return lengths, details, gquads

    #  计算特定的hairpins的各个features的值

    # 传入 字符串序列，茎长度，环上碱基数，最大错配个数，GC阈值

    # 通过累加计算出特定hairpin的feature值

    def _is_hairpin_pass(self, seq, stem, loop, max_mismatch, gc_high):
        # quickmode Haripin subroutine
        # max_mismatch = 2  by default  最大错配
        # gc_high = [0, 0.8]
        i = 0
        while i < len(seq) - (stem + loop + stem) + 1:
            j, mismatch, count_gc = 0, 0, 0.0
            while j < stem:
                if seq[i + j] in ['G', 'C']:
                    count_gc += 1.0
                if seq[i + (stem + loop + stem) - j - 1] in ['G', 'C']:
                    count_gc += 1.0
                if seq[i + j] != seq[i + (stem + loop + stem) - j - 1].translate(self.comp_table):
                    mismatch += 1
                if mismatch > max_mismatch:
                    break
                j += 1
            if j == stem:
                gc_content = (count_gc / 2.0) / stem
                if gc_content >= gc_high:
                    hairpin = seq[i:i + (stem + loop + stem)]
                    if hairpin not in self.h_info.keys():
                        if i < 60 or i > len(seq) - 61:
                            self.h_info[hairpin] = {"locations": [(self.begin + i, i + (stem))],
                                                    "types": ["terminal hairpin"]}
                            self.terminal_hairpins += 1
                        if gc_high > 0.50:
                            self.h_info[hairpin] = {"locations": [(self.begin + i, i + (stem))],
                                                    "types": ["strong hairpin"]}
                            self.strong_hairpins += 1
                            if gc_content > self.richest_hairpin:
                                self.richest_hairpin = gc_content
                        else:
                            self.h_info[hairpin] = {"locations": [(self.begin + i, i + (stem))],
                                                    "types": ["long hairpin"]}
                            if loop > 0 and loop < 49:
                                self.long_hairpins += 1
                            elif loop == 0:
                                self.palindromes += 1
                            else:
                                self.wide_hairpins += 1
                            if gc_content > self.richest_hairpin:
                                self.richest_hairpin = gc_content
                    else:
                        self.h_info[hairpin]["locations"].append((self.begin + i, i + (stem)))
                        if i < 60 or i > len(seq) - 61:
                            self.h_info[hairpin]["types"].append("terminal hairpin")
                            self.terminal_hairpins += 1
                        if gc_high > 0.50:
                            self.h_info[hairpin]["types"].append("strong hairpin")
                            self.strong_hairpins += 1
                            if gc_content > self.richest_hairpin:
                                self.richest_hairpin = gc_content
                        else:
                            self.h_info[hairpin]["types"].append("long hairpin")
                            if loop > 0 and loop < 49:
                                self.long_hairpins += 1
                            elif loop == 0:
                                self.palindromes += 1
                            else:
                                self.wide_hairpins += 1
                            if gc_content > self.richest_hairpin:
                                self.richest_hairpin = gc_content
            i += 1

    # 通过调用_is_hairpin_pass 函数，计算出特定的hairpin 的各个feature值
    # 函数中定义了specific 的具体数值

    def fast_hairpins(self, fast=False):

        seq = self.sequence1
        stem, max_mismatch, gc_high = 11, 2, 0.80
        for loop in range(3, 48 + 1):
            if len(seq) >= (stem + loop + stem):  # 是否现在的 seq长度 大于 loop 长度 + 2*stem
                self._is_hairpin_pass(seq, stem, loop, max_mismatch, gc_high)
            else:
                break
        stem, max_mismatch, gc_high = 17, 3, 0.0
        for loop in range(3, 100 + 1):
            if len(seq) >= (stem + loop + stem):
                self._is_hairpin_pass(seq, stem, loop, max_mismatch, gc_high)
            else:
                break
        if len(seq) > 500:
            stem, max_mismatch, gc_high = 22, 5, 0.0
            for loop in range(100, 500 + 1):
                if len(seq) >= (stem + loop + stem):
                    self._is_hairpin_pass(seq, stem, loop, max_mismatch, gc_high)
                else:
                    break

        # Palindromes
        stem, loop, max_mismatch, gc_high = 11, 0, 1, 0.0
        if len(seq) >= (stem + loop + stem):
            self._is_hairpin_pass(seq, stem, loop, max_mismatch, gc_high)


if __name__ == "__main__":
    seq = 'TGCATGATCTACGTGCGTCACATGCACGCGTACCAACCAGAAAAGCCAACCTGCGGGTTGGCTTTTTTATGCAATCACTTCTCTGTTGGCACGAAAAGGGCAATAAGATTTACGGATTACTATCTTGACATGAAGTGTTAGACGTCATATAATCGTGGTGTTCCGCCAGCCTGCCAACCGCTTCAGATCCAGAAATGGAAAGTTGAAGTGAGGCAGGTCCGGTAGCAACTCGAAAGAGTGAGAAAAGAGGGGAGCGGGAAACCGCTCCCCTTTTTTCGTTTAGACTACTATTGTCGTTTATTATCGCAACAGAGGGAAGTTCACTGACCTATTGACATAGGCAAGCCAGTATAGTATAATCACATACGCAGGCGAGTGAGCATAATCTTACCGAACTAGGAATAGTAAGTGGTAAGAAGGCCTGACCGTAATAAGCCTGAAAAGGCGACCAAAAAGGGGGGATTTTATCTCCCCTTTAATTTTTCCTGAATAAGCACTGTTGATAATCGCAATCTGTCTCTTCGTGAAAAGTAGCTTGACAATCGCTGTCTACGTGAATATAATGAATTTTCAGCCCCCGGGTCGCCATCCATTTTGGCGTCGAAAGACGAAGTAAAATGAAGGCGAGACCGATATCAACTGGAAGCAGTGTCTGGTAGTCCTGGTAAGACGCGAACAGCGTCGCATCAGGCATATTGCCAACTAGAGACTCCCTGTATCGTTGAAAAGTTGGAACACTGTGAATCCTATTACTGATTGACATGACTCTCCAGCTGTGCTATAATTGTACTATCCATCGCAAGACAGATGTGTTTAAGAGCTAGAAATAGCACGTTTAAATAAGGCTAGTCCGTTTTCAACTTGAAAAAGTGATAACAAAGCCGGGTAATTCCCGGCTTTGTTGTATCGTGAACGACACTACTATTTCTTACGAGATACTTATTCTGGAAGCAACGGTAGAAATCATCCTTAGCGAAAGCTAAGGATTTTTTTTATCTGTTACACTGCGCCATATGCTCAGATTCAGTAGACCGCTGTTGTAGTAATGCAGACACTTGCGGTCCATCTCG'
    assess = seq_assess_hairpin(seq, begin=0)
    assess.run()
    print(assess.long_hairpins)  # 0
    print(assess.longest_hairpin)  # 15
    print(assess.strong_hairpins)  # 0
    print(assess.richest_hairpin)  # 0
    print(assess.terminal_hairpins)  # 0
    print(assess.wide_hairpins)  # 16
    print(assess.palindromes)  # 0