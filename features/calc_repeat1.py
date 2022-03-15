#!/home/zhangjq/anaconda3/envs/september/bin/python
# coding=utf-8
'''
@File    :   calc_repeat1.py
@Author  :   renshuai
@Version :   1.0
@Contact :   rens@tib.cas.cn
@License :   (C)Copyright 2022-2023, renshuai
'''

import sys 
sys.path.append("..") 
from utils.RepeatFinder import call_repeatFinder

# from bisect import bisect


def calc_repeat(seq):
    longest_repeat = 0
    freq_repeat = 0
    repeat_10 = 0
    repeat_15 = 0
    repeat_20 = 0
    repeat_25 = 0
    repeat_40 = 0
    repeat_large = 0

    seq = [seq]
    # print("Processing repeat")
    RepeatDict = call_repeatFinder(seq, k_low=9, k_high=100, verb=False)
    # print RepeatDict
    # RepeatDict = {'TGCCCTCGC': {'direct': {0: [2508, 2204]}}, 'GAACTGCGTC': {'direct': {0: [2131, 3277]}}, 'ATTCTGGAT': {'invert': {0: [702]}, 'direct': {0: [2944]}}, 'CCCTTTTTAT': {'direct': {0: [763, 1487]}}, 'GCTGGCGCAG': {'invert': {0: [3094]}, 'direct': {0: [3519]}}, 'AGCATGAAA': {'direct': {0: [2402, 2599]}}, 'TGCTGCAGA': {'direct': {0: [3134, 2855]}}, 'CGCCAGCAG': {'direct': {0: [3097, 158]}}, 'CTGCCCGCC': {'direct': {0: [74, 1987]}}, 'GGCGCTGGC': {'invert': {0: [1999]}, 'direct': {0: [2754, 3639]}}, 'AATGGAAAT': {'invert': {0: [504]}, 'direct': {0: [2889]}}, 'ATTTTTCTG': {'direct': {0: [1406, 510]}}, 'GTATTCCTG': {'direct': {0: [3425, 443]}}, 'TTCGATATT': {'invert': {0: [951]}, 'direct': {0: [1655]}}, 'ACACAGCAGC': {'direct': {0: [1891, 350]}}, 'CGGGCAGGA': {'invert': {0: [844]}, 'direct': {0: [2310]}}, 'GGCGCTGGCG': {'direct': {0: [2754, 3639]}}, 'CCGCGCTGGC': {'direct': {0: [3632, 2327]}}, 'TCAGCGCGG': {'direct': {0: [454, 3383]}}, 'GCGCGGGCGGC': {'invert': {0: [2535]}, 'direct': {0: [2705]}}, 'AACCTGCTG': {'direct': {0: [2740, 404]}}, 'GCAGGAACTG': {'invert': {0: [3088]}, 'direct': {0: [3273]}}, 'GGCCGGAGCC': {'invert': {0: [2254]}, 'direct': {0: [2589]}}, 'CGCGCTGGC': {'direct': {0: [2328, 3633, 3516]}}, 'ACGGCGGCG': {'invert': {0: [61]}, 'direct': {0: [3659]}}, 'GGCGGGCAG': {'invert': {0: [74]}, 'direct': {0: [2308]}}, 'CCGGCGCTG': {'direct': {0: [2752, 2830]}}, 'TCGCCGGAA': {'invert': {0: [618]}, 'direct': {0: [2224]}}, 'GCCGCAGCC': {'invert': {0: [183]}, 'direct': {0: [2658]}}, 'ATGAACTGAC': {'direct': {0: [2186, 3215]}}, 'TCGGGAAAA': {'invert': {0: [396]}, 'direct': {0: [2480]}}, 'GCTGATGGC': {'invert': {0: [3019]}, 'direct': {0: [3186]}}, 'CGGTACAGC': {'direct': {0: [2969, 261]}}, 'GCAGGAAAA': {'invert': {0: [841]}, 'direct': {0: [3144]}}, 'GCGAATTTT': {'invert': {0: [2086]}, 'direct': {0: [3195]}}, 'CCTCACCGGT': {'direct': {0: [680, 127]}}, 'TGGCGCTGG': {'invert': {0: [2000]}, 'direct': {0: [2564, 3638]}}, 'GAGGAAACG': {'direct': {0: [1043, 1415]}}, 'CGCGCTGGCGC': {'direct': {0: [3633, 3516]}}, 'CAGGCGCGC': {'invert': {0: [725]}, 'direct': {0: [2362]}}, 'CGGGAGCTGG': {'invert': {0: [253]}, 'direct': {0: [3247]}}, 'ATTCCCGCC': {'invert': {0: [28]}, 'direct': {0: [2992]}}, 'GCAGTATCT': {'direct': {0: [584, 2150]}}, 'TGGCGCTGGC': {'invert': {0: [1999]}, 'direct': {0: [3638]}}, 'CGCCTGCCA': {'direct': {0: [728, 3074]}}, 'GGATTCAGATG': {'invert': {0: [419]}, 'direct': {0: [1020]}}, 'AGCGGCGCG': {'invert': {0: [2160]}, 'direct': {0: [3173]}}, 'GCGATTCCC': {'invert': {0: [31]}, 'direct': {0: [2296]}}, 'GAAAGCGCT': {'invert': {0: [381]}, 'direct': {0: [1959]}}, 'GCCAGCAGC': {'invert': {0: [2106]}, 'direct': {0: [3098]}}, 'ATACTGCCC': {'invert': {0: [582]}, 'direct': {0: [2200]}}, 'CGCGGCGGT': {'invert': {0: [3037]}, 'direct': {0: [3387]}}, 'GCGCTCTGC': {'invert': {0: [377]}, 'direct': {0: [2579]}}, 'AGTGCCGCC': {'direct': {0: [153, 2532]}}, 'GCGCTGGCGC': {'direct': {0: [3634, 2755, 3517]}}, 'GGCGGCGCAT': {'invert': {0: [58]}, 'direct': {0: [2710]}}, 'CTGGCGCTGG': {'direct': {0: [2563, 3637]}}, 'AAGCTGGTGTAG': {'direct': {0: [314, 885]}}, 'GCGGCGAGC': {'invert': {0: [2654]}, 'direct': {0: [3662]}}, 'CCGGAAGCG': {'invert': {0: [2174]}, 'direct': {0: [2227]}}, 'TTTATGAATT': {'invert': {0: [1437]}, 'direct': {0: [1721]}}, 'GCTGGAACTG': {'direct': {0: [3468, 3252]}}, 'CGCTGGCAG': {'invert': {0: [3077]}, 'direct': {0: [3238]}}, 'TTTCTCACC': {'direct': {0: [2018, 2978]}}, 'TGGCCGGAG': {'invert': {0: [340]}, 'direct': {0: [2588]}}, 'TATGTCCCT': {'direct': {0: [1306, 1386]}}, 'CACGCTTCC': {'invert': {0: [474]}, 'direct': {0: [2172]}}, 'GAGGCTATC': {'direct': {0: [3649, 1669]}}, 'CTACATTCC': {'invert': {0: [1273]}, 'direct': {0: [2988]}}, 'GGCGCTGAC': {'direct': {0: [2832, 3054]}}, 'CACCATATC': {'direct': {0: [2023, 735]}}, 'ATCGAAACA': {'invert': {0: [1652]}, 'direct': {0: [2215]}}, 'CGGCGGGCAG': {'invert': {0: [1987]}, 'direct': {0: [2307]}}, 'CCTGTTTCA': {'invert': {0: [360]}, 'direct': {0: [1762]}}, 'TACCGCCGC': {'invert': {0: [2965]}, 'direct': {0: [3036]}}, 'TGATGAACTG': {'invert': {0: [1182]}, 'direct': {0: [2184]}}, 'AACCGGAAG': {'invert': {0: [2176]}, 'direct': {0: [2870]}}, 'GCCAGCCAG': {'direct': {0: [3457, 661]}}, 'CGACCCGGT': {'invert': {0: [200]}, 'direct': {0: [2382]}}, 'GGGGATTCAG': {'invert': {0: [422]}, 'direct': {0: [3156]}}, 'TTTTAATTC': {'direct': {0: [1433, 1555]}}, 'GGCTACTGG': {'direct': {0: [2633, 2558]}}, 'AATACGATG': {'direct': {0: [1297, 555]}}, 'TCCGCGCTG': {'invert': {0: [455]}, 'direct': {0: [2326]}}, 'GACTGGGCGA': {'invert': {0: [46]}, 'direct': {0: [3560]}}, 'GAACCAGCT': {'direct': {0: [250, 1034]}}, 'CGGTCGCTG': {'direct': {0: [1810, 3127]}}, 'CGAAACATCG': {'direct': {0: [2217, 1050]}}, 'GCGCTGGCG': {'direct': {0: [3640, 3634, 2755, 3517]}}, 'CACTGGCTG': {'invert': {0: [663]}, 'direct': {0: [1873]}}, 'CAGGTGCTG': {'direct': {0: [2851, 3463]}}, 'TGGATACGCTGGA': {'invert': {0: [694]}, 'direct': {0: [2912]}}, 'ATTTTTATG': {'direct': {0: [1728, 1718]}}, 'CTGCTGGCC': {'direct': {0: [2584, 2107]}}, 'TTCGCAGCG': {'direct': {0: [242, 2090]}}, 'TGACGGCACAG': {'direct': {0: [1169, 2354]}}, 'TGCTGGCCGG': {'direct': {0: [2585, 3596]}}}

    for i in RepeatDict.keys():
        # print(i) # repeat parts. e.g. 'AGTCCGTTA'

        # ===================================
        # 计算 repeat_10, repeat_15, repeat_20, repeat_25, repeat_40, repeat_large
        # ===================================
        #     #     lambda x,y:
        #     thres_list=[7,10,15,20,25]
        #     repeat_list
        #     lambda : map(repeatlist:> funcution lis>)

        if 7 < len(i) <= 10:
            for j in RepeatDict[i].keys():
                repeat_10 += len(list(RepeatDict[i][j][0]))
        elif 10 < len(i) <= 15:
            for j in RepeatDict[i].keys():
                repeat_15 += len(list(RepeatDict[i][j][0]))
        elif 15 < len(i) <= 20:
            for j in RepeatDict[i].keys():
                repeat_20 += len(list(RepeatDict[i][j][0]))
        elif 20 < len(i) <= 25:
            for j in RepeatDict[i].keys():
                repeat_25 += len(list(RepeatDict[i][j][0]))
        elif 25 < len(i) <= 40:
            for j in RepeatDict[i].keys():
                repeat_40 += len(list(RepeatDict[i][j][0]))
        elif len(i) > 40:
            for j in RepeatDict[i].keys():
                repeat_large += len(list(RepeatDict[i][j][0]))
        # ===================================
        # 计算 longest_repeat
        # ===================================
        if len(i) > longest_repeat and len(i):
            longest_repeat = len(i)

        # ===================================
        # 计算 freq_repeat
        # ===================================
        temp_freq_repeat = 0
        for j in RepeatDict[i].keys():
            # temp_freq_repeat += len(RepeatDict[i][j].values())
            temp_freq_repeat += len(RepeatDict[i][j][0])
            # print "{} : {} -------------- {}".format(i, RepeatDict[i][j], len(RepeatDict[i][j][0]))
            if temp_freq_repeat > freq_repeat:
                freq_repeat = temp_freq_repeat

            # if len(RepeatDict[i][j].values()) > freq_repeat:
            #     freq_repeat = len(RepeatDict[i][j].values()[0])

    result = {"longest_repeat": longest_repeat,
              "repeat_10": repeat_10,
              "repeat_15": repeat_15,
              "repeat_20": repeat_20,
              "repeat_25": repeat_25,
              "repeat_40": repeat_40,
              "repeat_large": repeat_large,
              "freq_repeat": freq_repeat
              }
    return result


# PATH = '/home/liaoxp/Notes/ren_s/SSC_FEATURE/Baci_repeat1.csv'
# df = pd.read_csv('Baci_repeat1.csv').iloc[:,0:9]
# # print(df)
# for index, row in df.iterrows():
#     print(index)
#     #     print(row)
#     # print(row[0]) # 'AT..CG'
#     result = calc_repeat(row[0])  # data = {'repeat_25': 0, 'freq_repeat': 4, 'repeat_40': 30, 'repeat_large': 29, 'longest_repeat': 85, 'repeat_20': 11, 'repeat_10': 17, 'repeat_15': 0}
#     #     print(df.iloc[index,1])
#     df.iloc[index, 1] = result["longest_repeat"]
#     df.iloc[index, 2] = result["repeat_10"]
#     df.iloc[index, 3] = result["repeat_15"]
#     df.iloc[index, 4] = result["repeat_20"]
#     df.iloc[index, 5] = result["repeat_25"]
#     df.iloc[index, 6] = result["repeat_40"]
#     df.iloc[index, 7] = result["repeat_large"]
#     df.iloc[index, 8] = result["freq_repeat"]
# #
# df.to_csv("new_Baci_repeat_feature1.csv")

# seq1 = 'CCACCCGGCAATATCGCAAACCGGATGGTTAACACAGATATATGTTTATGTTGTAATAACTAGGTAAGCTTAAGATAAGGAGGAAAGACATATGGAGAAAAAAATCACTGGATATACCACCGTTGATATATCCCAATGGCATCGTAAAGAACATTTTGAGGCATTTCAGTCAGTTGCTCAATGTACCTATAACCAGACCGTTCAGCTGGATATTACGGCCTTTTTAAAGACCGTAAAGAAAAATAAGCACAAGTTTTATCCGGCCTTTATTCACATTCTTGCCCGCCTGATGAATGCTCATCCGGAATTCCGTATGGCAATGAAAGACGGTGAGCTGGTGATATGGGATAGTGTTCACCCTTGTTACACCGTTTTCCATGAGCAAACTGAAACGTTTTCATCGCTCTGGAGTGAATACCACGACGATTTCCGGCAGTTTCTACACATATATTCGCAAGATGTGGCGTGTTACGGTGAAAACCTGGCCTATTTCCCTAAAGGGTTTATTGAGAATATGTTTTTCGTCTCAGCCAATCCCTGGGTGAGTTTCACCAGTTTTGATTTAAACGTGGCCAATATGGACAACTTCTTCGCCCCCGTTTTCACCATGGGCAAATATTATACGCAAGGCGACAAGGTGCTGATGCCGCTGGCGATTCAGGTTCATCATGCCGTTTGTGATGGCTTCCATGTCGGCAGAATGCTTAATGAATTACAACAGTACTGCGATGAGTGGCAGGGCGGGGCGTAATGACTAAGCTTCCGCCTGGTAGCCACTATCTTCCTGCGTGGGATT'
# print(calc_repeat(seq1))
# seq2 = 'ATGAAATTCACGATTCAAAAAGATCGTCTTGTTGAAAGTGTCCAAGATGTATTAAAAGCAGTTTCATCCAGAACCACGATTCCCATTCTGACTGGTATTAAAATTGTTGCATCAGATGATGGAGTATCCTTTACAGGGAGTGACTCAGATATTTCTATTGAATCCTTCATTCCAAAAGAAGAAGGAGATAAAGAAATCGTCACTATTGAACAGCCCGGAAGCATCGTTTTACAGGCTCGCTTTTTTAGTGAAATTGTAAAAAAATTGCCGATGGCAACTGTAGAAATTGAAGTCCAAAATCAGTATTTGACGATTATCCGTTCTGGTAAAGCTGAATTTAATCTAAACGGACTGGATGCTGATGAATATCCGCACTTGCCGCAGATTGAAGAGCATCATGCGATTCAGATCCCAACTGATTTGTTAAAAAATCTAATCAGACAAACAGTATTTGCAGTGTCCACCTCAGAAACACGCCCTATCTTGACAGGTGTAAACTGGAAAGTGGAGCAAAGTGAATTATTATGCACTGCAACGGATAGCCACCGTCTTGCATTAAGAAAGGCGAAACTTGATATTCCAGAAGACAGATCTTATAACGTCGTGATTCCGGGAAAAAGTTTAACTGAACTCAGCAAGATTTTAGATGACAACCAGGAACTTGTAGATATCGTCATCACAGAAACCCAAGTTCTGTTTAAAGCGAAAAACGTCTTGTTCTTCTCACGGCTTCTGGACGGGAATTATCCAGACACAACCAGCCTGATTCCGCAAGACAGCAAAACAGAAATCATTGTGAACACAAAAGAATTCCTTCAGGCCATTGATCGTGCATCTCTTTTAGCTAGAGAGGGACGCAACAACGTTGTAAAACTGTCCGCAAAACCGGCTGAATCCATTGAAATTTCTTCCAATTCGCCAGAAATCGGTAAAGTTGTGGAAGCAATTGTTGCGGATCAAATTGAAGGTGAGGAATTAAATATCTCTTTTAGTCCAAAATATATGCTGGATGCACTAAAGGTGCTTGAAGGAGCAGAAATACGCGTAAGCTTTACAGGCGCAATGAGACCTTTCTTAATTCGCACGCCGAATGATGAAACGATTGTACAGCTTATCCTTCCTGTCAGAACCTATTAA'
# print(calc_repeat(seq2))
