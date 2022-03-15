#!/home/zhangjq/anaconda3/envs/september/bin/python
'''
@File    :   calc_repeat1.py
@Author  :   renshuai
@Version :   1.0
@Contact :   rens@tib.cas.cn
@License :   (C)Copyright 2022-2023, renshuai
'''

import re

__all__ = [
    "calc_BfaI",
    "calc_MseI",
    "calc_ApaI",
    "calc_NlaIII",
    "calc_AatI",
    "calc_AatII",
    "calc_Acc16I",
    "calc_AccII",
    "calc_AccIII",
    "calc_AclI",
    "calc_AcvI",
    "calc_AfaI",
    "calc_AfeI",
    "calc_AflII",
    "calc_AgeI",
    "calc_AhlI",
    "calc_Alw441",
    "calc_AluI",
    "calc_Aor51HI",
    "calc_ApaLI",
    "calc_AscI",
    "calc_AseI",
    "calc_Asp718I",
    "calc_AsuII",
    "calc_AvaI",
    "calc_AviII",
    "calc_AvrII",
    "calc_BalI",
    "calc_BamHI",
    "calc_BanIII",
    "calc_BbeI",
    "calc_BbrPI",
    "calc_BbuI",
    "calc_BcuI",
    "calc_BclI",
    "calc_BfrI",
    "calc_BfrBI",
    "calc_BglII",
    "calc_BlnI",
    "calc_BseCI",
    "calc_BsePI",
    "calc_BseX3I",
    "calc_BshTI",
    "calc_Bsp1407I",
    "calc_Bsp19I",
    "calc_BspDI",
    "calc_BspEI",
    "calc_BsrGI",
    "calc_BssHII",
    "calc_BstUI",
    "calc_ClaI",
    "calc_DpnII",
    "calc_DraI",
    "calc_EagI",
    "calc_EcoRI",
    "calc_EcoRV",
    "calc_EgeI",
    "calc_FseI",
    "calc_FspI",
    "calc_HaeIII",
    "calc_HincII",
    "calc_HindIII",
    "calc_HinfI",
    "calc_HpaI",
    "calc_HpaII",
    "calc_KasI",
    "calc_KpnI",
    "calc_MboI",
    "calc_MfeI",
    "calc_MluI",
    "calc_MscI",
    "calc_MspI",
    "calc_NaeI",
    "calc_NarI",
    "calc_NcoI",
    "calc_NdeI",
    "calc_NdeII",
    "calc_NgoMIV",
    "calc_NheI",
    "calc_NotI",
    "calc_NruI",
    "calc_NsiI",
    "calc_PacI",
    "calc_PciI",
    "calc_PhoI",
    "calc_PmeI",
    "calc_PmlI",
    "calc_PsiI",
    "calc_PstI",
    "calc_PvuI",
    "calc_PvuII",
    "calc_RsaI",
    "calc_SacI",
    "calc_SacII",
    "calc_SalI",
    "calc_SbfI",
    "calc_ScaI",
    "calc_SfoI",
    "calc_SmaI",
    "calc_SnaBI",
    "calc_SpeI",
    "calc_SphI",
    "calc_SspI",
    "calc_SstI",
    "calc_SstII",
    "calc_StuI",
    "calc_SwaI",
    "calc_TaqI",
    "calc_TliI",
    "calc_VspI",
    "calc_XbaI",
    "calc_XhoI",
    "calc_XmaI",
    "calc_ALL_Enzyme_sites",
]

def calc_AatI(seq):
    rule = "AGGCCT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_AatII(seq):
    rule = "GACGTC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_Acc16I(seq):
    rule = "TGCGCA"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_AccII(seq):
    rule = "CGCG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_AccIII(seq):
    rule = "TCCGGA"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_AclI(seq):
    rule = "AACGTT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_AcvI(seq):
    rule = "CACGTG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_AfaI(seq):
    rule = "GTAC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_AfeI(seq):
    rule = "AGCGCT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_AflII(seq):
    rule = "CTTAAG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_AgeI(seq):
    rule = "ACCGGT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_AhlI(seq):
    rule = "ACTAGT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_Alw441(seq):
    rule = "GTGCAC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_AluI(seq):
    rule = "AGCT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_Aor51HI(seq):
    rule = "AGCGCT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_ApaI(seq):
    rule = "GGGCCC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_ApaLI(seq):
    rule = "GTGCAC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_AscI(seq):
    rule = "GGCGCGCC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_AseI(seq):
    rule = "ATTAAT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_Asp718I(seq):
    rule = "GGTACC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_AsuII(seq):
    rule = "TTCGAA"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_AvaI(seq):
    rule = "C[CT]CG[AG]G"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_AviII(seq):
    rule = "TGCGCA"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_AvrII(seq):
    rule = "CCTAGG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_BalI(seq):
    rule = "TGGCCA"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_BamHI(seq):
    rule = "GGATCC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_BanIII(seq):
    rule = "ATCGAT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_BbeI(seq):
    rule = "GGCGCC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_BbrPI(seq):
    rule = "CACGTG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_BbuI(seq):
    rule = "GCATGC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_BcuI(seq):
    rule = "ACTAGT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_BclI(seq):
    rule = "TGATCA"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_BfaI(seq):
    rule = "CTAG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_BfrI(seq):
    rule = "CTTAAG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_BfrBI(seq):
    rule = "ATGCAT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_BglII(seq):
    rule = "AGATCT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_BlnI(seq):
    rule = "CCTAGG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_BseCI(seq):
    rule = "ATCGAT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_BsePI(seq):
    rule = "GCGCGC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_BseX3I(seq):
    rule = "CGGCCG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_BshTI(seq):
    rule = "ACCGGT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_Bsp1407I(seq):
    rule = "TGTACA"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_Bsp19I(seq):
    rule = "CCATGG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_BspDI(seq):
    rule = "ATCGAT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_BspEI(seq):
    rule = "TCCGGA"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_BsrGI(seq):
    rule = "TGTACA"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_BssHII(seq):
    rule = "GCGCGC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_BstUI(seq):
    rule = "CGCG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_ClaI(seq):
    rule = "ATCGAT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_DpnII(seq):
    rule = "GATC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_DraI(seq):
    rule = "TTTAAA"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_EagI(seq):
    rule = "CGGCCG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_EcoRI(seq):
    rule = "GAATTC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_EcoRV(seq):
    rule = "GATATC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_EgeI(seq):
    rule = "GGCGCC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_FseI(seq):
    rule = "GGCCGGCC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_FspI(seq):
    rule = "TGCGCA"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_HaeIII(seq):
    rule = "GGCC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_HincII(seq):
    rule = "GT[CT][AG]AC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_HindIII(seq):
    rule = "AAGCTT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_HinfI(seq):
    rule = "GA[ATGC]TC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_HpaI(seq):
    rule = "GTTAAC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_HpaII(seq):
    rule = "CCGG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_KasI(seq):
    rule = "GGCGCC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_KpnI(seq):
    rule = "GGTACC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_MboI(seq):
    rule = "GATC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_MfeI(seq):
    rule = "CAATTG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_MluI(seq):
    rule = "ACGCGT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_MscI(seq):
    rule = "TGGCCA"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_MseI(seq):
    rule = "TTAA"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_MspI(seq):
    rule = "CCGG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_NaeI(seq):
    rule = "GCCGGC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_NarI(seq):
    rule = "GGCGCC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_NcoI(seq):
    rule = "CCATGG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_NdeI(seq):
    rule = "CATATG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_NdeII(seq):
    rule = "GATC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_NgoMIV(seq):
    rule = "GCCGGC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_NheI(seq):
    rule = "GCTAGC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_NlaIII(seq):
    rule = "CATG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_NotI(seq):
    rule = "GCGGCCGC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_NruI(seq):
    rule = "TCGCGA"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_NsiI(seq):
    rule = "ATGCAT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_PacI(seq):
    rule = "TTAATTAA"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_PciI(seq):
    rule = "ACATGT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_PhoI(seq):
    rule = "GGCC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_PmeI(seq):
    rule = "GTTTAAAC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_PmlI(seq):
    rule = "CACGTG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_PsiI(seq):
    rule = "TTATAA"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_PstI(seq):
    rule = "CTGCAG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_PvuI(seq):
    rule = "CGATCG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_PvuII(seq):
    rule = "CAGCTG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_RsaI(seq):
    rule = "GTAC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_SacI(seq):
    rule = "GAGCTC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_SacII(seq):
    rule = "CCGCGG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_SalI(seq):
    rule = "GTCGAC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_SbfI(seq):
    rule = "CCTGCAGG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_ScaI(seq):
    rule = "AGTACT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_SfoI(seq):
    rule = "GGCGCC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_SmaI(seq):
    rule = "CCCGGG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_SnaBI(seq):
    rule = "TACGTA"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_SpeI(seq):
    rule = "ACTAGT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_SphI(seq):
    rule = "GCATGC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_SspI(seq):
    rule = "AATATT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_SstI(seq):
    rule = "GAGCTC"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_SstII(seq):
    rule = "CCGCGG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_StuI(seq):
    rule = "AGGCCT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_SwaI(seq):
    rule = "ATTTAAAT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_TaqI(seq):
    rule = "TCGA"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_TliI(seq):
    rule = "CTCGAG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_VspI(seq):
    rule = "ATTAAT"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter
def calc_XbaI(seq):
    rule = "TCTAGA"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter

def calc_XhoI(seq):
    rule = "CTCGAG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter

def calc_XmaI(seq):
    rule = "CCCGGG"
    counter = 0   
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter

def calc_ALL_Enzyme_sites(seq):
    sum = 0
    sum += calc_BfaI(seq)
    sum += calc_MseI(seq)
    sum += calc_ApaI(seq)
    sum += calc_NlaIII(seq)
    sum += calc_AatI(seq)
    sum += calc_AatII(seq)
    sum += calc_Acc16I(seq)
    sum += calc_AccII(seq)
    sum += calc_AccIII(seq)
    sum += calc_AclI(seq)
    sum += calc_AcvI(seq)
    sum += calc_AfaI(seq)
    sum += calc_AfeI(seq)
    sum += calc_AflII(seq)
    sum += calc_AgeI(seq)
    sum += calc_AhlI(seq)
    sum += calc_Alw441(seq)
    sum += calc_AluI(seq)
    sum += calc_Aor51HI(seq)
    sum += calc_ApaLI(seq)
    sum += calc_AscI(seq)
    sum += calc_AseI(seq)
    sum += calc_Asp718I(seq)
    sum += calc_AsuII(seq)
    sum += calc_AvaI(seq)
    sum += calc_AviII(seq)
    sum += calc_AvrII(seq)
    sum += calc_BalI(seq)
    sum += calc_BamHI(seq)
    sum += calc_BanIII(seq)
    sum += calc_BbeI(seq)
    sum += calc_BbrPI(seq)
    sum += calc_BbuI(seq)
    sum += calc_BcuI(seq)
    sum += calc_BclI(seq)
    sum += calc_BfrI(seq)
    sum += calc_BfrBI(seq)
    sum += calc_BglII(seq)
    sum += calc_BlnI(seq)
    sum += calc_BseCI(seq)
    sum += calc_BsePI(seq)
    sum += calc_BseX3I(seq)
    sum += calc_BshTI(seq)
    sum += calc_Bsp1407I(seq)
    sum += calc_Bsp19I(seq)
    sum += calc_BspDI(seq)
    sum += calc_BspEI(seq)
    sum += calc_BsrGI(seq)
    sum += calc_BssHII(seq)
    sum += calc_BstUI(seq)
    sum += calc_ClaI(seq)
    sum += calc_DpnII(seq)
    sum += calc_DraI(seq)
    sum += calc_EagI(seq)
    sum += calc_EcoRI(seq)
    sum += calc_EcoRV(seq)
    sum += calc_EgeI(seq)
    sum += calc_FseI(seq)
    sum += calc_FspI(seq)
    sum += calc_HaeIII(seq)
    sum += calc_HincII(seq)
    sum += calc_HindIII(seq)
    sum += calc_HinfI(seq)
    sum += calc_HpaI(seq)
    sum += calc_HpaII(seq)
    sum += calc_KasI(seq)
    sum += calc_KpnI(seq)
    sum += calc_MboI(seq)
    sum += calc_MfeI(seq)
    sum += calc_MluI(seq)
    sum += calc_MscI(seq)
    sum += calc_MspI(seq)
    sum += calc_NaeI(seq)
    sum += calc_NarI(seq)
    sum += calc_NcoI(seq)
    sum += calc_NdeI(seq)
    sum += calc_NdeII(seq)
    sum += calc_NgoMIV(seq)
    sum += calc_NheI(seq)
    sum += calc_NotI(seq)
    sum += calc_NruI(seq)
    sum += calc_NsiI(seq)
    sum += calc_PacI(seq)
    sum += calc_PciI(seq)
    sum += calc_PhoI(seq)
    sum += calc_PmeI(seq)
    sum += calc_PmlI(seq)
    sum += calc_PsiI(seq)
    sum += calc_PstI(seq)
    sum += calc_PvuI(seq)
    sum += calc_PvuII(seq)
    sum += calc_RsaI(seq)
    sum += calc_SacI(seq)
    sum += calc_SacII(seq)
    sum += calc_SalI(seq)
    sum += calc_SbfI(seq)
    sum += calc_ScaI(seq)
    sum += calc_SfoI(seq)
    sum += calc_SmaI(seq)
    sum += calc_SnaBI(seq)
    sum += calc_SpeI(seq)
    sum += calc_SphI(seq)
    sum += calc_SspI(seq)
    sum += calc_SstI(seq)
    sum += calc_SstII(seq)
    sum += calc_StuI(seq)
    sum += calc_SwaI(seq)
    sum += calc_TaqI(seq)
    sum += calc_TliI(seq)
    sum += calc_VspI(seq)
    sum += calc_XbaI(seq)
    sum += calc_XhoI(seq)
    sum += calc_XmaI(seq)
    return sum