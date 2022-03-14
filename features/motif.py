import re
# import collections
# import pandas as pd
# import numpy as np
# import scipy as sp
# from IPython.core.display import display
# import seaborn as sns
# import warnings

__all__ = [
    "calc_hTelo",
    "calc_c_myc",
    "calc_bcl_2",
    "calc_Rb",
    "calc_RET",
    "calc_VEGF_A",
    "calc_c_ki_ras",
    "calc_c_kit",
    "calc_PDGF_A",
    "calc_c_myb",
    "calc_hTERT",
    "calc_HIF_1a",
    "calc_c_jun",
    "calc_ILPR",
    "calc_n_MYC",
]
def calc_hTelo(seq):
    rule = "CCCTAACCCTAACCCTAACCCT"
    counter = 0
    
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter

def calc_c_myc(seq):
    rule = "CCCCACCTTCCCCACCCTCCCCACCCTCCCC"
    counter = 0
    
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter

def calc_bcl_2(seq):
    rule = "CAGCCCCGCTCCCGCCCCCTTCCTCCCGCGCCCGCCCCT"
    counter = 0
    
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter

def calc_Rb(seq):
    rule = "GCCGCCCAAAACCCCCCG"
    counter = 0
    
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter

def calc_RET(seq):
    rule = "CCGCCCCCGCCCCGCCCCGCCCCTA"
    counter = 0
    
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter

def calc_VEGF_A(seq):
    rule = "GACCCCGCCCCCGGCCCGCCCCGG"
    counter = 0
    
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter

def calc_c_ki_ras(seq):
    rule = "GCTCCCTCCCTCCCTCCTTCCCTCCCTCCC"
    counter = 0
    
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter

def calc_c_kit(seq):
    rule = "CCCTCCTCCCAGCGCCCACCCT"
    counter = 0
    
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter

def calc_PDGF_A(seq):
    rule = "CCGCGCCCCTCCCCCGCCCCCGCCCCCGCCCCCCCCCCCCC"
    counter = 0
    
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter

def calc_c_myb(seq):
    rule = "TCCTCCTCCTCCTTCTCCTCCTCCTCCGTGTCCTCCTCCTCC"
    counter = 0
    
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter

def calc_hTERT(seq):
    rule = "CCCCGCCCCGTCCCGACCCCTCCCGGGTCCCCGGCCCAGCCCCCACCGGGCCCTCCCAGCCCCTCCCC"
    counter = 0
    
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter

def calc_HIF_1a(seq):
    rule = "CGCGCTCCCGCCCCCTCTCCCCTCCCCGCGCGCCCGAGCGCGCCTCCGCCCTTGCCCGCCCCCTG"
    counter = 0
    
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter

def calc_c_jun(seq):
    rule = "TAACCCCCTCCCCCTCCCCCCTTTAAT"
    counter = 0
    
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter

def calc_ILPR(seq):
    rule = "TGTCCCCACACCCCTGTCCCCACACCCCTGT}"
    counter = 0
    
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter

def calc_n_MYC(seq):
    rule = "ACCCCCTGCATCTGCATGCCCCCTCCCACCCCCT"
    counter = 0
    
    i_motifs = re.finditer(r"%s" % (rule), seq)
    motif = [[n.start(), n.end()] for n in i_motifs]
    if len(motif) > 0:
        for i in motif:
            counter += 1
    return counter