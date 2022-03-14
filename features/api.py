#!/home/zhangjq/anaconda3/envs/september/bin/python
# -*-coding:utf-8 -*-
'''
@File    :   feat_calc.py
@Time    :   2022/02/21 17:46:14
@Author  :   zhangjq
@Version :   1.0
@Contact :   zhangjq@tib.cas.cn
@License :   (C)Copyright 2022-2023, zhangjq
@Desc    :   Enjoy your dinner
'''
import re
import pandas as pd
from pandarallel import pandarallel

from Bio import SeqIO

from calc_repeat1 import calc_repeat
from calc_repeat2 import seq_assess_repeat
from hairpins import seq_assess_hairpin
# from features.mix_gc import calc_delta_tm, calc_seq_gc_high, calc_seq_gc_low, calc_terminal_high, calc_terminal_low, calc_tm_high, calc_tm_low, calc_delta_gc, calc_g_quad_motifs, calc_i_motifs, calc_patternruns, calc_poly_run, count_gc
# from features.mix_gc import calc_patternruns, calc_g_quad_motifs,
from mix_gc import calc_tm_low, calc_seq_gc_low, calc_i_motifs, count_gc, calc_tm_high, calc_delta_gc, calc_terminal_low
# from features.Enzyme_cutting_site import calc_NlaIII,calc_ApaI,calc_MseI,calc_BfaI
from Enzyme_cutting_site import calc_EagI, calc_BbrPI, calc_AfaI, calc_BanIII, calc_BshTI, calc_BbeI, calc_AsuII, calc_ClaI, calc_AatI, calc_EgeI, calc_ApaLI, calc_ApaI, calc_FspI
# from CpG_island import *
# from motif import *

pandarallel.initialize(nb_workers=80)

def _feat_calculator(seq: str):
    """In the future, we will plan to optimize the API design of feature_calculator."""
    result = calc_repeat(seq)
    assess2 = seq_assess_hairpin(seq, begin=0)
    assess2.run()
    assess = seq_assess_repeat(seq, begin=0)
    assess.run()
    longest_repeat = result["longest_repeat"]

    return pd.Series(data={
            'richest_hairpin': float(assess2.richest_hairpin),
            'Tm_low': float(calc_tm_low(seq, 20)),
            'repeat_9_metric': float(assess.repeat9),
            'longest_repeat': float(longest_repeat),
            'FspI': float(calc_FspI(seq)),
            'high_total_density': float(assess.high_total_density),
            'GC_long_l': float(calc_seq_gc_low(seq, window_length=100)),
            'terminal_hairpins': float(assess2.terminal_hairpins),
            'palindromes': float(assess2.palindromes),
            'i_motifs': float(calc_i_motifs(seq)),
            'high_local_density_2': float(assess.high_local_density_2),
            'size': float(len(seq)),
            'ApaI': float(calc_ApaI(seq)),
            'repeat_20': float(result["repeat_20"]),
            'ApaLI': float(calc_ApaLI(seq)),
            'EgeI': float(calc_EgeI(seq)),
            'high_local_density_1': float(assess.high_local_density_1),
            'repeat_10': float(result["repeat_10"]),
            'ClaI': float(calc_ClaI(seq)),
            'high_specific_density': float(assess.high_specific_density),
            'repeat_40': float(result["repeat_40"]),
            'AatI': float(calc_AatI(seq)),
            'GC_short_l': float(calc_seq_gc_low(seq, window_length=20)),
            'AsuII': float(calc_AsuII(seq)),
            'BbeI': float(calc_BbeI(seq)),
            'BshTI': float(calc_BshTI(seq)),
            'BanIII': float(calc_BanIII(seq)),
            'AfaI': float(calc_AfaI(seq)),
            'GC_term_l': float(calc_terminal_low(seq, 50)),
            'total_GC': float(count_gc(seq)),
            'Tm_high': float(calc_tm_high(seq, 20)),
            'BbrPI': float(calc_BbrPI(seq)),
            'dGC': float(calc_delta_gc(seq, 20)),
            'EagI': float(calc_EagI(seq)),
            }
    )

def calc_feat(fasta: str, outdir=None):
    # seqs = [str(seq.seq) for seq in SeqIO.parse(fasta, 'fasta') if len(seq) < 3000]
    seqs = [str(seq.seq) for seq in SeqIO.parse(fasta, 'fasta')]
    assert len(seqs) > 0
    data = pd.DataFrame(list(pd.Series(seqs).parallel_map(lambda x: _feat_calculator(x))))
    if outdir is not None:
        # outdir = re.split('/.', fasta)[-2] # FIXME: BUG: ata
        data.to_csv(outdir, index=0)
        print(f"{outdir} saved.")

def write_csv_to_fna():
    dfdir = '../../CaseStudy/SSC_DC4318.csv'
    
    df = pd.read_csv(dfdir).loc[:, ['seq', 'label']]

    fasta = '../../data/dachang4318-with-label.fna'

    with open(fasta, 'a') as f:
        for index, row in df.iterrows():
            f.write(f"> {row['label']} \n")
            f.write(str(row['seq']))
            f.write('\n')
    print("OK")

def test_Ecoli():
    # fasta = '../../data/GCF_ASM584v2_sorted_1000.fna'
    fasta = '../../data/dachang4318-with-label.fna'
    outdir = 'dachang4318_sep_feat-with-3000bp.csv'
    calc_feat(fasta, outdir)
    
def test_10_seq():
    dfdir = './10sequence.xlsx'
    df = pd.read_excel(dfdir).loc[:, ['Sequence', 'label']]
    fastadir = '10-seq.fna'
    with open(fastadir, 'a') as f:
        for index, row in df.iterrows():
            f.write(f"> {row['label']} \n")
            f.write(str(row['Sequence']))
            f.write('\n')
    outdir = '10-seq-with-34-feature.csv'
    calc_feat(fastadir, outdir)


if __name__ == "__main__":
    test_10_seq()
    # test_Ecoli()