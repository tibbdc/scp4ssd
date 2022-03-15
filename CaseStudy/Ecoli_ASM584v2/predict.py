
import pandas as pd

from joblib import load
from Bio import SeqIO

import sys 
sys.path.append("") 

from scp4ssd.features.calc_repeat1 import *
from scp4ssd.features.calc_repeat2 import *
from scp4ssd.features.hairpins import *
# from features.mix_gc import calc_delta_tm, calc_seq_gc_high, calc_seq_gc_low, calc_terminal_high, calc_terminal_low, calc_tm_high, calc_tm_low, calc_delta_gc, calc_g_quad_motifs, calc_i_motifs, calc_patternruns, calc_poly_run, count_gc
# from features.mix_gc import calc_patternruns, calc_g_quad_motifs,
from scp4ssd.features.mix_gc import *
# from features.Enzyme_cutting_site import calc_NlaIII,calc_ApaI,calc_MseI,calc_BfaI
from scp4ssd.features.Enzyme_cutting_site import *
from scp4ssd.features.CpG_island import *
from scp4ssd.features.motif import *

import argparse

######################################################################################
# python predict.py --fasta example.fna --out example_out.csv                        #
######################################################################################

def _feat_calculator(seq):
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

def predict_fasta(fasta, return_prob=False):
    seqs = [str(seq.seq) for seq in SeqIO.parse(fasta, 'fasta')]
    assert len(seqs) > 0
    data = pd.DataFrame(list(pd.Series(seqs).map(lambda x: _feat_calculator(x))))
    clf = load('./scp4ssd.joblib')
    if return_prob:
        return clf.predict(data), clf.predict_proba(data)
    return clf.predict_proba(data)

def test_Ecoli_GCF_ASM584v2():
    fasta = './GCF_000005845.2_ASM584v2_cds_from_genomic.fna'
    print(predict_fasta(fasta))
    


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
    # load model && predict
    predict_fasta(args.fasta)
    # save result
    savedir = args.out

if __name__ == "__main__":
    test_Ecoli_GCF_ASM584v2()
    # main()