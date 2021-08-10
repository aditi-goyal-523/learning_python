import numpy as np 
import seqlib
import random
import argparse
import sys


"""
parser = argparse.ArgumentParser(description="Virterbi gene finder")
parser.add_argument("--file", required=True, type=str, metavar='<path>')
arg = parser.parse_args()

myfasta = seqlib.read_fasta(arg.file)

seq=myfasta[1]
"""

seq = " GGGGGGAAAAAAAA"






# probabilities
intron_prob = {"A": 0.3, "T": 0.3, "C": 0.2, "G": 0.2}
exon_prob = {"G": 0.3, "C": 0.3, "A": 0.2, "T": 0.2}
transition_prob = {"EE": 0.6, "EI": 0.4, "II": 0.7, "IE": 0.3}

length = len(seq)

row1 = [None] * length
row2 = [None] * length
matrix = [row1, row2]
matrix[0][0] = ("N", 0.5)
matrix[1][0] = ("N", 0.5)


def exon_score(nt, m, i):
    # pass in transition prob? dont make it global
    L = exon_prob[nt] * m[0][i - 1][1] * transition_prob["EE"]
    D = exon_prob[nt] * m[1][i - 1][1] * transition_prob["EI"]
    if L > D:
        return ("L", L)
    elif D > L:
        return ("D", L)
    elif D == L:
        print("coin toss")
        if random.randint(0, 1) == 0:
            return ("L", L)
        elif random.randint(0, 1) == 1:
            return ("D", D)


def intron_score(nt, m, i):
    L = intron_prob[nt] * m[1][i - 1][1] * transition_prob["II"]
    D = intron_prob[nt] * m[0][i - 1][1] * transition_prob["IE"]
    if L > D:
        return ("L", L)
    elif D > L:
        return ("D", L)
    elif D == L:
        if random.randint(0, 1) == 0:
            return ("L", L)
        elif random.randint(0, 1) == 1:
            return ("D", D)


for i in range(1, length):
    nt = seq[i]
    matrix[0][i] = exon_score(nt, matrix, i)
    matrix[1][i] = intron_score(nt, matrix, i)

for i in matrix:
    print(i)

