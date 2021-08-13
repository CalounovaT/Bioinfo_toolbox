#!/bin/env python

from parse_fasta import Fasta
from Bio import pairwise2

def edit_distance(s1, s2):

    # globalms(s1, s2, identity = 0, mismatch = -1, open_gap = -1, extend_gap = -1 
    score = pairwise2.align.globalms(s1, s2, 0, -1, -1,- 1)[0].score
    return -int(score)

def get_alignments(s1, s2):
    alignments = pairwise2.align.globalms(s1, s2, 0, -1, -1, -1)
    alignments_pairs = [(a.seqA, a.seqB) for a in alignments]    
    return alignments_pairs

if __name__ == "__main__":
    #print(edit_distance('clock','lacks'))
    #print(get_alignments('clock', 'lacks'))
    #print(edit_distance('AGGCGTDA','ACDFVRFRFRG'))i
    r = get_alignments('cdfdjc', 'lacks')
    print(r)
