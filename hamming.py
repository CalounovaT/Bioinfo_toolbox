#!/bin/env python
from parse_fasta import Fasta

def hamming_dist(s1, s2):
    if len(s1) != len(s2):
        raise Exception(f"Sequences have different lengths! s1: {len(s1)}, s2: {len(s2)}")
    else:
        diff_chars = [0 if c1 == c2 else 1 for c1, c2 in list(zip(s1,s2))]
        return(sum(diff_chars))

if __name__ == "__main__":
    print(hamming_dist('ACT','ACC'))
    f = Fasta("protein2.fa")
    s1 = f.get_sequence(1)
    s2 = f.get_sequence(2)
    print(s1, s2)
    print(hamming_dist(s1,s2))
