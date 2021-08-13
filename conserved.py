#!/bin/env python

from parse_msa import Msa

def position_conservation(msa):
    columns = len(msa.alignment[0])
    rows = len(msa.alignment[:,0])
    scores = []
    for i in range(columns):
        col = msa.alignment[:,i]
        aa = set(col)
        if '-' in aa:
            aa.remove('-')
        occurances = [(col.count(a),a)  for a in aa]
        occurances.sort(reverse=True)
        score = occurances, occurances[0][0]/rows
        scores.append((occurances[0][0]/rows,i))
    return scores   

def get_top_positions(msa, n_top=5):
    scores = position_conservation(msa)
    scores.sort(reverse=True)
    positions = [i for s,i in scores[:n_top]] 
    return positions
 
if __name__ == "__main__":
    msa = Msa("p53_mafft_clustal.txt")
    s = position_conservation(msa)
    print(get_top_positions(msa))
