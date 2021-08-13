#!/bin/env python

from Bio import AlignIO
from Bio.Align import substitution_matrices

# Read and parse MSA.
class Msa:
    def __init__(self, filename):
        self.alignment = AlignIO.read(filename, "clustal")
# Retrieve sequence by its position or ID
    def get_sequence_by_position(self, seq_position):
        return self.alignment[seq_position].seq

    def get_sequence_by_id(self, seq_id):
        s = [s for s in self.alignment if s.id == seq_id][0]
        return s.seq

# Retrieve given column from the MSA
    def get_column(self, column_index):
        c = [s.seq[column_index] for s in self.alignment]
        return "".join(c)

# Retrieve sum of pairs score of a column and whole MSA with respect to given scoring matrix.
    def get_score_column(self, column, matrix = substitution_matrices.load("PAM250")):
        score = 0
        n_rows = len(self.get_column(0))
        for i in range(n_rows):
            for j in range(i,n_rows):
                a = self.alignment[i,column]
                b = self.alignment[j,column]
                if a == '-':
                    a = '*'
                if b == '-':
                    b = '*'
                score += matrix[a,b]
        return score           

    def get_score_msa(self, matrix = substitution_matrices.load("PAM250")):
        n_columns = len(self.alignment[0].seq)
        scores = [self.get_score_column(i, matrix) for i in range(n_columns)]
        return sum(scores)

if __name__ == "__main__":
    msa = Msa("p53_mafft_clustal.txt")
    print(msa.get_sequence_by_position(0))
    print(msa.get_sequence_by_id('UniRef90_A0A151'))
    print(msa.get_column(8))
    print(msa.get_score_column(1))
    print(msa.get_score_msa())
    #print(substitution_matrices.load("PAM250"))
