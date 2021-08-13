#!/bin/env python
from Bio import SeqIO

class Fasta:
    def __init__(self, filename):
        fasta_file = SeqIO.parse(filename, "fasta")
        # create a list of record
        self.records = [record for record in fasta_file]

# return sequence given index of the record in the fasta file
    def get_sequence(self, index):
        record = self.records[index]
        return record.seq

# return decsription given index of the record in the fasta file
    def get_description(self, index):
        record = self.records[index]
        return record.description

# return length given index of the record in the fasta file
    def get_length(self, index):
        record = self.records[index]
        return len(record.seq)

# return subsequence given index of the record in the fasta file
    def get_subsequence(self, index, start = 0, end = -1):
        record = self.records[index]
        sequence = record.seq
        if end == -1:
            return sequence[start:]
        else:
            return sequence[start:end]

if __name__ == "__main__":
    f = Fasta("protein2.fa")
    print(f.get_sequence(0))
    print(f.get_description(1))
    print(f.get_length(1))
    print(f.get_subsequence(1,1,3))
