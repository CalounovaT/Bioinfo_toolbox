#!/bin/env python
import argparse
import sys
from edit import *
from parse_fasta import *
from parse_pdb import *
from structure_prop import *
from conserved import *
from hamming import *
from parse_msa import *

arguments = argparse.ArgumentParser()
arguments_sub = arguments.add_subparsers()

arguments_fasta = arguments_sub.add_parser('fasta', help='Parse a fasta file.')
arguments_fasta.set_defaults(action='fasta')
arguments_fasta.add_argument('file')
arguments_fasta.add_argument('--description')#, dest='description')#, type=int)
arguments_fasta.add_argument('--sequence')#, dest='sequence')#, type=int)
arguments_fasta.add_argument('--length')#, dest='length')#, type=int)
arguments_fasta.add_argument('--subsequence', nargs = 3)#, dest='sub')

arguments_ham = arguments_sub.add_parser('hamming', help = 'Count the hamming distance.')
arguments_ham.set_defaults(action="hamming")
arguments_ham.add_argument('seq1')
arguments_ham.add_argument('seq2')

arguments_ed = arguments_sub.add_parser('edit', help='Count the edit distance, optionally show alignments')
arguments_ed.set_defaults(action='edit')
arguments_ed.add_argument('seq1')
arguments_ed.add_argument('seq2')
arguments_ed.add_argument('--alignments', action='store_true')

arguments_pdb = arguments_sub.add_parser('pdb', help = 'Parse a PDB file.')
arguments_pdb.set_defaults(action='pdb')
arguments_pdb.add_argument('file')
arguments_pdb.add_argument('--information', action='store_true')
arguments_pdb.add_argument('--width', action='store_true')
arguments_pdb.add_argument('--atoms', nargs=2)
arguments_pdb.add_argument('--residues', nargs=2)

arguments_msa = arguments_sub.add_parser('msa', help = 'Parse a MSA (Clustal) file')
arguments_msa.set_defaults(action='msa')
arguments_msa.add_argument('file')
arguments_msa.add_argument('--sequenceid')
arguments_msa.add_argument('--sequencepos')
arguments_msa.add_argument('--column')
arguments_msa.add_argument('--spcolumn')
arguments_msa.add_argument('--spmsa', action='store_true')

arguments_cons = arguments_sub.add_parser('conservation', help = 'Get MSA conservation score for every position and top scoring positions')
arguments_cons.set_defaults(action="conservation")
arguments_cons.add_argument('file')
arguments_cons.add_argument('--conservation', action = 'store_true')
arguments_cons.add_argument('--toppositions')

arguments_struct = arguments_sub.add_parser('properties', help = 'Get PDB structure properties.')
arguments_struct.set_defaults(action='properties')
arguments_struct.add_argument('file')
arguments_struct.add_argument('--diameter', action ='store_true')
arguments_struct.add_argument('--surfaceaa',action = 'store_true')
arguments_struct.add_argument('--buriedaa', action ='store_true')
arguments_struct.add_argument('--ratio', action = 'store_true')
arguments_struct.add_argument('--surfacecounts', action = 'store_true')
arguments_struct.add_argument('--buriedcounts', action = 'store_true')
arguments_struct.add_argument('--surfacepolar', action = 'store_true')
arguments_struct.add_argument('--buriedpolar', action = 'store_true')

if len(sys.argv) < 2:
    arguments.print_help()
    sys.exit(0)
else:
    config = arguments.parse_args()
if config.action == 'fasta':
    f = config.file
    fasta_file = Fasta(f)

    if config.sequence:
        ind = int(config.sequence)
        print(fasta_file.get_sequence(ind))

    if config.description:
        ind = int(config.description)
        print(fasta_file.get_description(ind))    

    if config.length:
        ind = int(config.length)
        print(fasta_file.get_length(ind))

    if config.subsequence:
        index, start, end = config.subsequence
        print(fasta_file.get_subsequence(int(index), int(start), int(end)))

elif config.action == "hamming":
    s1 = config.seq1
    s2 = config.seq2
    print(hamming_dist(s1, s2))        

elif config.action == 'edit':
    s1 = config.seq1
    s2 = config.seq2
    print(edit_distance(s1,s2))
    if config.alignments:
        a = get_alignments(s1, s2) 
        for a1,a2 in a:
            print(a1)
            print(a2)
            print('#'* 100)

elif config.action == 'pdb':
    f = config.file
    pdb = Pdb(f)
    
    if config.information:
       print(pdb.get_information()) 
    if config.width:
        print(pdb.count_width())
    if config.atoms:
        ligand, distance = config.atoms
        print(pdb.get_close_atoms(ligand, distance))
    if config.residues:
        ligand, distance = config.residues
        print(pdb.get_close_residues(ligand, distance))

elif config.action == 'msa':
    f = config.file
    msa = Msa(f)
    
    if config.sequenceid:
        print(msa.get_sequence_by_id(config.sequenceid))
    if config.sequencepos:
        print(msa.get_sequence_by_position(int(config.sequencepos)))
    if config.column:
        print(msa.get_column(int(config.column)))
    if config.spcolumn:
        print(msa.get_score_column(int(config.spcolumn)))
    if config.spmsa:
        print(msa.get_score_msa())

elif config.action == "conservation":
    f = config.file
    msa = Msa(f)

    if config.conservation:
        scores = [score for score, pos in position_conservation(msa)]
        print(scores)
    if config.toppositions:
        print(get_top_positions(msa, int(config.toppositions)))

elif config.action == 'properties':
    f = config.file
    pdb = Pdb(f)
    
    if config.diameter:
        print(get_diameter(pdb))
    if config.surfaceaa:
        print(get_surface_aa(pdb))
    if config.buriedaa:
        print(get_buried_aa(pdb))
    if config.ratio:
        print(get_ratio(pdb))
    if config.surfacecounts:
        print(get_surface_counts(pdb))
    if config.buriedcounts:
        print(get_buried_counts(pdb))
    if config.surfacepolar:
        print(count_polar_on_surface(pdb))
    if config.buriedpolar:
        print(count_polar_buried(pdb))
else:
    pass
