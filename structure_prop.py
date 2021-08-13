#!/bin/env python

from parse_pdb import Pdb
from Bio.PDB import *
from Bio.PDB.ResidueDepth import get_surface, residue_depth


# Compute the diameter of the protein.
def get_diameter(pdb):
    return pdb.count_width()

# Compute the ratio of surface and buried amino acids.
def get_surface_aa(pdb):
    surface = get_surface(pdb.structure)
    surface_aa = [r for r in pdb.structure.get_residues() if residue_depth(r,surface) < 2 and r.get_resname() != 'HOH']
    return surface_aa

def get_buried_aa(pdb):
    surface = get_surface(pdb.structure)
    buried_aa = [r for r in pdb.structure.get_residues() if residue_depth(r,surface) > 2 and r.get_resname() != 'HOH']
    return buried_aa

def get_ratio(pdb):
    surface = len(get_surface_aa(pdb))
    buried = len(get_buried_aa(pdb))
    total = surface + buried
    return (surface/total, buried/total)

def get_buried_counts(pdb):
    buried_aa = get_buried_aa(pdb)
    buried_names = [aa.get_resname() for aa in buried_aa]
    aa_names = {aa.get_resname() for aa in buried_aa}
    counts = [(buried_names.count(aa),aa) for aa in aa_names]
    counts.sort(reverse=True)
    return counts

def get_surface_counts(pdb):
    surface_aa = get_surface_aa(pdb)
    surface_names = [aa.get_resname() for aa in surface_aa]
    aa_names = set(surface_names)
    counts = [(surface_names.count(aa),aa) for aa in aa_names]
    counts.sort(reverse=True)
    return counts

POLAR = ['ARG','ASN','ASP','CYS','GLN','GLU','HIS','LYS','SER','THR','TYR']

def count_polar_on_surface(pdb):
    surface_aa = get_surface_aa(pdb)
    surface_names = [aa.get_resname() for aa in surface_aa]
    polar_aa = [aa for aa in surface_names if aa in POLAR]
    return len(polar_aa)/len(surface_aa)

def count_polar_buried(pdb):
    buried_counts = get_buried_counts(pdb)
    total = 0
    num = 0
    for count, aa in buried_counts:
        total += count
        if aa in POLAR:
            num += count
    return num/total    
    

if __name__ == "__main__":
    p = Pdb('1B0B.pdb')
    #print(get_diameter(p))
    print('surface aa:',len(get_surface_aa(p)))
    print('buried aa:',len(get_buried_aa(p)))
    #print(len(list(p.structure.get_residues())))
    print('surface vs buried:',get_ratio(p))
    #print(get_buried_counts(p)) 
    print('polar on surface:',count_polar_on_surface(p))
    print('polar buried:',count_polar_buried(p))
