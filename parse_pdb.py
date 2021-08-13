#!/bin/env python

from Bio.PDB import *

class Pdb:
    def __init__(self, filename):
        parser = PDBParser()
        self.structure = parser.get_structure("struct", filename)
        self.models = list(self.structure.get_models())

# Obtain an object representing a structure (chain) within a model.
    def get_structure(self, s_id):
        struct = [s for s in self.structure.get_chains() if s.id == s_id][0]
        return struct

# Obtain an object representing residuum within a chain.
    def get_residue(self, res_id):
        r = [res for res in self.structure.get_residues() if res.id == res_id][0]
        return r

# Obtain an object representing an atom within a residue
    def get_atom(self, atom_id, residue):
        a = [atom_id for atom in residue if atom_id in atom.id][0]
        return a

# Obtain information about the stored structure (number of models, structures, residues, atoms).
    def get_information(self):
# number of models:
        models = [m for m in self.structure.get_models()]
        n_models = len(models)

# number of structures:
        structures = [s for s in self.structure.get_chains()]
        n_structures = len(structures)

# number of residues:
        residues = [r for r in self.structure.get_residues()]
        n_residues = len(residues)

# number of atoms:
        atoms = [a for a in self.structure.get_atoms()]
        n_atoms = len(atoms)

        info = f"Number of models: {n_models}, number of structures: {n_structures}, number of residues: {n_residues}, number of atoms: {n_atoms}"
        return info

# Comptue the width of the structure (maximum of distance of any two atoms).
    def count_width(self):
        max_width = 0
        atoms = self.structure.get_atoms()
        for a1 in atoms:
            for a2 in atoms:
                d = a1-a2
                if d > max_width:
                    max_width = d
        return max_width

# Obtain list of atoms being in given distance from given ligand (HETATM).        
    def get_close_atoms(self, ligand, distance):
        atoms = [a for a in self.structure.get_atoms() if round(a - ligand) == distance]
        return atoms

# Obtain list of residues being in given distance from given ligand (HETATM).
    def get_close_residues(self, ligand, distance):
        atoms = self.get_close_atoms(ligand, distance)
        residues = list({a.get_parent() for a in atoms})
        return residues
if __name__ == "__main__":
    p = Pdb('3VQJ.pdb')
    print(p.structure)
    print(p.models)
    print(p.get_information())
    print(p.count_width())
    r = [r for r in p.structure.get_atoms()][0]
    print(r)
    print(p.get_close_residues("<Atom C>",5))
