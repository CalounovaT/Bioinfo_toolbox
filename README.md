# Bioinfo_toolbox
## Prepare steps
### Create conda environment and install needed packages
create an environment
`conda create --name toolbox`

activate it
`conda activate toolbox`

add bioconda channel
`conda config --add channels bioconda`

install necessary packages
`conda install biopython msms`

### Clone this repository
`git clone "https://github.com/CalounovaT/Bioinfo_toolbox/"`

## Run the toolbox
run the main program
`./main.py`

### 1) Processing FASTA files
Obtain description of an entry given its index in a file (0 in the example below)

`./main.py fasta --description 0 proteins.fasta`

Obtain sequence of an entry given its index in a file (0 in the example below)

`./main.py fasta --sequence 0 proteins.fasta`

Return sequence length given its index in a file (0 in the example below)

`./main.py fasta --length 0 proteins.fasta`

Return subsequence of a sequence given its index in a file, subsequence start and end positions

`./main.py fasta --subsequence 0 1 6 proteins.fasta`

### 2) Measuring sequence similarity using Hamming distance
Measure Hamming distance of two sequences

`./main.py hamming 'ABC' 'ACC'`

### 3) Sequence alignment using edit distance
Measure the edit distance of two sequences

`./main.py edit 'ABC' 'AAC'`

Output their alignments

`./main.py edit --alignment 'ABC' 'AAC'`

### 4) Processing PDB files
Obtain information about the stored structure (number of models, structures, residues, atoms).

`./main.py pdb --information 1B0B.pdb`

Comptue the width of the structure (maximum of distance of any two atoms).

`./main.py pdb --width 1B0B.pdb`

Obtain list of atoms being in given distance from given ligand (HETATM). - give ligand and distance as arguments

`./main.py pdb --atoms ligand distance 1B0B.pdb`

Obtain list of residues being in given distance from given ligand (HETATM).

`./main.py pdb --residues ligand distance 1BOB.pdb`

### 5) Processing multiple sequence alignment
Retrieve sequence by its position

`./main.py msa --sequencepos 0 align.msa`

Retrieve sequence by its id

`./main.py msa --sequenceid "seq1" align.msa`

Retrieve given column from the MSA

`./main.py msa --column 0 align.msa`

Retrieve sum of pairs score of a column

`./main.py msa -spcolumn 0 align.msa`

Retrieve sum of pairs score of a msa

`./main.py msa -spmsa align.msa`

### 6) Conservation determination from multiple aligned sequences
compute conservation score of msa

`./main.py conservation --conservation align.msa`

 identify top N scoring positions in the msa
 
 `./main.py conservation --toppositions 5 align.msa`
 
 ### 7) Computing structure-related properties
Compute the diameter of the protein.

`./main.py properties --diameter 1B0B.pdb`

Compute the ratio of surface and buried amino acids.

`./main.py properties --ratio 1B0B.pdb`

Output data for a histogram of amino acids composition of buried and exposed amino acids.

` ./main.py properties --surfacecounts 1B0B.pdb`

`./main.py properties --buriedcounts 1B0B.pdb`

Quantify portion of polar amino acids in the core and on the surface of the protein.

`./main.py properties --surfacepolar 1B0B.pdb`

` ./main.py properties --buriedpolar 1B0B.pdb`
