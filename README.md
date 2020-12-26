# StructCalc
Program to calculate protein backbone structure from solid state NMR spectra and Rosetta bioinformatic restraints

The original project was completed between spring 18' and spring 19', and resulted in the following publication

Lapin, J., Nevzorov, A.A. (2019) Validation of protein backbone structures calculated from NMR angular restraints using Rosetta. Journal of Biomolecular NMR, 73 (5), 229-244

Following publication, the constituent python and C codes were combined into an interactive software program. The code is easily implemented from a linux command line.

1. create_project

Bash script that gets the user started by creating all necessary directories to organize the workflow

2. geninp.py

Interactive python script that reads in raw NMR spectrum to be fit, and any other informative input files so that structure calculations can commence.
Generates multiple files which will be used as input for the structure fitting program

3. config.conf

Configuration script with options for structure fitting program

4. structcalc.c

C source code for structure fitting to produce numerous candidate structures by fitting peptide plane orientations to the solid state NMR resonances.
Needs to read in Ramachandran restraints from the files in the "RAMA" directory. Outputs list of starting orientations and phi/psi dihedral angles
so that the backbone can be reconstructed. Must be compiled to an executeable before running.

5. structanal.py

Interactive python script that allows the user to interface with PyRosetta. Allows the filtering of structure results using custom scoring functions
made up of the individual Rosetta scoring terms.
