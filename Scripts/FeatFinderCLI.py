#!python
'''
FeatFinderCLI reads molecules as SMILES from the first column of a tab, comma or space
separated file and annotates the atoms of the molecules with their pharmacophore property.

Use 'FeatFinderCLI.py --help' for further information
'''
from rdkit.Chem import FeatFinderCLI
FeatFinderCLI.main()
