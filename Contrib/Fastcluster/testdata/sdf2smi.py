import argparse

from rdkit import Chem
from rdkit.Chem.rdmolfiles import SmilesWriter

parser = argparse.ArgumentParser()
parser.add_argument('inputfile', help="sdf filename for convert to smiles")
args = parser.parse_args()
sdf = Chem.SDMolSupplier(args.inputfile)
writer = SmilesWriter("converted.smi")

for mol in sdf:
  writer.write(mol)
writer.close()
