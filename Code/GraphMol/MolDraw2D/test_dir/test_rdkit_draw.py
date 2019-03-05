#!/usr/bin/env python

from rdkit import RDConfig
import os, sys
from rdkit import Chem
from rdkit.Chem import Draw, AllChem

suppl = Chem.SDMolSupplier()

if sys.argv[1].endswith('.sdf'):
  suppl = Chem.SDMolSupplier(sys.argv[1])
elif sys.argv[1].endswith('.smi'):
  suppl = Chem.SmilesMolSupplier(sys.argv[1])
else:
  print('Need a file ending in .sdf or .smi')
  exit(1)

for mol in suppl:
  print(mol.GetProp('_Name'))
  fn = mol.GetProp('_Name') + '.png'
  AllChem.Compute2DCoords(mol)
  Draw.MolToFile(mol, fn)
