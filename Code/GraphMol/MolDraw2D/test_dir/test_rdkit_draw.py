#!/usr/bin/env python

import os
import sys

from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, Draw

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
