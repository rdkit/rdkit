#!/usr/bin/env python3

import rdkit
from rdkit import Chem

in_smi = 'c1ccncc1'
mol = Chem.MolFromSmiles(in_smi)
out_smi = Chem.MolToSmiles(mol)
assert(in_smi == out_smi)
print(rdkit.__version__)
