#
#  Copyright (C) 2026 The RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""Writes reference_descriptors.tsv: raw rdkit.Chem.Descriptors values for a
few molecules, used by the NormalizedDescriptors tests.

Usage: python make_reference.py > reference_descriptors.tsv
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

SMILES = [
  'CC(=O)Oc1ccccc1C(=O)O',
  'CN1CCC[C@H]1c1cccnc1',
  'CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O',
  'O=C(O)C[C@](O)(CC(=O)O)C(=O)O',
  'Cn1cnc2c1c(=O)n(C)c(=O)n2C',
  'C1CC2(C1)CCNCC2',
  'CCN(CC)C(=O)[C@H]1CN(C)[C@@H]2Cc3c[nH]c4cccc(c34)C2=C1',
  'c1ccc2c(c1)ccc1ccccc12',
  'C[N+](C)(C)CCOP(=O)([O-])OCC',
  'ClC(Cl)(Cl)Br',
]


def main():
  names = [n for n, _ in Descriptors._descList]
  print('\t'.join(['smiles'] + names))
  for smi in SMILES:
    mol = Chem.MolFromSmiles(smi)
    vals = []
    for _, fn in Descriptors._descList:
      try:
        vals.append(repr(float(fn(mol))))
      except Exception:
        vals.append('nan')
    print('\t'.join([smi] + vals))


if __name__ == '__main__':
  main()
