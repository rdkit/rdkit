# Copyright (C) 2002-2026 greg Landrum and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Basic EState definitions

"""

import numpy

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def GetPrincipleQuantumNumber(atNum):
  """ Get principal quantum number for atom number """
  if atNum <= 2:
    return 1
  if atNum <= 10:
    return 2
  if atNum <= 18:
    return 3
  if atNum <= 36:
    return 4
  if atNum <= 54:
    return 5
  if atNum <= 86:
    return 6
  return 7


def EStateIndices(mol, force=True):
  """ returns a tuple of EState indices for the molecule

    Reference: Hall, Mohney and Kier. JCICS _31_ 76-81 (1991)

  """
  if not force and hasattr(mol, '_eStateIndices'):
    return mol._eStateIndices
  res = numpy.array(rdMolDescriptors.CalcEStateIndices(mol), dtype=numpy.float64)
  mol._eStateIndices = res
  return res


EStateIndices.version = '1.0.0'


def MaxEStateIndex(mol, force=1):
  return rdMolDescriptors.CalcMaxEStateIndex(mol)


MaxEStateIndex.version = "1.0.0"


def MinEStateIndex(mol, force=1):
  return rdMolDescriptors.CalcMinEStateIndex(mol)


MinEStateIndex.version = "1.0.0"


def MaxAbsEStateIndex(mol, force=1):
  return rdMolDescriptors.CalcMaxAbsEStateIndex(mol)


MaxAbsEStateIndex.version = "1.0.0"


def MinAbsEStateIndex(mol, force=1):
  return rdMolDescriptors.CalcMinAbsEStateIndex(mol)


MinAbsEStateIndex.version = "1.0.0"


def _exampleCode():
  """ Example code for calculating E-state indices """
  smis = ['CCCC', 'CCCCC', 'CCCCCC', 'CC(N)C(=O)O', 'CC(N)C(=O)[O-].[Na+]']
  for smi in smis:
    m = Chem.MolFromSmiles(smi)
    print(smi)
    inds = EStateIndices(m)
    print('\t', inds)


if __name__ == '__main__':  # pragma: nocover
  _exampleCode()
