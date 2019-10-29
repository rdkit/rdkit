# $Id$
#
# Copyright (C) 2002-2006 greg Landrum and Rational Discovery LLC
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


def GetPrincipleQuantumNumber(atNum):
  """ Get principal quantum number for atom number """
  if atNum <= 2:
    return 1
  elif atNum <= 10:
    return 2
  elif atNum <= 18:
    return 3
  elif atNum <= 36:
    return 4
  elif atNum <= 54:
    return 5
  elif atNum <= 86:
    return 6
  else:
    return 7


def EStateIndices(mol, force=True):
  """ returns a tuple of EState indices for the molecule

    Reference: Hall, Mohney and Kier. JCICS _31_ 76-81 (1991)

  """
  if not force and hasattr(mol, '_eStateIndices'):
    return mol._eStateIndices

  tbl = Chem.GetPeriodicTable()
  nAtoms = mol.GetNumAtoms()
  Is = numpy.zeros(nAtoms, numpy.float)
  for i in range(nAtoms):
    at = mol.GetAtomWithIdx(i)
    atNum = at.GetAtomicNum()
    d = at.GetDegree()
    if d > 0:
      h = at.GetTotalNumHs()
      dv = tbl.GetNOuterElecs(atNum) - h
      N = GetPrincipleQuantumNumber(atNum)
      Is[i] = (4. / (N * N) * dv + 1) / d
  dists = Chem.GetDistanceMatrix(mol, useBO=0, useAtomWts=0)
  dists += 1
  accum = numpy.zeros(nAtoms, numpy.float)
  for i in range(nAtoms):
    for j in range(i + 1, nAtoms):
      p = dists[i, j]
      if p < 1e6:
        tmp = (Is[i] - Is[j]) / (p * p)
        accum[i] += tmp
        accum[j] -= tmp

  res = accum + Is
  mol._eStateIndices = res
  return res


EStateIndices.version = '1.0.0'


def MaxEStateIndex(mol, force=1):
  return max(EStateIndices(mol, force))


MaxEStateIndex.version = "1.0.0"


def MinEStateIndex(mol, force=1):
  return min(EStateIndices(mol, force))


MinEStateIndex.version = "1.0.0"


def MaxAbsEStateIndex(mol, force=1):
  return max([abs(x) for x in EStateIndices(mol, force)])


MaxAbsEStateIndex.version = "1.0.0"


def MinAbsEStateIndex(mol, force=1):
  return min([abs(x) for x in EStateIndices(mol, force)])


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
