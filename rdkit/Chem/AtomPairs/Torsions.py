# $Id$
#
#  Copyright (C) 2004-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""
Contains an implementation of Topological-torsion fingerprints, as
described in:

R. Nilakantan, N. Bauman, J. S. Dixon, R. Venkataraghavan;
"Topological Torsion: A New Molecular Descriptor for SAR Applications.
Comparison with Other Descriptors" JCICS 27, 82-85 (1987).

The fingerprints can be accessed through the following functions:
- GetTopologicalTorsionFingerprint
- GetHashedTopologicalTorsionFingerprint
- GetTopologicalTorsionFingerprintAsIntVect (identical to GetTopologicalTorsionFingerprint)
- GetTopologicalTorsionFingerprintAsIds

"""
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.AtomPairs import Utils
from rdkit.Chem.rdMolDescriptors import (
    GetHashedTopologicalTorsionFingerprint, GetTopologicalTorsionFingerprint)

GetTopologicalTorsionFingerprintAsIntVect = rdMolDescriptors.GetTopologicalTorsionFingerprint


def pyScorePath(mol, path, size, atomCodes=None):
  """ Returns a score for an individual path.

  >>> from rdkit import Chem
  >>> m = Chem.MolFromSmiles('CCCCC')
  >>> c1 = Utils.GetAtomCode(m.GetAtomWithIdx(0), 1)
  >>> c2 = Utils.GetAtomCode(m.GetAtomWithIdx(1), 2)
  >>> c3 = Utils.GetAtomCode(m.GetAtomWithIdx(2), 2)
  >>> c4 = Utils.GetAtomCode(m.GetAtomWithIdx(3), 1)
  >>> t = c1 | (c2 << rdMolDescriptors.AtomPairsParameters.codeSize) | (c3 << (rdMolDescriptors.AtomPairsParameters.codeSize * 2)) | (c4 << (rdMolDescriptors.AtomPairsParameters.codeSize * 3))
  >>> pyScorePath(m, (0, 1, 2, 3), 4) == t
  1

  The scores are path direction independent:

  >>> pyScorePath(m, (3, 2, 1, 0), 4) == t
  1

  >>> m = Chem.MolFromSmiles('C=CC(=O)O')
  >>> c1 = Utils.GetAtomCode(m.GetAtomWithIdx(0), 1)
  >>> c2 = Utils.GetAtomCode(m.GetAtomWithIdx(1), 2)
  >>> c3 = Utils.GetAtomCode(m.GetAtomWithIdx(2), 2)
  >>> c4 = Utils.GetAtomCode(m.GetAtomWithIdx(4), 1)
  >>> t = c1 | (c2 << rdMolDescriptors.AtomPairsParameters.codeSize) | (c3 << (rdMolDescriptors.AtomPairsParameters.codeSize * 2)) | (c4 << (rdMolDescriptors.AtomPairsParameters.codeSize * 3))
  >>> pyScorePath(m, (0, 1, 2, 4), 4) == t
  1

  """
  codes = [None] * size
  for i in range(size):
    if i == 0 or i == size - 1:
      sub = 1
    else:
      sub = 2
    if not atomCodes:
      codes[i] = Utils.GetAtomCode(mol.GetAtomWithIdx(path[i]), sub)
    else:
      codes[i] = atomCodes[path[i]] - sub

  # "canonicalize" the code vector:
  beg = 0
  end = len(codes) - 1
  while beg < end:
    if codes[beg] == codes[end]:
      beg += 1
      end -= 1
    else:
      if codes[beg] > codes[end]:
        codes.reverse()
      break

  accum = 0
  codeSize = rdMolDescriptors.AtomPairsParameters.codeSize
  for i, code in enumerate(codes):
    accum |= code << (codeSize * i)
  return accum


def ExplainPathScore(score, size=4):
  """

  >>> from rdkit import Chem
  >>> m = Chem.MolFromSmiles('C=CC')
  >>> score=pyScorePath(m, (0, 1, 2), 3)
  >>> ExplainPathScore(score, 3)
  (('C', 1, 0), ('C', 2, 1), ('C', 1, 1))

  Again, it's order independent:

  >>> score = pyScorePath(m, (2, 1, 0), 3)
  >>> ExplainPathScore(score, 3)
  (('C', 1, 0), ('C', 2, 1), ('C', 1, 1))

  >>> m = Chem.MolFromSmiles('C=CO')
  >>> score=pyScorePath(m, (0, 1, 2), 3)
  >>> ExplainPathScore(score, 3)
  (('C', 1, 1), ('C', 2, 1), ('O', 1, 0))

  >>> m = Chem.MolFromSmiles('OC=CO')
  >>> score=pyScorePath(m, (0, 1, 2, 3), 4)
  >>> ExplainPathScore(score, 4)
  (('O', 1, 0), ('C', 2, 1), ('C', 2, 1), ('O', 1, 0))

  >>> m = Chem.MolFromSmiles('CC=CO')
  >>> score=pyScorePath(m, (0, 1, 2, 3), 4)
  >>> ExplainPathScore(score, 4)
  (('C', 1, 0), ('C', 2, 1), ('C', 2, 1), ('O', 1, 0))


  >>> m = Chem.MolFromSmiles('C=CC(=O)O')
  >>> score=pyScorePath(m, (0, 1, 2, 3), 4)
  >>> ExplainPathScore(score, 4)
  (('C', 1, 1), ('C', 2, 1), ('C', 3, 1), ('O', 1, 1))
  >>> score=pyScorePath(m, (0, 1, 2, 4), 4)
  >>> ExplainPathScore(score, 4)
  (('C', 1, 1), ('C', 2, 1), ('C', 3, 1), ('O', 1, 0))


  >>> m = Chem.MolFromSmiles('OOOO')
  >>> score=pyScorePath(m, (0, 1, 2), 3)
  >>> ExplainPathScore(score, 3)
  (('O', 1, 0), ('O', 2, 0), ('O', 2, 0))
  >>> score=pyScorePath(m, (0, 1, 2, 3), 4)
  >>> ExplainPathScore(score, 4)
  (('O', 1, 0), ('O', 2, 0), ('O', 2, 0), ('O', 1, 0))

  """
  codeSize: int = rdMolDescriptors.AtomPairsParameters.codeSize
  codeMask = (1 << codeSize) - 1
  res = [None] * size
  for i in range(size):
    if i == 0 or i == size - 1:
      sub = 1
    else:
      sub = 2
    code = score & codeMask
    score = score >> codeSize
    symb, nBranch, nPi = Utils.ExplainAtomCode(code)
    res[i] = (symb, nBranch + sub, nPi)
  return tuple(res)


def GetTopologicalTorsionFingerprintAsIds(mol, targetSize=4):
  nonZeroElements = GetTopologicalTorsionFingerprint(mol, targetSize).GetNonzeroElements()
  frequencies = sorted(nonZeroElements.items())
  res = []
  for k, v in frequencies:
    res.extend([k] * v)
  return res


# ------------------------------------
#
#  doctest boilerplate
#
def _runDoctests(verbose=None):  # pragma: nocover
  import doctest
  import sys
  failed, _ = doctest.testmod(optionflags=doctest.ELLIPSIS, verbose=verbose)
  sys.exit(failed)


if __name__ == '__main__':  # pragma: nocover
  _runDoctests()
