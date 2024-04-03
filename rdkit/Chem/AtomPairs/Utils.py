#
#  Copyright (C) 2004-2017 Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

import math
import warnings

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def ExplainAtomCode(code, branchSubtract=0, includeChirality=False):
  """

    **Arguments**:

      - the code to be considered

      - branchSubtract: (optional) the constant that was subtracted off
        the number of neighbors before integrating it into the code.
        This is used by the topological torsions code.

      - includeChirality: (optional) Determines whether or not chirality
        was included when generating the atom code.

    >>> m = Chem.MolFromSmiles('C=CC(=O)O')
    >>> code = GetAtomCode(m.GetAtomWithIdx(0))
    >>> ExplainAtomCode(code)
    ('C', 1, 1)
    >>> code = GetAtomCode(m.GetAtomWithIdx(1))
    >>> ExplainAtomCode(code)
    ('C', 2, 1)
    >>> code = GetAtomCode(m.GetAtomWithIdx(2))
    >>> ExplainAtomCode(code)
    ('C', 3, 1)
    >>> code = GetAtomCode(m.GetAtomWithIdx(3))
    >>> ExplainAtomCode(code)
    ('O', 1, 1)
    >>> code = GetAtomCode(m.GetAtomWithIdx(4))
    >>> ExplainAtomCode(code)
    ('O', 1, 0)

    we can do chirality too, that returns an extra element in the tuple:

    >>> m = Chem.MolFromSmiles('C[C@H](F)Cl')
    >>> code = GetAtomCode(m.GetAtomWithIdx(1))
    >>> ExplainAtomCode(code)
    ('C', 3, 0)
    >>> code = GetAtomCode(m.GetAtomWithIdx(1),includeChirality=True)
    >>> ExplainAtomCode(code,includeChirality=True)
    ('C', 3, 0, 'R')

    note that if we don't ask for chirality, we get the right answer even if
    the atom code was calculated with chirality:

    >>> ExplainAtomCode(code)
    ('C', 3, 0)

    non-chiral atoms return '' in the 4th field:

    >>> code = GetAtomCode(m.GetAtomWithIdx(0),includeChirality=True)
    >>> ExplainAtomCode(code,includeChirality=True)
    ('C', 1, 0, '')

    Obviously switching the chirality changes the results:

    >>> m = Chem.MolFromSmiles('C[C@@H](F)Cl')
    >>> code = GetAtomCode(m.GetAtomWithIdx(1),includeChirality=True)
    >>> ExplainAtomCode(code,includeChirality=True)
    ('C', 3, 0, 'S')

    """
  typeMask = (1 << rdMolDescriptors.AtomPairsParameters.numTypeBits) - 1
  branchMask = (1 << rdMolDescriptors.AtomPairsParameters.numBranchBits) - 1
  piMask = (1 << rdMolDescriptors.AtomPairsParameters.numPiBits) - 1
  chiMask = (1 << rdMolDescriptors.AtomPairsParameters.numChiralBits) - 1

  nBranch = int(code & branchMask)
  code = code >> rdMolDescriptors.AtomPairsParameters.numBranchBits

  nPi = int(code & piMask)
  code = code >> rdMolDescriptors.AtomPairsParameters.numPiBits

  typeIdx = int(code & typeMask)
  if typeIdx < len(rdMolDescriptors.AtomPairsParameters.atomTypes):
    atomNum = rdMolDescriptors.AtomPairsParameters.atomTypes[typeIdx]
    atomSymbol = Chem.GetPeriodicTable().GetElementSymbol(atomNum)
  else:
    atomSymbol = 'X'

  if not includeChirality:
    return (atomSymbol, nBranch, nPi)

  code = code >> rdMolDescriptors.AtomPairsParameters.numTypeBits
  chiDict = {0: '', 1: 'R', 2: 'S'}
  chiCode = int(code & chiMask)
  return (atomSymbol, nBranch, nPi, chiDict[chiCode])


GetAtomCode = rdMolDescriptors.GetAtomPairAtomCode


def NumPiElectrons(atom):
  """ DEPRECATED: please use rdMolDescriptors.GetNumPiElectrons instead.
  """

  warnings.warn(
    "The Chem.AtomPairs.Utils.NumPiElectrons is deprecated; please use rdMolDescriptors.GetNumPiElectrons instead.",
    DeprecationWarning, stacklevel=2)

  res = 0
  if atom.GetIsAromatic():
    res = 1
  elif atom.GetHybridization() != Chem.HybridizationType.SP3:
    # the number of pi electrons is just the number of
    # unsaturations (valence - degree):
    res = atom.GetExplicitValence() - atom.GetNumExplicitHs()
    if res < atom.GetDegree():
      raise ValueError("explicit valence exceeds atom degree")
    res -= atom.GetDegree()
  return res


def BitsInCommon(v1, v2):
  """ Returns the number of bits in common between two vectors

    **Arguments**:

      - two vectors (sequences of bit ids)

    **Returns**: an integer

    **Notes**

      - the vectors must be sorted

      - duplicate bit IDs are counted more than once

    >>> BitsInCommon( (1,2,3,4,10), (2,4,6) )
    2

    Here's how duplicates are handled:

    >>> BitsInCommon( (1,2,2,3,4), (2,2,4,5,6) )
    3

    """
  res = 0
  v2Pos = 0
  nV2 = len(v2)
  for val in v1:
    while v2Pos < nV2 and v2[v2Pos] < val:
      v2Pos += 1
    if v2Pos >= nV2:
      break
    if v2[v2Pos] == val:
      res += 1
      v2Pos += 1
  return res


def DiceSimilarity(v1, v2, bounds=None):
  """ Implements the DICE similarity metric.
     This is the recommended metric in both the Topological torsions
     and Atom pairs papers.

    **Arguments**:

      - two vectors (sequences of bit ids)

    **Returns**: a float.

    **Notes**

      - the vectors must be sorted


    >>> DiceSimilarity( (1,2,3), (1,2,3) )
    1.0
    >>> DiceSimilarity( (1,2,3), (5,6) )
    0.0
    >>> DiceSimilarity( (1,2,3,4), (1,3,5,7) )
    0.5
    >>> DiceSimilarity( (1,2,3,4,5,6), (1,3) )
    0.5

    Note that duplicate bit IDs count multiple times:

    >>> DiceSimilarity( (1,1,3,4,5,6), (1,1) )
    0.5

    but only if they are duplicated in both vectors:

    >>> DiceSimilarity( (1,1,3,4,5,6), (1,) )==2./7
    True

    edge case

    >>> DiceSimilarity( (), () )
    0.0

    and bounds check

    >>> DiceSimilarity( (1,1,3,4), (1,1))
    0.666...
    >>> DiceSimilarity( (1,1,3,4), (1,1), bounds=0.3)
    0.666...
    >>> DiceSimilarity( (1,1,3,4), (1,1), bounds=0.33)
    0.666...
    >>> DiceSimilarity( (1,1,3,4,5,6), (1,1), bounds=0.34)
    0.0

    """
  denom = 1.0 * (len(v1) + len(v2))
  if not denom:
    res = 0.0
  else:
    if bounds and (min(len(v1), len(v2)) / denom) < bounds:
      numer = 0.0
    else:
      numer = 2.0 * BitsInCommon(v1, v2)
    res = numer / denom
  return res


def Dot(v1, v2):
  """ Returns the Dot product between two vectors:

    **Arguments**:

      - two vectors (sequences of bit ids)

    **Returns**: an integer

    **Notes**

      - the vectors must be sorted

      - duplicate bit IDs are counted more than once

    >>> Dot( (1,2,3,4,10), (2,4,6) )
    2

    Here's how duplicates are handled:

    >>> Dot( (1,2,2,3,4), (2,2,4,5,6) )
    5
    >>> Dot( (1,2,2,3,4), (2,4,5,6) )
    2
    >>> Dot( (1,2,2,3,4), (5,6) )
    0
    >>> Dot( (), (5,6) )
    0

    """
  res = 0
  nV1 = len(v1)
  nV2 = len(v2)
  i = 0
  j = 0
  while i < nV1:
    v1Val = v1[i]
    v1Count = 1
    i += 1
    while i < nV1 and v1[i] == v1Val:
      v1Count += 1
      i += 1
    while j < nV2 and v2[j] < v1Val:
      j += 1
    if j < nV2 and v2[j] == v1Val:
      v2Count = 1
      j += 1
      while j < nV2 and v2[j] == v1Val:
        v2Count += 1
        j += 1
      commonCount = min(v1Count, v2Count)
      res += commonCount * commonCount
    elif j >= nV2:
      break
  return res


def CosineSimilarity(v1, v2):
  """ Implements the Cosine similarity metric.
     This is the recommended metric in the LaSSI paper

    **Arguments**:

      - two vectors (sequences of bit ids)

    **Returns**: a float.

    **Notes**

      - the vectors must be sorted

    >>> print('%.3f'%CosineSimilarity( (1,2,3,4,10), (2,4,6) ))
    0.516
    >>> print('%.3f'%CosineSimilarity( (1,2,2,3,4), (2,2,4,5,6) ))
    0.714
    >>> print('%.3f'%CosineSimilarity( (1,2,2,3,4), (1,2,2,3,4) ))
    1.000
    >>> print('%.3f'%CosineSimilarity( (1,2,2,3,4), (5,6,7) ))
    0.000
    >>> print('%.3f'%CosineSimilarity( (1,2,2,3,4), () ))
    0.000

    """
  d1 = Dot(v1, v1)
  d2 = Dot(v2, v2)
  denom = math.sqrt(d1 * d2)
  if not denom:
    res = 0.0
  else:
    numer = Dot(v1, v2)
    res = numer / denom
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
