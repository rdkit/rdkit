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
from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import math

def ExplainAtomCode(code,branchSubtract=0):
  """

  **Arguments**:

    - the code to be considered

    - branchSubtract: (optional) the constant that was subtracted off
      the number of neighbors before integrating it into the code.  
      This is used by the topological torsions code.
      

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
  
  """
  typeMask = (1<<rdMolDescriptors.AtomPairsParameters.numTypeBits)-1
  branchMask = (1<<rdMolDescriptors.AtomPairsParameters.numBranchBits)-1
  piMask = (1<<rdMolDescriptors.AtomPairsParameters.numPiBits)-1

  nBranch = int(code&branchMask)
  #print(code,end='')
  code = code>>rdMolDescriptors.AtomPairsParameters.numBranchBits
  nPi = int(code&piMask)
  #print(code,end='')
  code = code>>rdMolDescriptors.AtomPairsParameters.numPiBits
  #print(code,end='')
  typeIdx=int(code&typeMask)
  if typeIdx<len(rdMolDescriptors.AtomPairsParameters.atomTypes):
    atomNum = rdMolDescriptors.AtomPairsParameters.atomTypes[typeIdx]
    atomSymbol=Chem.GetPeriodicTable().GetElementSymbol(atomNum)
  else:
    atomSymbol='X'
  return (atomSymbol,nBranch,nPi)

GetAtomCode=rdMolDescriptors.GetAtomPairAtomCode

def NumPiElectrons(atom):
  """ Returns the number of electrons an atom is using for pi bonding

  >>> m = Chem.MolFromSmiles('C=C')
  >>> NumPiElectrons(m.GetAtomWithIdx(0))
  1

  >>> m = Chem.MolFromSmiles('C#CC')
  >>> NumPiElectrons(m.GetAtomWithIdx(0))
  2
  >>> NumPiElectrons(m.GetAtomWithIdx(1))
  2

  >>> m = Chem.MolFromSmiles('O=C=CC')
  >>> NumPiElectrons(m.GetAtomWithIdx(0))
  1
  >>> NumPiElectrons(m.GetAtomWithIdx(1))
  2
  >>> NumPiElectrons(m.GetAtomWithIdx(2))
  1
  >>> NumPiElectrons(m.GetAtomWithIdx(3))
  0

  FIX: this behaves oddly in these cases:
  >>> m = Chem.MolFromSmiles('S(=O)(=O)')
  >>> NumPiElectrons(m.GetAtomWithIdx(0))
  2

  >>> m = Chem.MolFromSmiles('S(=O)(=O)(O)O')
  >>> NumPiElectrons(m.GetAtomWithIdx(0))
  0

  In the second case, the S atom is tagged as sp3 hybridized.

  """
  
  res = 0
  if atom.GetIsAromatic():
    res = 1
  elif atom.GetHybridization() != Chem.HybridizationType.SP3:
    # the number of pi electrons is just the number of
    # unsaturations (valence - degree):
    res = atom.GetExplicitValence()  - atom.GetNumExplicitHs()
    if res<atom.GetDegree():
      raise ValueError("explicit valence exceeds atom degree")
    res -= atom.GetDegree()
  return res


def BitsInCommon(v1,v2):
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
    while v2Pos<nV2 and v2[v2Pos]<val:
      v2Pos+=1
    if v2Pos >= nV2:
      break
    if v2[v2Pos]==val:
      res += 1
      v2Pos += 1
  return res


def DiceSimilarity(v1,v2,bounds=None):
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

  """
  denom = 1.0*(len(v1)+len(v2))
  if not denom:
    res = 0.0
  else:
    if bounds and (min(len(v1),len(v2))/denom) < bounds:
      numer = 0.0
    else:
      numer = 2.0*BitsInCommon(v1,v2)
    res = numer/denom

  return res


def Dot(v1,v2):
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
  while i<nV1:
    v1Val = v1[i]
    v1Count = 1
    i+=1
    while i<nV1 and v1[i]==v1Val:
      v1Count += 1
      i+=1
    while j<nV2 and v2[j]<v1Val:
      j+=1
    if j < nV2 and v2[j]==v1Val:
      v2Val = v2[j]
      v2Count = 1
      j+=1
      while j<nV2 and v2[j]==v1Val:
        v2Count+=1
        j+=1
      commonCount=min(v1Count,v2Count)
      res += commonCount*commonCount
    elif j>=nV2:
      break
  return res

def CosineSimilarity(v1,v2):
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
  d1 = Dot(v1,v1)
  d2 = Dot(v2,v2)
  denom = math.sqrt(d1*d2)
  if not denom:
    res = 0.0
  else:
    numer = Dot(v1,v2)
    res = numer/denom
  return res




#------------------------------------
#
#  doctest boilerplate
#
def _test():
  import doctest,sys
  return doctest.testmod(sys.modules["__main__"])


if __name__ == '__main__':
  import sys
  failed,tried = _test()
  sys.exit(failed)
