# $Id$
#
#  Copyright (C) 2004-2007 Greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" Contains an implementation of Atom-pair fingerprints, as
described in:

R.E. Carhart, D.H. Smith, R. Venkataraghavan;
"Atom Pairs as Molecular Features in Structure-Activity Studies:
Definition and Applications" JCICS 25, 64-73 (1985).

"""
from DataStructs.SparseIntVect import SparseIntVect
import Chem
from Chem.AtomPairs import Utils
import DataStructs

_maxPathLen=31
numPathBits=5
numFpBits=numPathBits+2*Utils.codeSize
fpLen=1L<<numFpBits

def ScorePair(at1,at2,dist):
  """ Returns a score for an individual atom pair.

  >>> m = Chem.MolFromSmiles('CCCCC')
  >>> c1 = Utils.GetAtomCode(m.GetAtomWithIdx(0))
  >>> c2 = Utils.GetAtomCode(m.GetAtomWithIdx(1))
  >>> c3 = Utils.GetAtomCode(m.GetAtomWithIdx(2))
  >>> t = 1 | min(c1,c2)<<numPathBits | max(c1,c2)<<(Utils.codeSize+numPathBits)
  >>> ScorePair(m.GetAtomWithIdx(0),m.GetAtomWithIdx(1),1)==t
  1
  >>> ScorePair(m.GetAtomWithIdx(1),m.GetAtomWithIdx(0),1)==t
  1
  >>> t = 2 | min(c1,c3)<<numPathBits | max(c1,c3)<<(Utils.codeSize+numPathBits)
  >>> ScorePair(m.GetAtomWithIdx(0),m.GetAtomWithIdx(2),2)==t
  1

  """
  code1 = Utils.GetAtomCode(at1)
  code2 = Utils.GetAtomCode(at2)
  accum = dist % _maxPathLen
  accum |= min(code1,code2) << numPathBits
  accum |= max(code1,code2) << (Utils.codeSize+numPathBits)
  return accum

def ExplainPairScore(score):
  """ 
  >>> m = Chem.MolFromSmiles('C=CC')
  >>> score = ScorePair(m.GetAtomWithIdx(0),m.GetAtomWithIdx(1),1)
  >>> ExplainPairScore(score)
  (('C', 1, 1), 1, ('C', 2, 1))
  >>> score = ScorePair(m.GetAtomWithIdx(0),m.GetAtomWithIdx(2),2)
  >>> ExplainPairScore(score)
  (('C', 1, 0), 2, ('C', 1, 1))
  >>> score = ScorePair(m.GetAtomWithIdx(1),m.GetAtomWithIdx(2),1)
  >>> ExplainPairScore(score)
  (('C', 1, 0), 1, ('C', 2, 1))
  >>> score = ScorePair(m.GetAtomWithIdx(2),m.GetAtomWithIdx(1),1)
  >>> ExplainPairScore(score)
  (('C', 1, 0), 1, ('C', 2, 1))

  """
  codeMask = (1<<Utils.codeSize)-1
  pathMask = (1<<numPathBits)-1
  dist = score&pathMask

  score = score>>numPathBits
  code1 = score&codeMask
  score = score>>Utils.codeSize
  code2 = score&codeMask

  res = Utils.ExplainAtomCode(code1),dist,Utils.ExplainAtomCode(code2)
  return res

def GetAtomPairFingerprintAsCounts(mol):
  """ Returns the Atom-pair fingerprint for a molecule as
  a tuple of on-bit IDs.

  **Arguments**:

    - mol: a molecule

  **Returns**: a tuple of ints

  >>> m = Chem.MolFromSmiles('CCC')
  >>> v = [ ScorePair(m.GetAtomWithIdx(0),m.GetAtomWithIdx(1),1),
  ...       ScorePair(m.GetAtomWithIdx(0),m.GetAtomWithIdx(2),2),
  ...       ScorePair(m.GetAtomWithIdx(1),m.GetAtomWithIdx(2),1),
  ...     ]
  >>> v.sort()
  >>> v = tuple(v)
  >>> GetAtomPairFingerprintAsCounts(m)==v
  1
  
  """
  res = []
  dists = Chem.GetDistanceMatrix(mol)
  nAtoms = mol.GetNumAtoms()
  for i in range(nAtoms):
    at1 = mol.GetAtomWithIdx(i)
    for j in range(i+1,nAtoms):
      dist = dists[i,j]
      if dist<_maxPathLen:
        at2 = mol.GetAtomWithIdx(j)
      res.append(ScorePair(at1,at2,int(dist)))
  res.sort()
  return tuple(res)

def GetAtomPairFingerprintAsIntVect(mol):
  """ Returns the Atom-pair fingerprint for a molecule as
  a SparseIntVect

  **Arguments**:

    - mol: a molecule

  **Returns**: a SparseIntVect


  >>> m = Chem.MolFromSmiles('CCC')
  >>> v = [ ScorePair(m.GetAtomWithIdx(0),m.GetAtomWithIdx(1),1),
  ...       ScorePair(m.GetAtomWithIdx(0),m.GetAtomWithIdx(2),2),
  ...       ScorePair(m.GetAtomWithIdx(1),m.GetAtomWithIdx(2),1),
  ...     ]
  >>> val = SparseIntVect(fpLen)
  >>> val.InitFromSequence(v)
  >>> fp = GetAtomPairFingerprintAsIntVect(m)
  >>> val==fp
  True
  
  """
  res = SparseIntVect(fpLen)
  counts = GetAtomPairFingerprintAsCounts(mol)
  for v in counts:
    res[v]+=1
  return res

def GetAtomPairFingerprintAsBitVect(mol):
  """ Returns the Atom-pair fingerprint for a molecule as
  a SparseBitVect. Note that this doesn't match the standard
  definition of atom pairs, which uses counts of the
  pairs, not just their presence.

  **Arguments**:

    - mol: a molecule

  **Returns**: a SparseBitVect

  >>> m = Chem.MolFromSmiles('CCC')
  >>> v = [ ScorePair(m.GetAtomWithIdx(0),m.GetAtomWithIdx(1),1),
  ...       ScorePair(m.GetAtomWithIdx(0),m.GetAtomWithIdx(2),2),
  ...     ]
  >>> v.sort()
  >>> fp = GetAtomPairFingerprintAsBitVect(m)
  >>> list(fp.GetOnBits())==v
  True
  
  """
  if numFpBits>31:
    raise ValueError,"fingerprint is too long to be encoded as a BitVect."
  res = DataStructs.SparseBitVect(fpLen)
  counts = GetAtomPairFingerprintAsCounts(mol)
  for idx in counts:
    res.SetBit(idx)
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
  
  

