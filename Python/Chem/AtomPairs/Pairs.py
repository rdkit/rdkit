# $Id: Pairs.py 5007 2006-02-22 15:14:41Z glandrum $
#
#  Copyright (C) 2004-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" Contains an implementation of Atom-pair fingerprints, as
described in:

R.E. Carhart, D.H. Smith, R. Venkataraghavan;
"Atom Pairs as Molecular Features in Structure-Activity Studies:
Definition and Applications" JCICS 25, 64-73 (1985).

"""
import Chem
import Utils

_maxPathLen=31
numPathBits=5

def ScorePair(at1,at2,dist):
  """ Returns a score for an individual atom pair.

  >>> m = Chem.MolFromSmiles('CCCCC')
  >>> c1 = Utils.GetAtomCode(m.GetAtomWithIdx(0))
  >>> c2 = Utils.GetAtomCode(m.GetAtomWithIdx(1))
  >>> c3 = Utils.GetAtomCode(m.GetAtomWithIdx(2))
  >>> t = 1 | min(c1,c2)<<Utils.codeSize | max(c1,c2)<<(2*Utils.codeSize)
  >>> ScorePair(m.GetAtomWithIdx(0),m.GetAtomWithIdx(1),1)==t
  1
  >>> t = 2 | min(c1,c3)<<Utils.codeSize | max(c1,c3)<<(2*Utils.codeSize)
  >>> ScorePair(m.GetAtomWithIdx(0),m.GetAtomWithIdx(2),2)==t
  1

  """
  code1 = Utils.GetAtomCode(at1)
  code2 = Utils.GetAtomCode(at2)
  accum = dist % _maxPathLen
  accum |= min(code1,code2) << Utils.codeSize
  accum |= max(code1,code2) << 2*Utils.codeSize
  return accum

def GetAtomPairFingerprint(mol):
  """ Returns the Atom-pair fingerprint for a molecule as
  a tuple of on-bit IDs.

  **Arguments**:

    - mol: a molecule

  **Returns**: a sorted tuple of integers and long integers

  >>> m = Chem.MolFromSmiles('CCC')
  >>> v = [ ScorePair(m.GetAtomWithIdx(0),m.GetAtomWithIdx(1),1),
  ...       ScorePair(m.GetAtomWithIdx(0),m.GetAtomWithIdx(2),2),
  ...       ScorePair(m.GetAtomWithIdx(1),m.GetAtomWithIdx(2),1),
  ...     ]
  >>> v.sort()
  >>> v = tuple(v)
  >>> GetAtomPairFingerprint(m)==v
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
  
  

