# $Id$
#
#  Copyright (C) 2004-2007 Greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Contains an implementation of Atom-pair fingerprints, as
described in:

R.E. Carhart, D.H. Smith, R. Venkataraghavan;
"Atom Pairs as Molecular Features in Structure-Activity Studies:
Definition and Applications" JCICS 25, 64-73 (1985).

"""
from rdkit.DataStructs import IntSparseIntVect
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.AtomPairs import Utils
from rdkit import DataStructs

from rdkit.Chem.rdMolDescriptors import GetAtomPairFingerprint,GetHashedAtomPairFingerprint
GetAtomPairFingerprintAsIntVect=rdMolDescriptors.GetAtomPairFingerprint

numPathBits=rdMolDescriptors.AtomPairsParameters.numPathBits
_maxPathLen=(1<<numPathBits)-1
numFpBits=numPathBits+2*rdMolDescriptors.AtomPairsParameters.codeSize
fpLen=1<<numFpBits

def pyScorePair(at1,at2,dist,atomCodes=None):
  """ Returns a score for an individual atom pair.

  >>> m = Chem.MolFromSmiles('CCCCC')
  >>> c1 = Utils.GetAtomCode(m.GetAtomWithIdx(0))
  >>> c2 = Utils.GetAtomCode(m.GetAtomWithIdx(1))
  >>> c3 = Utils.GetAtomCode(m.GetAtomWithIdx(2))
  >>> t = 1 | min(c1,c2)<<numPathBits | max(c1,c2)<<(rdMolDescriptors.AtomPairsParameters.codeSize+numPathBits)
  >>> pyScorePair(m.GetAtomWithIdx(0),m.GetAtomWithIdx(1),1)==t
  1
  >>> pyScorePair(m.GetAtomWithIdx(1),m.GetAtomWithIdx(0),1)==t
  1
  >>> t = 2 | min(c1,c3)<<numPathBits | max(c1,c3)<<(rdMolDescriptors.AtomPairsParameters.codeSize+numPathBits)
  >>> pyScorePair(m.GetAtomWithIdx(0),m.GetAtomWithIdx(2),2)==t
  1
  >>> pyScorePair(m.GetAtomWithIdx(0),m.GetAtomWithIdx(2),2,
  ...  atomCodes=(Utils.GetAtomCode(m.GetAtomWithIdx(0)),Utils.GetAtomCode(m.GetAtomWithIdx(2))))==t
  1

  """
  if not atomCodes:
    code1 = Utils.GetAtomCode(at1)
    code2 = Utils.GetAtomCode(at2)
  else:
    code1,code2=atomCodes
  accum = int(dist) % _maxPathLen
  accum |= min(code1,code2) << numPathBits
  accum |= max(code1,code2) << (rdMolDescriptors.AtomPairsParameters.codeSize+numPathBits)
  return accum

def ExplainPairScore(score):
  """ 
  >>> m = Chem.MolFromSmiles('C=CC')
  >>> score = pyScorePair(m.GetAtomWithIdx(0),m.GetAtomWithIdx(1),1)
  >>> ExplainPairScore(score)
  (('C', 1, 1), 1, ('C', 2, 1))
  >>> score = pyScorePair(m.GetAtomWithIdx(0),m.GetAtomWithIdx(2),2)
  >>> ExplainPairScore(score)
  (('C', 1, 0), 2, ('C', 1, 1))
  >>> score = pyScorePair(m.GetAtomWithIdx(1),m.GetAtomWithIdx(2),1)
  >>> ExplainPairScore(score)
  (('C', 1, 0), 1, ('C', 2, 1))
  >>> score = pyScorePair(m.GetAtomWithIdx(2),m.GetAtomWithIdx(1),1)
  >>> ExplainPairScore(score)
  (('C', 1, 0), 1, ('C', 2, 1))

  """
  codeMask = (1<<rdMolDescriptors.AtomPairsParameters.codeSize)-1
  pathMask = (1<<numPathBits)-1
  dist = score&pathMask

  score = score>>numPathBits
  code1 = score&codeMask
  score = score>>rdMolDescriptors.AtomPairsParameters.codeSize
  code2 = score&codeMask

  res = Utils.ExplainAtomCode(code1),dist,Utils.ExplainAtomCode(code2)
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
  >>> v = [ pyScorePair(m.GetAtomWithIdx(0),m.GetAtomWithIdx(1),1),
  ...       pyScorePair(m.GetAtomWithIdx(0),m.GetAtomWithIdx(2),2),
  ...     ]
  >>> v.sort()
  >>> fp = GetAtomPairFingerprintAsBitVect(m)
  >>> list(fp.GetOnBits())==v
  True
  
  """
  res = DataStructs.SparseBitVect(fpLen)
  fp = rdMolDescriptors.GetAtomPairFingerprint(mol)
  for val in fp.GetNonzeroElements().keys():
    res.SetBit(val)
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
  
  

