# $Id$
#
#  Copyright (C) 2004-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" Contains an implementation of Topological-torsion fingerprints, as
described in:

R. Nilakantan, N. Bauman, J. S. Dixon, R. Venkataraghavan;
"Topological Torsion: A New Molecular Descriptor for SAR Applications.
Comparison with Other Descriptors" JCICS 27, 82-85 (1987).

"""
import warnings
#warnings.filterwarnings("ignore",'',FutureWarning)
import Chem
import Utils

def ScorePath(mol,path,size,atomCodes=None):
  """ Returns a score for an individual path.

  >>> m = Chem.MolFromSmiles('CCCCC')
  >>> c1 = long(Utils.GetAtomCode(m.GetAtomWithIdx(0),1))
  >>> c2 = long(Utils.GetAtomCode(m.GetAtomWithIdx(1),2))
  >>> c3 = long(Utils.GetAtomCode(m.GetAtomWithIdx(2),2))
  >>> c4 = long(Utils.GetAtomCode(m.GetAtomWithIdx(3),1))
  >>> t = c1 | (c2 << Utils.codeSize) | (c3 << (Utils.codeSize*2)) | (c4 << (Utils.codeSize*3))
  >>> ScorePath(m,(0,1,2,3),4)==t
  1

  The scores are path direction independent:
  >>> ScorePath(m,(3,2,1,0),4)==t
  1


  >>> m = Chem.MolFromSmiles('C=CC(=O)O')
  >>> c1 = long(Utils.GetAtomCode(m.GetAtomWithIdx(0),1))
  >>> c2 = long(Utils.GetAtomCode(m.GetAtomWithIdx(1),2))
  >>> c3 = long(Utils.GetAtomCode(m.GetAtomWithIdx(2),2))
  >>> c4 = long(Utils.GetAtomCode(m.GetAtomWithIdx(4),1))
  >>> t = c1 | (c2 << Utils.codeSize) | (c3 << (Utils.codeSize*2)) | (c4 << (Utils.codeSize*3))
  >>> ScorePath(m,(0,1,2,4),4)==t
  1

  """
  codes = [None]*size
  for i in range(size):
    if i==0 or i==(size-1):
      sub = 1
    else:
      sub = 2
    if not atomCodes:
      codes[i] = Utils.GetAtomCode(mol.GetAtomWithIdx(path[i]),sub)
    else:
      base = atomCodes[path[i]]
      codes[i]=base-sub

  # "canonicalize" the code vector:
  beg=0
  end = len(codes)-1
  while(beg < end):
    if codes[beg] > codes[end]:
      codes.reverse()
      break
    elif codes[beg]==codes[end]:
      beg += 1
      end -= 1
    else:
      break
  accum = 0L
  for i in range(size):
    accum |= long(codes[i]) << (Utils.codeSize*i)
  return accum

def ExplainPathScore(score,size=4):
  """ 

  >>> m = Chem.MolFromSmiles('C=CC')
  >>> score=ScorePath(m,(0,1,2),3)
  >>> ExplainPathScore(score,3)
  (('C', 1, 0), ('C', 2, 1), ('C', 1, 1))

  Again, it's order independent:
  >>> score=ScorePath(m,(2,1,0),3)
  >>> ExplainPathScore(score,3)
  (('C', 1, 0), ('C', 2, 1), ('C', 1, 1))


  >>> m = Chem.MolFromSmiles('C=CO')
  >>> score=ScorePath(m,(0,1,2),3)
  >>> ExplainPathScore(score,3)
  (('C', 1, 1), ('C', 2, 1), ('O', 1, 0))

  >>> m = Chem.MolFromSmiles('OC=CO')
  >>> score=ScorePath(m,(0,1,2,3),4)
  >>> ExplainPathScore(score,4)
  (('O', 1, 0), ('C', 2, 1), ('C', 2, 1), ('O', 1, 0))

  >>> m = Chem.MolFromSmiles('CC=CO')
  >>> score=ScorePath(m,(0,1,2,3),4)
  >>> ExplainPathScore(score,4)
  (('C', 1, 0), ('C', 2, 1), ('C', 2, 1), ('O', 1, 0))


  >>> m = Chem.MolFromSmiles('C=CC(=O)O')
  >>> score=ScorePath(m,(0,1,2,3),4)
  >>> ExplainPathScore(score,4)
  (('C', 1, 1), ('C', 2, 1), ('C', 3, 1), ('O', 1, 1))
  >>> score=ScorePath(m,(0,1,2,4),4)
  >>> ExplainPathScore(score,4)
  (('C', 1, 1), ('C', 2, 1), ('C', 3, 1), ('O', 1, 0))


  >>> m = Chem.MolFromSmiles('OOOO')
  >>> score=ScorePath(m,(0,1,2),3)
  >>> ExplainPathScore(score,3)
  (('O', 1, 0), ('O', 2, 0), ('O', 2, 0))
  >>> score=ScorePath(m,(0,1,2,3),4)
  >>> ExplainPathScore(score,4)
  (('O', 1, 0), ('O', 2, 0), ('O', 2, 0), ('O', 1, 0))





  """
  codeMask=(1<<Utils.codeSize)-1
  res=[None]*size
  #print '>>>>>>>>>>>',score,size,codeMask
  for i in range(size):
    if i==0 or i==(size-1):
      sub = 1
    else:
      sub = 2
    code = score&codeMask
    #print i,code,score
    score = score>>Utils.codeSize
    symb,nBranch,nPi = Utils.ExplainAtomCode(code)
    expl = symb,nBranch+sub,nPi
    res[i] = expl
  return tuple(res)

def GetTopologicalTorsionFingerprint(mol,targetSize=4):
  """ Returns the topological torsion fingerprint for a molecule as
  a tuple of on-bit IDs.

  **Arguments**:

    - mol: a molecule

    - targetSize: targetSize (in atoms) of the torsions to find.

  **Returns**: a sorted tuple of integers and long integers


  >>> m = Chem.MolFromSmiles('CC(C)C')

  No torsions (no paths of length 4):
  >>> GetTopologicalTorsionFingerprint(m)
  ()

  A single path:
  >>> m = Chem.MolFromSmiles('CCCC')
  >>> c1 = long(Utils.GetAtomCode(m.GetAtomWithIdx(0),1))
  >>> c2 = long(Utils.GetAtomCode(m.GetAtomWithIdx(1),2))
  >>> c3 = long(Utils.GetAtomCode(m.GetAtomWithIdx(2),2))
  >>> c4 = long(Utils.GetAtomCode(m.GetAtomWithIdx(3),1))
  >>> t = c1 | (c2 << Utils.codeSize) | (c3 << (Utils.codeSize*2)) | (c4 << (Utils.codeSize*3))
  >>> GetTopologicalTorsionFingerprint(m)==(t,)
  1

  Two paths, both the same:
  >>> m = Chem.MolFromSmiles('CCCCC')
  >>> c4 = long(Utils.GetAtomCode(m.GetAtomWithIdx(3),1))
  >>> t = c1 | (c2 << Utils.codeSize) | (c3 << (Utils.codeSize*2)) | (c4 << (Utils.codeSize*3))
  >>> GetTopologicalTorsionFingerprint(m)==(t,t)
  1

  Real molecules:
  >>> m1 = Chem.MolFromSmiles('CCN(CCO)CC')
  >>> m2 = Chem.MolFromSmiles('CCN(CCO)CCO')
  >>> m3 = Chem.MolFromSmiles('OCCN(CCO)CCO')
  >>> fp1 = GetTopologicalTorsionFingerprint(m1)
  >>> fp2 = GetTopologicalTorsionFingerprint(m2)
  >>> fp3 = GetTopologicalTorsionFingerprint(m3)
  >>> Utils.DiceSimilarity(fp1,fp1)
  1.0
  >>> print '%.3f'%Utils.DiceSimilarity(fp1,fp2)
  0.667
  >>> print '%.3f'%Utils.DiceSimilarity(fp1,fp3)
  0.375
  >>> print '%.3f'%Utils.DiceSimilarity(fp2,fp3)
  0.706
  >>> print '%.3f'%Utils.DiceSimilarity(fp3,fp1)
  0.375

  """
  atomScores=[]
  for i in range(mol.GetNumAtoms()):
    atomScores.append(Utils.GetAtomCode(mol.GetAtomWithIdx(i)))

  paths = Chem.FindAllPathsOfLengthN(mol,targetSize,useBonds=0)
  nPaths = len(paths)
  res = [ScorePath(mol,x,targetSize,atomCodes=atomScores) for x in paths]
  #for x in paths:
  #  print [int(y) for y in x],ScorePath(mol,x,targetSize)

  res.sort()  
  return tuple(res)
  
def GetTopologicalTorsionFingerprintAsIntVect(mol,targetSize=4):
  """

  Real molecules:
  >>> m1 = Chem.MolFromSmiles('CCN(CCO)CC')
  >>> m2 = Chem.MolFromSmiles('CCN(CCO)CCO')
  >>> m3 = Chem.MolFromSmiles('OCCN(CCO)CCO')
  >>> fp1 = GetTopologicalTorsionFingerprintAsIntVect(m1)
  >>> fp2 = GetTopologicalTorsionFingerprintAsIntVect(m2)
  >>> fp3 = GetTopologicalTorsionFingerprintAsIntVect(m3)
  
  >>> from DataStructs import SparseIntVect
  >>> SparseIntVect.DiceSimilarity(fp1,fp1)
  1.0
  >>> print '%.3f'%SparseIntVect.DiceSimilarity(fp1,fp2)
  0.667
  >>> print '%.3f'%SparseIntVect.DiceSimilarity(fp1,fp3)
  0.375
  >>> print '%.3f'%SparseIntVect.DiceSimilarity(fp2,fp3)
  0.706
  >>> print '%.3f'%SparseIntVect.DiceSimilarity(fp3,fp1)
  0.375

  """
  from DataStructs import LongSparseIntVect as ResultType
  nBits = Utils.codeSize*targetSize
  sz = (1L<<nBits)-1
  res = ResultType(sz)
  ttv = GetTopologicalTorsionFingerprint(mol,targetSize=targetSize)
  res.UpdateFromSequence(ttv)
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
  
  

