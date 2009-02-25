# $Id$
#
# Copyright (C) 2003-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" functionality for finding pharmacophore matches in molecules


  See Docs/Chem/Pharm2D.triangles.jpg for an illustration of the way
  pharmacophores are broken into triangles and labelled.

  See Docs/Chem/Pharm2D.signatures.jpg for an illustration of bit
  numbering

"""
from rdkit import Chem
from rdkit.Chem.Pharm2D import Utils

import types
import exceptions
class MatchError(exceptions.Exception):
  pass

_verbose = 0
def GetAtomsMatchingBit(sig,bitIdx,mol,dMat=None,justOne=0,matchingAtoms=None):
  """ Returns a list of lists of atom indices for a bit

    **Arguments**

      - sig: a signature

      - bitIdx: the bit to be queried

      - mol: the molecule to be examined

      - dMat: (optional) the distance matrix of the molecule

      - justOne: (optional) if this is nonzero, only the first match
        will be returned.

      - matchingAtoms: (optional) if this is nonzero, it should
        contain a sequence of sequences with the indices of atoms in
        the molecule which match each of the patterns used by the
        signature.

    **Returns**

      a list of tuples with the matching atoms
  """
  assert sig.GetShortestPathsOnly(),'not implemented for non-shortest path signatures'
  nPts,pattCombo,scaffold = sig.GetBitInfo(bitIdx)
  if _verbose:
    print 'info:',nPts
    print '\t',pattCombo
    print '\t',scaffold
  
  # find the atoms that match each pattern
  choices = []
  for pattIdx in pattCombo:
    tmp = None
    if matchingAtoms is not None:
      tmp = matchingAtoms[pattIdx]
    else:
      patt = sig.GetPattern(pattIdx)
      tmp = mol.GetSubstructMatches(patt)
    if tmp:
      choices.append(tmp)
    else:
      # one of the patterns didn't find a match, we
      #  can return now
      if _verbose: print 'no match found for pattern:',pattIdx
      return []
 
  if _verbose:
    print 'choices:'
    print choices

  if dMat is None:
    useBO = sig.GetIncludeBondOrder()
    dMat = Chem.GetDistanceMatrix(mol,useBO)

  matches = []
  distsToCheck = Utils.nPointDistDict[nPts]

  protoPharmacophores = Utils.GetAllCombinations(choices,noDups=1)

  res = []
  for protoPharm in protoPharmacophores:
    if _verbose: print 'protoPharm:',protoPharm
    for i in range(len(distsToCheck)):
      dLow,dHigh = sig.GetBin(scaffold[i])
      a1,a2 = distsToCheck[i]
      #
      # FIX: this is making all kinds of assumptions about
      #  things being single-atom matches (or at least that
      #  only the first atom matters
      #
      idx1,idx2 = protoPharm[a1][0],protoPharm[a2][0]
      dist = dMat[idx1,idx2]
      if _verbose: print '\t dist: %d->%d = %d (%d,%d)'%(idx1,idx2,dist,dLow,dHigh)
      if dist < dLow or dist >= dHigh:
        break
    else:
      if _verbose: print 'Found one'
      # we found it
      protoPharm.sort()
      protoPharm = tuple(protoPharm)
      if protoPharm not in res:
        res.append(protoPharm)
        if justOne: break
  return res

if __name__ == '__main__':
  from rdkit import Chem
  from rdkit.Chem.Pharm2D import SigFactory,Generate
      
  factory = SigFactory.SigFactory()
  factory.SetBins([(1,2),(2,5),(5,8)])
  factory.SetPatternsFromSmarts(['O','N'])
  factory.SetMinCount(2)
  factory.SetMaxCount(3)
  sig = factory.GetSignature()
  
  mol = Chem.MolFromSmiles('OCC(=O)CCCN')
  Generate.Gen2DFingerprint(mol,sig)
  print 'onbits:',list(sig.GetOnBits())

  _verbose=0
  for bit in sig.GetOnBits():
    ats = GetAtomsMatchingBit(sig,bit,mol)
    print '\tBit %d: '%(bit),ats

    
  print '--------------------------'
  sig = factory.GetSignature()
  sig.SetIncludeBondOrder(1)
  Generate.Gen2DFingerprint(mol,sig)
  print 'onbits:',list(sig.GetOnBits())

  for bit in sig.GetOnBits():
    ats = GetAtomsMatchingBit(sig,bit,mol)
    print '\tBit %d: '%(bit),ats

  
  
