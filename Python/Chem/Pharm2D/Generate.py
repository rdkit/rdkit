# $Id: Generate.py 5007 2006-02-22 15:14:41Z glandrum $
#
# Copyright (C) 2002-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" generation of 2D pharmacophores

**Notes**

  - The terminology for this gets a bit rocky, so here's a glossary of what
    terms used here mean:

      1) *N-point pharmacophore* a combination of N features along with
         distances betwen them.

      2) *N-point proto-pharmacophore*: a combination of N feature
         definitions without distances.  Each N-point
         proto-pharmacophore defines a manifold of potential N-point
         pharmacophores.

      3) *N-point scaffold*: a collection of the distances defining
         an N-point pharmacophore without feature identities.


  See Docs/Chem/Pharm2D.triangles.jpg for an illustration of the way
  pharmacophores are broken into triangles and labelled.

  See Docs/Chem/Pharm2D.signatures.jpg for an illustration of bit
  numbering

"""
from Chem.Pharm2D import Utils,Signature,SigFactory

_verbose = 0

def _ShortestPathsMatch(match,featureSet,sig,dMat):
  """  Internal use only

  """
  if _verbose:
    print 'match:',match
  nPts = len(match)
  distsToCheck =  Utils.nPointDistDict[nPts]
  nDists = len(distsToCheck)
  dist = [0]*nDists
  #minD,maxD = sig.GetMinDist(),sig.GetMaxDist()
  minD,maxD = sig.bins[0][0],sig.bins[-1][1]
  
  for i in range(nDists):
    pt0,pt1 = distsToCheck[i]
    idx1,idx2 = match[pt0][0],match[pt1][0]
    # FIX: treats all matches as single-atom
    # FIX: is it ok that this is an int?
    d = int(dMat[idx1,idx2])
    # do a quick distance filter
    if d == 0 or d <= minD or d > maxD:
      return
    dist[i] = d
  if _verbose:
    print '\t',dist,sig.GetMinDist(),sig.GetMaxDist()
  sig.SetBit(featureSet,dist,checkPatts=0)

  

def Gen2DFingerprint(mol,sig,perms=None,dMat=None):
  """ generates a 2D fingerprint for a molecule using the
   parameters in _sig_

   **Arguments**

     - mol: the molecule for which the signature should be generated

     - sig: the signature object which will be filled.

     - perms: (optional) a sequence of permutation indices limiting which
       pharmacophore combinations are allowed

     - dMat: (optional) the distance matrix to be used

   **Notes**  

     - no preprocessing is carried out for _sig_. It *must* be pre-initialized.
       Any bits which are already set will not be changed by this operation.

   
  """
  if isinstance(sig,SigFactory.SigFactory):
    sig = sig.GetSignature()
  nPatts = sig.GetNumPatterns()
  minCount = sig.GetMinCount()
  maxCount = sig.GetMaxCount()

  # generate the molecule's distance matrix, if required
  if dMat is None:
    import Chem
    useBO = sig.GetIncludeBondOrder()
    dMat = Chem.GetDistanceMatrix(mol,useBO)

  # generate the permutations, if required
  if perms is None:
    perms = []
    for count in range(minCount,maxCount+1):
      perms += Utils.GetIndexCombinations(nPatts,count)

  # generate the matches:
  pattMatches = [None]*nPatts
  for i in range(nPatts):
    patt = sig.GetPattern(i)
    pattMatches[i] = mol.GetSubstructMatches(patt)

  for perm in perms:
    # the permutation is a combination of feature indices
    #   defining the feature set for a proto-pharmacophore

    # Get a set of matches at each index of
    #  the proto-pharmacophore.
    matchPerms = [pattMatches[x] for x in perm]
    if _verbose:
      print '\n->Perm: %s'%(str(perm))
      print '    matchPerms: %s'%(str(matchPerms))
    
    # Get all unique combinations of those possible matches:
    matchesToMap = Utils.UniquifyCombinations(Utils.GetAllCombinations(matchPerms))
    
    for match in matchesToMap:
      if sig.GetShortestPathsOnly():
        _ShortestPathsMatch(match,perm,sig,dMat)

  return sig
if __name__ == '__main__':
  import Chem
  def test1():
    sig = Signature.Pharm2DSig()
    sig.SetPatternsFromSmarts(['O','N'])
    sig.SetBins([(1,2),(2,4),(4,6)])
    sig.SetMinCount(2)
    sig.SetMaxCount(3)
    sig.Init()

    mol = Chem.MolFromSmiles('C1OCOCC1O')
    Gen2DFingerprint(mol,sig)

    print list(sig.GetOnBits())
    print sig._size

  def test2():
    sig = Signature.Pharm2DSig()
    sig.SetPatternsFromSmarts(['[OD1]','[OD2]','[CD3]',])
    sig.SetBins([(1,3),(3,8)])
    sig.SetMinCount(3)
    sig.SetMaxCount(3)
    sig.Init()
    mol = Chem.MolFromSmiles('OCCC1COCCO1')
    dMat= Chem.GetDistanceMatrix(mol)

    #Signature._verbose = 1
    Gen2DFingerprint(mol,sig)

    #print sig._binCombos
    print list(sig.GetOnBits())
    print sig._size
  Signature._verbose=1
  _verbose = 1
  test1()
