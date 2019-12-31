#
# Copyright (C) 2002-2008 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" generation of 2D pharmacophores

**Notes**

  - The terminology for this gets a bit rocky, so here's a glossary of what
    terms used here mean:

      1) *N-point pharmacophore* a combination of N features along with
         distances between them.

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


from rdkit.Chem.Pharm2D import Utils, SigFactory
from rdkit.RDLogger import logger

logger = logger()

_verbose = 0


def _ShortestPathsMatch(match, featureSet, sig, dMat, sigFactory):
  """  Internal use only

  """
  if _verbose:
    print('match:', match)
  nPts = len(match)
  distsToCheck = Utils.nPointDistDict[nPts]
  nDists = len(distsToCheck)
  dist = [0] * nDists
  bins = sigFactory.GetBins()
  minD, maxD = bins[0][0], bins[-1][1]

  for i in range(nDists):
    pt0, pt1 = distsToCheck[i]
    minSeen = maxD
    for idx1 in match[pt0]:
      for idx2 in match[pt1]:
        minSeen = min(minSeen, dMat[idx1, idx2])
        if minSeen == 0 or minSeen < minD:
          return
    # FIX: this won't be an int if we're using the bond order.
    d = int(minSeen)
    # do a quick distance filter
    if d == 0 or d < minD or d >= maxD:
      return
    dist[i] = d

  idx = sigFactory.GetBitIdx(featureSet, dist, sortIndices=False)
  if _verbose:
    print('\t', dist, minD, maxD, idx)

  if sigFactory.useCounts:
    sig[idx] = sig[idx] + 1
  else:
    sig.SetBit(idx)
  return idx


def Gen2DFingerprint(mol, sigFactory, perms=None, dMat=None, bitInfo=None):
  """ generates a 2D fingerprint for a molecule using the
   parameters in _sig_

   **Arguments**

     - mol: the molecule for which the signature should be generated

     - sigFactory : the SigFactory object with signature parameters
       NOTE: no preprocessing is carried out for _sigFactory_.
             It *must* be pre-initialized.

     - perms: (optional) a sequence of permutation indices limiting which
       pharmacophore combinations are allowed

     - dMat: (optional) the distance matrix to be used

     - bitInfo: (optional) used to return the atoms involved in the bits

  """
  if not isinstance(sigFactory, SigFactory.SigFactory):
    raise ValueError('bad factory')
  featFamilies = sigFactory.GetFeatFamilies()
  if _verbose:
    print('* feat famillies:', featFamilies)
  nFeats = len(featFamilies)
  minCount = sigFactory.minPointCount
  maxCount = sigFactory.maxPointCount
  if maxCount > 3:
    logger.warning(' Pharmacophores with more than 3 points are not currently supported.\n' +
                   'Setting maxCount to 3.')
    maxCount = 3

  # generate the molecule's distance matrix, if required
  if dMat is None:
    from rdkit import Chem
    useBO = sigFactory.includeBondOrder
    dMat = Chem.GetDistanceMatrix(mol, useBO)

  # generate the permutations, if required
  if perms is None:
    perms = []
    for count in range(minCount, maxCount + 1):
      perms += Utils.GetIndexCombinations(nFeats, count)

  # generate the matches:
  featMatches = sigFactory.GetMolFeats(mol)
  if _verbose:
    print('  featMatches:', featMatches)

  sig = sigFactory.GetSignature()
  for perm in perms:
    # the permutation is a combination of feature indices
    #   defining the feature set for a proto-pharmacophore
    featClasses = [0]
    for i in range(1, len(perm)):
      if perm[i] == perm[i - 1]:
        featClasses.append(featClasses[-1])
      else:
        featClasses.append(featClasses[-1] + 1)

    # Get a set of matches at each index of
    #  the proto-pharmacophore.
    matchPerms = [featMatches[x] for x in perm]
    if _verbose:
      print('\n->Perm: %s' % (str(perm)))
      print('    matchPerms: %s' % (str(matchPerms)))

    # Get all unique combinations of those possible matches:
    matchesToMap = Utils.GetUniqueCombinations(matchPerms, featClasses)
    for i, entry in enumerate(matchesToMap):
      entry = [x[1] for x in entry]
      matchesToMap[i] = entry
    if _verbose:
      print('    mtM:', matchesToMap)

    for match in matchesToMap:
      if sigFactory.shortestPathsOnly:
        idx = _ShortestPathsMatch(match, perm, sig, dMat, sigFactory)
        if idx is not None and bitInfo is not None:
          l = bitInfo.get(idx, [])
          l.append(match)
          bitInfo[idx] = l
  return sig
