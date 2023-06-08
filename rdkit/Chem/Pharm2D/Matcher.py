#
# Copyright (C) 2003-2008 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" functionality for finding pharmacophore matches in molecules


  See Docs/Chem/Pharm2D.triangles.jpg for an illustration of the way
  pharmacophores are broken into triangles and labelled.

  See Docs/Chem/Pharm2D.signatures.jpg for an illustration of bit
  numbering

"""
from rdkit import Chem
from rdkit.Chem.Pharm2D import Utils


class MatchError(Exception):
  pass


_verbose = 0


def GetAtomsMatchingBit(sigFactory, bitIdx, mol, dMat=None, justOne=0, matchingAtoms=None):
  """ Returns a list of lists of atom indices for a bit

    **Arguments**

      - sigFactory: a SigFactory

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
  assert sigFactory.shortestPathsOnly, 'not implemented for non-shortest path signatures'
  nPts, featCombo, scaffold = sigFactory.GetBitInfo(bitIdx)
  if _verbose:
    print('info:', nPts)
    print('\t', featCombo)
    print('\t', scaffold)

  if matchingAtoms is None:
    matchingAtoms = sigFactory.GetMolFeats(mol)

  # find the atoms that match each features
  # fams = sigFactory.GetFeatFamilies()
  choices = []
  for featIdx in featCombo:
    tmp = matchingAtoms[featIdx]
    if tmp:
      choices.append(tmp)
    else:
      # one of the patterns didn't find a match, we
      #  can return now
      if _verbose:
        print('no match found for feature:', featIdx)
      return []

  if _verbose:
    print('choices:')
    print(choices)

  if dMat is None:
    dMat = Chem.GetDistanceMatrix(mol, sigFactory.includeBondOrder)

  distsToCheck = Utils.nPointDistDict[nPts]

  protoPharmacophores = Utils.GetAllCombinations(choices, noDups=1)

  res = []
  for protoPharm in protoPharmacophores:
    if _verbose:
      print('protoPharm:', protoPharm)
    for i in range(len(distsToCheck)):
      dLow, dHigh = sigFactory.GetBins()[scaffold[i]]
      a1, a2 = distsToCheck[i]
      #
      # FIX: this is making all kinds of assumptions about
      #  things being single-atom matches (or at least that
      #  only the first atom matters
      #
      idx1, idx2 = protoPharm[a1][0], protoPharm[a2][0]
      dist = dMat[idx1, idx2]
      if _verbose:
        print(f'\t dist: {idx1}->{idx2} = {dist} ({dLow}, {dHigh})')
      if dist < dLow or dist >= dHigh:
        break
    else:
      if _verbose:
        print('Found one')
      # we found it
      protoPharm.sort()
      protoPharm = tuple(protoPharm)
      if protoPharm not in res:
        res.append(protoPharm)
        if justOne:
          break
  return res


def _exampleCode():
  import os

  from rdkit import RDConfig
  from rdkit.Chem import ChemicalFeatures
  from rdkit.Chem.Pharm2D import Generate, SigFactory

  fdefFile = os.path.join(RDConfig.RDCodeDir, 'Chem', 'Pharm2D', 'test_data', 'BaseFeatures.fdef')
  featFactory = ChemicalFeatures.BuildFeatureFactory(fdefFile)
  factory = SigFactory.SigFactory(featFactory)
  factory.SetBins([(1, 2), (2, 5), (5, 8)])
  factory.Init()

  mol = Chem.MolFromSmiles('OCC(=O)CCCN')
  sig = Generate.Gen2DFingerprint(mol, factory)
  print('onbits:', list(sig.GetOnBits()))

  _verbose = 0
  for bit in sig.GetOnBits():
    atoms = GetAtomsMatchingBit(factory, bit, mol)
    print(f'\tBit {bit}: ', atoms)

  print('finished')


if __name__ == '__main__':  # pragma: nocover
  _exampleCode()
