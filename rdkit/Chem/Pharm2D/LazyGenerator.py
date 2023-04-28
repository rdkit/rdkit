#
# Copyright (C) 2003-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" lazy generator of 2D pharmacophore signature data

"""

from rdkit.Chem.Pharm2D import Matcher, SigFactory

raise NotImplementedError('not finished yet')


class Generator(object):
  """

  Important attributes:

   - mol: the molecules whose signature is being worked with

   - sigFactory : the SigFactory object with signature parameters
            NOTE: no preprocessing is carried out for _sigFactory_.
                  It *must* be pre-initialized.

   **Notes**

     -
  """

  def __init__(self, sigFactory, mol, dMat=None, bitCache=True):
    """ constructor

      **Arguments**

       - sigFactory: a signature factory, see class docs

       - mol: a molecule, see class docs

       - dMat: (optional) a distance matrix for the molecule.  If this
         is not provided, one will be calculated

       - bitCache: (optional) if nonzero, a local cache of which bits
         have been queried will be maintained.  Otherwise things must
         be recalculate each time a bit is queried.

    """
    if not isinstance(sigFactory, SigFactory.SigFactory):
      raise ValueError('bad factory')

    self.sigFactory = sigFactory
    self.mol = mol

    if dMat is None:
      useBO = sigFactory.includeBondOrder
      dMat = Chem.GetDistanceMatrix(mol, useBO)

    self.dMat = dMat

    if bitCache:
      self.bits = {}
    else:
      self.bits = None

    featFamilies = [
      fam for fam in sigFactory.featFactory.GetFeatureFamilies() if fam not in sigFactory.skipFeats
    ]
    nFeats = len(featFamilies)
    featMatches = {}
    for fam in featFamilies:
      featMatches[fam] = []
    feats = sigFactory.featFactory.GetFeaturesForMol(mol)
    for feat in feats:
      if feat.GetFamily() not in sigFactory.skipFeats:
        featMatches[feat.GetFamily()].append(feat.GetAtomIds())
    featMatches = [None] * nFeats
    for i in range(nFeats):
      featMatches[i] = sigFactory.featFactory.GetMolFeature()
    self.pattMatches = pattMatches

  def GetBit(self, idx):
    """ returns a bool indicating whether or not the bit is set

    """
    if idx < 0 or idx >= self.sig.GetSize():
      raise IndexError('Index %d invalid' % (idx))
    if self.bits is not None and idx in self.bits:
      return self.bits[idx]

    tmp = Matcher.GetAtomsMatchingBit(self.sig, idx, self.mol, dMat=self.dMat, justOne=1,
                                      matchingAtoms=self.pattMatches)
    if not tmp or len(tmp) == 0:
      res = 0
    else:
      res = 1

    if self.bits is not None:
      self.bits[idx] = res
    return res

  def __len__(self):
    """ allows class to support len()

    """
    return self.sig.GetSize()

  def __getitem__(self, itm):
    """ allows class to support random access.
      Calls self.GetBit()

    """
    return self.GetBit(itm)


if __name__ == '__main__':
  import random
  import time

  from rdkit import Chem, RDConfig
  from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D

  factory = Gobbi_Pharm2D.factory
  nToDo = 100
  inD = open(RDConfig.RDDataDir + "/NCI/first_5K.smi", 'r').readlines()[:nToDo]
  mols = [None] * len(inD)
  for i in range(len(inD)):
    smi = inD[i].split('\t')[0]
    smi.strip()
    mols[i] = Chem.MolFromSmiles(smi)

  sig = factory.GetSignature()

  nBits = 300
  random.seed(23)
  bits = [random.randint(0, sig.GetSize() - 1) for x in range(nBits)]

  print('Using the Lazy Generator')
  t1 = time.time()
  for i in range(len(mols)):
    if not i % 10:
      print('done mol %d of %d' % (i, len(mols)))
    gen = Generator(factory, mols[i])
    for bit in bits:
      v = gen[bit]
  t2 = time.time()
  print('\tthat took %4.2f seconds' % (t2 - t1))

  print('Generating and checking signatures')
  t1 = time.time()
  for i in range(len(mols)):
    if not i % 10:
      print('done mol %d of %d' % (i, len(mols)))
    sig = Generate.Gen2DFingerprint(mols[i], factory)
    for bit in bits:
      v = sig[bit]
  t2 = time.time()
  print('\tthat took %4.2f seconds' % (t2 - t1))
