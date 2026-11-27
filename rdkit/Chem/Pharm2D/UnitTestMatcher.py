#
#  Copyright (C) 2002-2008  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import os.path
import unittest
from io import StringIO

from rdkit import Chem, RDConfig
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D, Matcher, SigFactory
from rdkit.TestRunner import redirect_stdout


class TestCase(unittest.TestCase):

  def setUp(self):
    fdefFile = os.path.join(RDConfig.RDCodeDir, 'Chem', 'Pharm2D', 'test_data', 'BaseFeatures.fdef')
    featFactory = ChemicalFeatures.BuildFeatureFactory(fdefFile)
    self.factory = SigFactory.SigFactory(featFactory, minPointCount=2, maxPointCount=3)
    self.factory.SetBins([(0, 2), (2, 5), (5, 8)])
    self.factory.Init()

  def test1_simple(self):
    mol = Chem.MolFromSmiles('OCC(=O)CCCN')
    self.factory.skipFeats = ['Donor']
    self.factory.Init()
    self.assertEqual(self.factory.GetSigSize(), 510)
    Generate._verbose = False
    sig = Generate.Gen2DFingerprint(mol, self.factory)
    Generate._verbose = False
    tgt = (1, 2, 11, 52, 117)
    onBits = sig.GetOnBits()
    self.assertEqual(tuple(onBits), tgt)
    self.assertEqual(len(onBits), len(tgt))

    bitMatches = (
      [((0, ), (3, ))],
      [((0, ), (7, )), ((3, ), (7, ))],
      [((0, ), (3, ), (7, ))],
    )
    for i, bit in enumerate(onBits):
      matches = Matcher.GetAtomsMatchingBit(self.factory, bit, mol)
      # print bit,matches
      # tgt = bitMatches[i]
      # self.assertEqual(matches,tgt)

  def test2Bug28(self):
    smi = r'Cc([s]1)nnc1SCC(\CS2)=C(/C([O-])=O)N3C(=O)[C@H]([C@@H]23)NC(=O)C[n]4cnnn4'
    mol = Chem.MolFromSmiles(smi)
    factory = Gobbi_Pharm2D.factory
    factory.SetBins([(2, 3), (3, 4), (4, 5), (5, 8), (8, 100)])
    sig = Generate.Gen2DFingerprint(mol, factory)
    onBits = sig.GetOnBits()
    for bit in onBits:
      atoms = Matcher.GetAtomsMatchingBit(factory, bit, mol, justOne=1)
      self.assertTrue(len(atoms))

  def testMultiAtomFeatureDistance(self):
    mol = Chem.MolFromSmiles('c1ccc(F)cc1')
    families = self.factory.GetFeatFamilies()
    pair = [families.index('Acceptor'), families.index('Aromatic')]
    near = self.factory.GetBitIdx(pair, [1])
    far = self.factory.GetBitIdx(pair, [4])
    bitInfo = {}
    fp = Generate.Gen2DFingerprint(mol, self.factory, bitInfo=bitInfo)
    self.assertTrue(fp[near])
    self.assertFalse(fp[far])
    expected = [tuple(sorted(bitInfo[near][0]))]
    self.assertEqual(Matcher.GetAtomsMatchingBit(self.factory, near, mol), expected)
    self.assertEqual(Matcher.GetAtomsMatchingBit(self.factory, far, mol), [])

  def testOverlappingFeatures(self):
    mol = Chem.MolFromSmiles('c1ccc2ccccc2c1')
    aromatic = self.factory.GetFeatFamilies().index('Aromatic')
    bit = self.factory.GetBitIdx([aromatic, aromatic], [1])
    self.assertFalse(Generate.Gen2DFingerprint(mol, self.factory)[bit])
    self.assertEqual(Matcher.GetAtomsMatchingBit(self.factory, bit, mol), [])

  def testBondOrderWeightedDistance(self):
    mol = Chem.MolFromSmiles('N#N')
    acceptor = self.factory.GetFeatFamilies().index('Acceptor')
    bit = self.factory.GetBitIdx([acceptor, acceptor], [1])
    self.assertTrue(Generate.Gen2DFingerprint(mol, self.factory)[bit])
    self.assertEqual(Matcher.GetAtomsMatchingBit(self.factory, bit, mol), [((0,), (1,))])
    self.factory.includeBondOrder = True
    self.assertFalse(Generate.Gen2DFingerprint(mol, self.factory)[bit])
    self.assertEqual(Matcher.GetAtomsMatchingBit(self.factory, bit, mol), [])

  def testRepeatedFeatureTriangle(self):
    mol = Chem.MolFromSmiles('Fc1ccc(CO)cc1')
    families = self.factory.GetFeatFamilies()
    featCombo = (families.index('Acceptor'), families.index('Acceptor'),
                 families.index('Aromatic'))
    bitInfo = {}
    fp = Generate.Gen2DFingerprint(mol, self.factory, bitInfo=bitInfo)
    self.assertEqual(sum(self.factory.GetBitInfo(i)[1] == featCombo for i in fp.GetOnBits()), 1)
    # Exchanging the two acceptors must not produce additional bit matches.
    for bit in range(self.factory.GetSigSize()):
      if self.factory.GetBitInfo(bit)[1] == featCombo:
        expected = [tuple(sorted(match)) for match in bitInfo.get(bit, [])]
        self.assertEqual(Matcher.GetAtomsMatchingBit(self.factory, bit, mol), expected)

  def testFractionalUpperBound(self):
    self.factory.SetBins([(0, 1.5), (1.5, 3.5)])
    self.factory.Init()
    mol = Chem.MolFromSmiles('NCCN')
    acceptor = self.factory.GetFeatFamilies().index('Acceptor')
    bit = self.factory.GetBitIdx([acceptor, acceptor], [3])
    dMat = Chem.GetDistanceMatrix(mol).copy()
    for distance, matches in [(3.25, [((0,), (3,))]), (3.5, []), (7., [])]:
      with self.subTest(distance=distance):
        dMat[0, 3] = dMat[3, 0] = distance
        fp = Generate.Gen2DFingerprint(mol, self.factory, perms=[[acceptor, acceptor]], dMat=dMat)
        self.assertEqual(bool(fp[bit]), bool(matches))
        self.assertEqual(Matcher.GetAtomsMatchingBit(self.factory, bit, mol, dMat=dMat), matches)

  def test3Roundtrip(self):
    # longer-running Bug 28 test
    nToDo = 20
    with open(os.path.join(RDConfig.RDDataDir, 'NCI', 'first_5K.smi'), 'r') as inF:
      inD = inF.readlines()[:nToDo]
    factory = Gobbi_Pharm2D.factory
    factory.SetBins([(2, 3), (3, 4), (4, 5), (5, 8), (8, 100)])
    for line in inD:
      smi = line.split('\t')[0]
      mol = Chem.MolFromSmiles(smi)
      sig = Generate.Gen2DFingerprint(mol, factory)
      onBits = sig.GetOnBits()
      for bit in onBits:
        atoms = Matcher.GetAtomsMatchingBit(factory, bit, mol, justOne=1)
        assert len(atoms), f'bit {bit} failed to match for smi {smi}'

  def test_exampleCode(self):
    # We make sure that the example code runs
    f = StringIO()
    with redirect_stdout(f):
      Matcher._exampleCode()
    self.assertIn('finished', f.getvalue())


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
