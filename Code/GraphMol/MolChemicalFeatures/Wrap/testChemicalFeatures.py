import os
import unittest

from rdkit import Chem, RDConfig
from rdkit.Chem import ChemicalFeatures

try:
  from rdkit.Chem import AllChem
  haveAllChem = True
except ImportError:
  haveAllChem = False


def lstFeq(l1, l2, tol=1.e-4):
  if (len(list(l1)) != len(list(l2))):
    return 0
  for i in range(len(list(l1))):
    if not feq(l1[i], l2[i], tol):
      return 0
  return 1


def feq(v1, v2, tol2=1e-4):
  return abs(v1 - v2) <= tol2


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def testBasic(self):
    cfac = ChemicalFeatures.BuildFeatureFactory(
      os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'MolChemicalFeatures', 'test_data',
                   'featDef.txt'))
    self.assertTrue(cfac.GetNumFeatureDefs() == 2)

    fNames = cfac.GetFeatureFamilies()
    self.assertTrue(len(fNames) == 2)
    self.assertTrue(fNames[0] == 'HBondDonor')
    self.assertTrue(fNames[1] == 'HBondAcceptor')

    mol = Chem.MolFromSmiles(
      "COCN |(-1.22855,-0.651312,-0.097783;-0.706645,0.599392,0.182439;0.631075,0.659854,-0.17709;1.30412,-0.607935,0.0924339)|"
    )

    self.assertTrue(cfac.GetNumMolFeatures(mol) == 3)
    for i in range(cfac.GetNumMolFeatures(mol)):
      self.assertTrue(cfac.GetMolFeature(mol, i))
    # check that the recompute argument works:
    self.assertTrue(cfac.GetMolFeature(mol, 0))
    for i in range(cfac.GetNumMolFeatures(mol)):
      self.assertTrue(cfac.GetMolFeature(mol, i, "", False))
    self.assertRaises(IndexError, lambda: cfac.GetMolFeature(mol, 3))

    feats = cfac.GetFeaturesForMol(mol)
    self.assertTrue(len(feats) == 3)
    fTypes = ['HBondDonor', 'HBondAcceptor', 'HBondAcceptor']

    positions = [[1.3041, -0.6079, 0.0924], [-0.7066, 0.5994, 0.1824], [1.3041, -0.6079, 0.0924]]
    targetAids = [[3], [1], [3]]
    for i, feat in enumerate(feats):
      self.assertEqual(feat.GetFamily(), fTypes[i])
      pos = list(feat.GetPos())
      aids = list(feat.GetAtomIds())
      self.assertEqual(aids, targetAids[i])
      self.assertTrue(lstFeq(pos, positions[i]))
      nmol = feat.GetMol()
      self.assertEqual(Chem.MolToSmiles(nmol), "COCN")
      ncfac = feat.GetFactory()
      self.assertEqual(ncfac.GetNumFeatureDefs(), 2)
      self.assertEqual(feat.GetActiveConformer(), -1)

  def testIncludeOnly(self):
    cfac = ChemicalFeatures.BuildFeatureFactory(
      os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'MolChemicalFeatures', 'test_data',
                   'featDef.txt'))
    self.assertTrue(cfac.GetNumFeatureDefs() == 2)

    mol = Chem.MolFromSmiles(
      "COCN |(-1.22855,-0.651312,-0.097783;-0.706645,0.599392,0.182439;0.631075,0.659854,-0.17709;1.30412,-0.607935,0.0924339)|"
    )

    self.assertTrue(cfac.GetNumMolFeatures(mol, includeOnly="HBondAcceptor") == 2)
    self.assertTrue(cfac.GetNumMolFeatures(mol, includeOnly="HBondDonor") == 1)
    self.assertTrue(cfac.GetNumMolFeatures(mol, includeOnly="Bogus") == 0)

    self.assertRaises(IndexError, lambda: cfac.GetMolFeature(mol, 1, includeOnly="HBondDonor"))
    self.assertRaises(IndexError, lambda: cfac.GetMolFeature(mol, 2, includeOnly="HBondAcceptor"))
    f = cfac.GetMolFeature(mol, 0, includeOnly="HBondDonor")
    self.assertTrue(f.GetFamily() == 'HBondDonor')

    feats = cfac.GetFeaturesForMol(mol, includeOnly="HBondAcceptor")
    self.assertTrue(len(feats) == 2)
    feats = cfac.GetFeaturesForMol(mol, includeOnly="HBondDonor")
    self.assertTrue(len(feats) == 1)
    feats = cfac.GetFeaturesForMol(mol, includeOnly="Bogus")
    self.assertTrue(len(feats) == 0)

  def testStringParse(self):
    fdefBlock = \
"""DefineFeature HDonor1 [N,O;!H0]
    Family HBondDonor
    Weights 1.0
EndFeature
DefineFeature HAcceptor1 [N,O;H0]
    Family HBondAcceptor
    Weights 1.0
EndFeature
"""
    cfac = ChemicalFeatures.BuildFeatureFactoryFromString(fdefBlock)
    self.assertTrue(cfac.GetNumFeatureDefs() == 2)

  def testStringParse2(self):
    fdefBlock = \
"""DefineFeature HDonor1 [N,O;!H0]\r
    Family HBondDonor\r
    Weights 1.0\r
EndFeature\r
DefineFeature HAcceptor1 [N,O;H0]\r
    Family HBondAcceptor\r
    Weights 1.0\r
EndFeature\r
"""
    cfac = ChemicalFeatures.BuildFeatureFactoryFromString(fdefBlock)
    self.assertTrue(cfac.GetNumFeatureDefs() == 2)

  def testParseErrorHandling(self):
    fdefBlock = \
"""DefineFeature HDonor1 [N,O;!HQ]
    Family HBondDonor
    Weights 1.0
EndFeature
"""
    self.assertRaises(ValueError, lambda: ChemicalFeatures.BuildFeatureFactoryFromString(fdefBlock))
    fdefBlock = \
"""DefineFeature HDonor1 [N,O;!H0]
    Family HBondDonor
    Weights 1.0
"""
    self.assertRaises(ValueError, lambda: ChemicalFeatures.BuildFeatureFactoryFromString(fdefBlock))

    self.assertRaises(IOError, lambda: ChemicalFeatures.BuildFeatureFactory('noSuchFile.txt'))

  def testAtomMatch(self):
    fdefBlock = \
"""
DefineFeature HAcceptor1 [#7,#8]
    Family HBondAcceptor
    Weights 1.0
EndFeature
DefineFeature Arom1 a1aaaaa1
    Family Aromatic
    Weights 1.0,1.0,1.0,1.0,1.0,1.0
EndFeature
"""
    cfac = ChemicalFeatures.BuildFeatureFactoryFromString(fdefBlock)
    self.assertTrue(cfac.GetNumFeatureDefs() == 2)
    mol = Chem.MolFromSmiles('n1ccccc1')
    feats = cfac.GetFeaturesForMol(mol)
    self.assertTrue(len(feats) == 2)
    m = ChemicalFeatures.GetAtomMatch(feats)
    self.assertFalse(m)

    mol = Chem.MolFromSmiles('c1ccccc1N')
    feats = cfac.GetFeaturesForMol(mol)
    self.assertTrue(len(feats) == 2)
    m = ChemicalFeatures.GetAtomMatch(feats)
    self.assertTrue(len(m) == 2)

  def testIssue231(self):
    fdefs = """
DefineFeature HDonor1 [N,O;!H0]
  Family HBondDonor
  Weights 1.0
EndFeature
DefineFeature HAcceptor1 [N,O;H0]
  Family HBondAcceptor
  Weights 1.0
EndFeature
"""
    cfac = ChemicalFeatures.BuildFeatureFactoryFromString(fdefs)

    m = Chem.MolFromSmiles(
      'NCCC=O |(1.3336,-0.801443,-0.33565;0.753873,0.440597,0.0298162;-0.72858,0.533983,-0.182908;-1.49821,-0.46216,0.577961;-2.71702,-0.472209,0.474729)|'
    )
    feats = cfac.GetFeaturesForMol(m)
    for feat in feats:
      feat.GetPos()
    m = None
    for feat in feats:
      feat.GetPos()

  def testGithub2603(self):
    cfac = ChemicalFeatures.BuildFeatureFactory(
      os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef"))
    m = Chem.MolFromSmiles('OCc1ccccc1CN')
    feats = cfac.GetFeaturesForMol(m)
    self.assertEqual(feats[0].GetFamily(), 'Donor')
    cfac = None
    self.assertEqual(feats[0].GetFamily(), 'Donor')

  @unittest.skipIf(not haveAllChem, "AllChem not available in nanobind yet")
  def testGithub2530(self):
    cfac = ChemicalFeatures.BuildFeatureFactory(
      os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef"))
    m = Chem.MolFromSmiles(
      'COC1CCC1 |(2.38103,0.292767,0.259223;1.17409,-0.100802,0.882739;0.12366,0.138431,0.027286;-1.05585,0.888072,0.548495;-1.88253,-0.389559,0.484647;-0.708019,-1.10224,-0.133692)|'
    )
    feats = cfac.GetFeaturesForMol(m)
    feat_pos = feats[0].GetPos(-1)
    feat_pos_default = feats[0].GetPos()
    self.assertEqual(feat_pos[0], feat_pos_default[0])
    self.assertEqual(feat_pos[1], feat_pos_default[1])
    self.assertEqual(feat_pos[2], feat_pos_default[2])

    # another conformer
    m2 = Chem.MolFromSmiles(
      'COC1CCC1 |(2.38405,-0.455684,-0.194253;1.27804,-0.157448,0.590366;0.0931051,-0.0947927,-0.0913265;-0.719418,1.1751,0.106878;-1.93955,0.274584,0.118445;-1.03964,-0.801813,0.63462)|'
    )
    m2.AddConformer(m.GetConformer(), assignId=True)

    feats_0 = cfac.GetFeaturesForMol(m2, confId=0)
    feats_1 = cfac.GetFeaturesForMol(m2, confId=1)
    # raise ValueError(list(feats_0))
    self.assertNotEqual(feats_1[0], feats_0[0])


if __name__ == '__main__':
  unittest.main()
