import os
import unittest

from rdkit import Chem, RDConfig
from rdkit.Chem import ChemicalFeatures


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

    mol = Chem.MolFromSmiles('COCN |(-1.22855,-0.651312,-0.097783;-0.706645,0.599392,0.182439;'
                             '0.631075,0.659854,-0.17709;1.30412,-0.607935,0.0924339)|')

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

    mol = Chem.MolFromSmiles('COCN |(-1.30552,-0.589627,-0.120346;-0.642675,0.563645,0.297155;'
                             '0.602117,0.598768,-0.290102;1.34608,-0.572786,0.113293)|')

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

    m = Chem.MolFromSmiles('NCCC=O |(1.59083,0.424048,0.598431;1.15895,-0.185728,-0.633195;'
                           '-0.109831,-0.972839,-0.525873;-1.2816,-0.214521,-0.0919036;'
                           '-1.35836,0.94904,0.187585)|')
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

  def testGithub2530(self):
    cfac = ChemicalFeatures.BuildFeatureFactory(
      os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef"))
    m = Chem.MolFromSmiles('COC1CCC1 |(1.58784,0.882463,-0.310971;'
                           '1.37181,-0.430838,0.0331765;0.0549765,-0.775729,0.263463;'
                           '-0.921924,-0.600844,-0.852632;'
                           '-1.35129,0.714273,-0.218004;-0.741412,0.210675,1.08497)|')
    feats = cfac.GetFeaturesForMol(m)
    feat_pos = feats[0].GetPos(-1)
    feat_pos_default = feats[0].GetPos()
    self.assertEqual(feat_pos[0], feat_pos_default[0])
    self.assertEqual(feat_pos[1], feat_pos_default[1])
    self.assertEqual(feat_pos[2], feat_pos_default[2])

    # Conformers generation:
    m2 = Chem.MultiConfMolFromSDF(
      os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'MolChemicalFeatures', 'test_data',
                   'testGH2530.sdf'), removeHs=False)

    feats_0 = cfac.GetFeaturesForMol(m2, confId=-1)
    feats_5 = cfac.GetFeaturesForMol(m2, confId=5)
    self.assertNotEqual(feats_5[0], feats_0[0])
    self.assertNotEqual(feats_5[1], feats_0[1])
    self.assertNotEqual(feats_5[2], feats_0[2])


if __name__ == '__main__':
  unittest.main()
