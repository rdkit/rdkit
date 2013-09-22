from rdkit import DataStructs
from rdkit import RDConfig
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures,rdDistGeom
from rdkit import Geometry
import unittest,os


def lstFeq(l1, l2, tol=1.e-4):
  if (len(list(l1)) != len(list(l2))):
    return 0
  for i in range(len(list(l1))):
    if not feq(l1[i], l2[i], tol):
      return 0
  return 1

def feq(v1,v2,tol2=1e-4):
  return abs(v1-v2)<=tol2

class TestCase(unittest.TestCase):
  def setUp(self):
    pass

  def testBasic(self):
    cfac = ChemicalFeatures.BuildFeatureFactory(os.path.join(RDConfig.RDBaseDir,
                                                             'Code','GraphMol','MolChemicalFeatures','test_data','featDef.txt'))
    self.failUnless(cfac.GetNumFeatureDefs() == 2)

    fNames = cfac.GetFeatureFamilies()
    self.failUnless(len(fNames) == 2)
    self.failUnless(fNames[0] == 'HBondDonor')
    self.failUnless(fNames[1] == 'HBondAcceptor')

    mol = Chem.MolFromSmiles("COCN")
    rdDistGeom.EmbedMolecule(mol, 30, 100)

    self.failUnless(cfac.GetNumMolFeatures(mol) == 3)
    for i in range(cfac.GetNumMolFeatures(mol)):
      self.failUnless(cfac.GetMolFeature(mol,i))
    # check that the recompute argument works:
    self.failUnless(cfac.GetMolFeature(mol,0))
    for i in range(cfac.GetNumMolFeatures(mol)):
      self.failUnless(cfac.GetMolFeature(mol,i,"",False))
    self.failUnlessRaises(IndexError,lambda : cfac.GetMolFeature(mol,3))


    feats = cfac.GetFeaturesForMol(mol)
    self.failUnless(len(feats) == 3)
    fTypes = ['HBondDonor', 'HBondAcceptor', 'HBondAcceptor']

    positions = [[1.3041, -0.6079, 0.0924],
                 [-0.7066, 0.5994, 0.1824],
                 [1.3041, -0.6079, 0.0924]]
    targetAids = [[3], [1], [3]]
    i = 0
    for feat in feats:
      self.failUnless(feat.GetFamily() == fTypes[i])
      pos = list(feat.GetPos())
      aids = list(feat.GetAtomIds())
      self.failUnless(aids == targetAids[i])
      self.failUnless(lstFeq(pos, positions[i]))
      nmol = feat.GetMol()
      self.failUnless(Chem.MolToSmiles(nmol) == "COCN")
      ncfac = feat.GetFactory()
      self.failUnless(ncfac.GetNumFeatureDefs() == 2)
      i += 1


  def testIncludeOnly(self):
    cfac = ChemicalFeatures.BuildFeatureFactory(os.path.join(RDConfig.RDBaseDir,
                                                             'Code','GraphMol','MolChemicalFeatures','test_data','featDef.txt'))
    self.failUnless(cfac.GetNumFeatureDefs() == 2)

    mol = Chem.MolFromSmiles("COCN")
    rdDistGeom.EmbedMolecule(mol)

    self.failUnless(cfac.GetNumMolFeatures(mol,includeOnly="HBondAcceptor") == 2)
    self.failUnless(cfac.GetNumMolFeatures(mol,includeOnly="HBondDonor") == 1)
    self.failUnless(cfac.GetNumMolFeatures(mol,includeOnly="Bogus") == 0)

    self.failUnlessRaises(IndexError,lambda : cfac.GetMolFeature(mol,1,includeOnly="HBondDonor"))
    self.failUnlessRaises(IndexError,lambda : cfac.GetMolFeature(mol,2,includeOnly="HBondAcceptor"))
    f = cfac.GetMolFeature(mol,0,includeOnly="HBondDonor")
    self.failUnless(f.GetFamily()=='HBondDonor')

    feats = cfac.GetFeaturesForMol(mol,includeOnly="HBondAcceptor")
    self.failUnless(len(feats) == 2)
    feats = cfac.GetFeaturesForMol(mol,includeOnly="HBondDonor")
    self.failUnless(len(feats) == 1)
    feats = cfac.GetFeaturesForMol(mol,includeOnly="Bogus")
    self.failUnless(len(feats) == 0)

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
    self.failUnless(cfac.GetNumFeatureDefs() == 2)


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
    self.failUnless(cfac.GetNumFeatureDefs() == 2)

  def testParseErrorHandling(self):
    fdefBlock = \
"""DefineFeature HDonor1 [N,O;!HQ]
    Family HBondDonor
    Weights 1.0
EndFeature
"""
    self.failUnlessRaises(ValueError,
                          lambda:ChemicalFeatures.BuildFeatureFactoryFromString(fdefBlock))
    fdefBlock = \
"""DefineFeature HDonor1 [N,O;!H0]
    Family HBondDonor
    Weights 1.0
"""
    self.failUnlessRaises(ValueError,
                          lambda:ChemicalFeatures.BuildFeatureFactoryFromString(fdefBlock))

    self.failUnlessRaises(IOError,
                          lambda:ChemicalFeatures.BuildFeatureFactory('noSuchFile.txt'))

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
    self.failUnless(cfac.GetNumFeatureDefs() == 2)
    mol = Chem.MolFromSmiles('n1ccccc1')
    feats = cfac.GetFeaturesForMol(mol)
    self.failUnless(len(feats)==2)
    m = ChemicalFeatures.GetAtomMatch(feats)
    self.failIf(m)

    mol = Chem.MolFromSmiles('c1ccccc1N')
    feats = cfac.GetFeaturesForMol(mol)
    self.failUnless(len(feats)==2)
    m = ChemicalFeatures.GetAtomMatch(feats)
    self.failUnless(len(m)==2)

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

    m = Chem.MolFromSmiles('O=CCCN')
    rdDistGeom.EmbedMolecule(m)
    feats = cfac.GetFeaturesForMol(m)
    for feat in feats:
      feat.GetPos()
    m = None
    for feat in feats:
      feat.GetPos()

if __name__ == '__main__':
  unittest.main()
