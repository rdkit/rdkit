#
#  Copyright (C) 2004-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import os
import unittest

from rdkit import Chem
from rdkit import Geometry
from rdkit import RDConfig
from rdkit.Chem import ChemicalFeatures, AllChem
from rdkit.Chem.Pharm3D import Pharmacophore


def feq(n1, n2, tol=1e-5):
  return abs(n1 - n2) <= tol


class TestCase(unittest.TestCase):

  def setUp(self):
    self.fdefBlock = """
                   DefineFeature HAcceptor1 [N,O;H0]
                      Family HBondAcceptor
                      Weights 1.0
                   EndFeature
                   DefineFeature HDonor1 [N,O;!H0]
                      Family HBondDonor
                      Weights 1.0
                   EndFeature
                   DefineFeature Aromatic1 c1ccccc1
                      Family Aromatic
                      Weights 1.0,1.0,1.0,1.0,1.0,1.0
                   EndFeature\n"""

    self.featFactory = ChemicalFeatures.BuildFeatureFactoryFromString(self.fdefBlock)
    self.feats = [ChemicalFeatures.FreeChemicalFeature('HBondAcceptor', 'HAcceptor1',
                                                       Geometry.Point3D(0.0, 0.0, 0.0)),
                  ChemicalFeatures.FreeChemicalFeature('HBondDonor', 'HDonor1',
                                                       Geometry.Point3D(2.65, 0.0, 0.0)),
                  ChemicalFeatures.FreeChemicalFeature('Aromatic', 'Aromatic1',
                                                       Geometry.Point3D(5.12, 0.908, 0.0)), ]
    self.pcophore = Pharmacophore.Pharmacophore(self.feats)
    self.radii = [1.2, 1.2, 1.6]
    self.explicit_pcophore = Pharmacophore.ExplicitPharmacophore(self.feats,
                                                                 self.radii)

  def test1Basics(self):
    pcophore = self.pcophore
    self.assertTrue(len(pcophore.getFeatures()) == 3)
    self.assertTrue(pcophore.getFeature(0))
    self.assertTrue(pcophore.getFeature(1))
    self.assertTrue(pcophore.getFeature(2))
    self.assertRaises(IndexError, pcophore.getFeature, 3)
    # print()
    # print(str(pcophore))
    self.assertIn('Aromatic', str(pcophore))

  def test2BoundSetting(self):
    pcophore = self.pcophore

    pcophore.setUpperBound(0, 1, 3.0)
    self.assertTrue(feq(pcophore.getUpperBound(0, 1), 3.0))
    self.assertTrue(feq(pcophore.getUpperBound(1, 0), 3.0))
    pcophore.setUpperBound(1, 0, 5.0)
    self.assertTrue(feq(pcophore.getUpperBound(0, 1), 5.0))
    self.assertTrue(feq(pcophore.getUpperBound(1, 0), 5.0))
    self.assertRaises(IndexError, pcophore.setUpperBound, 0, 3, 2.0)
    self.assertRaises(ValueError, pcophore.setUpperBound, 0, 3, 2.0, checkBounds=True)
    self.assertRaises(IndexError, pcophore.setUpperBound, 3, 0, 2.0)
    self.assertRaises(ValueError, pcophore.setUpperBound, 3, 0, 2.0, checkBounds=True)

    nfeatures = len(pcophore._feats)
    self.assertTrue(pcophore._checkBounds(0, 0))
    self.assertTrue(pcophore._checkBounds(0, nfeatures - 1))
    self.assertTrue(pcophore._checkBounds(nfeatures - 1, 0))
    self.assertTrue(pcophore._checkBounds(nfeatures - 1, nfeatures - 1))
    self.assertRaises(ValueError, pcophore._checkBounds, -1, 0)
    self.assertRaises(ValueError, pcophore._checkBounds, 0, -1)
    self.assertRaises(ValueError, pcophore._checkBounds, nfeatures, 0)
    self.assertRaises(ValueError, pcophore._checkBounds, 0, nfeatures)

    pcophore.setLowerBound(0, 1, 2.0)
    self.assertTrue(feq(pcophore.getLowerBound(0, 1), 2.0))
    self.assertTrue(feq(pcophore.getLowerBound(1, 0), 2.0))
    pcophore.setLowerBound(1, 0, 3.0)
    self.assertTrue(feq(pcophore.getLowerBound(0, 1), 3.0))
    self.assertTrue(feq(pcophore.getLowerBound(1, 0), 3.0))
    self.assertRaises(IndexError, pcophore.setLowerBound, 0, 3, 2.0)
    self.assertRaises(ValueError, pcophore.setLowerBound, 0, 3, 2.0, checkBounds=True)
    self.assertRaises(IndexError, pcophore.setLowerBound, 3, 0, 2.0)
    self.assertRaises(ValueError, pcophore.setLowerBound, 3, 0, 2.0, checkBounds=True)

  def test3Bound2DSetting(self):
    pcophore = self.pcophore

    pcophore.setUpperBound2D(0, 1, 3)
    self.assertTrue(pcophore.getUpperBound2D(0, 1) == 3)
    self.assertTrue(pcophore.getUpperBound2D(1, 0) == 3)
    pcophore.setUpperBound2D(1, 0, 5)
    self.assertTrue(pcophore.getUpperBound2D(0, 1) == 5)
    self.assertTrue(pcophore.getUpperBound2D(1, 0) == 5)
    self.assertRaises(IndexError, pcophore.setUpperBound2D, 0, 3, 2)
    self.assertRaises(ValueError, pcophore.setUpperBound2D, 0, 3, 2, checkBounds=True)
    self.assertRaises(IndexError, pcophore.setUpperBound2D, 3, 0, 2)
    self.assertRaises(ValueError, pcophore.setUpperBound2D, 3, 0, 2, checkBounds=True)

    pcophore.setLowerBound2D(0, 1, 3)
    self.assertTrue(pcophore.getLowerBound2D(0, 1) == 3)
    self.assertTrue(pcophore.getLowerBound2D(1, 0) == 3)
    pcophore.setLowerBound2D(1, 0, 5)
    self.assertTrue(pcophore.getLowerBound2D(0, 1) == 5)
    self.assertTrue(pcophore.getLowerBound2D(1, 0) == 5)
    self.assertRaises(IndexError, pcophore.setLowerBound2D, 0, 3, 2)
    self.assertRaises(ValueError, pcophore.setLowerBound2D, 0, 3, 2, checkBounds=True)
    self.assertRaises(IndexError, pcophore.setLowerBound2D, 3, 0, 2)
    self.assertRaises(ValueError, pcophore.setLowerBound2D, 3, 0, 2, checkBounds=True)

  def test4Github252(self):
    fdef = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
    feat_factory = ChemicalFeatures.BuildFeatureFactory(fdef)

    m1 = Chem.MolFromSmiles('Cc1ccccc1')
    feats = feat_factory.GetFeaturesForMol(m1)
    self.assertRaises(RuntimeError, lambda: Pharmacophore.Pharmacophore(feats))

    AllChem.Compute2DCoords(m1)
    Pharmacophore.Pharmacophore(feats)

  def testApplyRadiiToBounds(self):
    """
    Ensure that applying radii to pharmacophore feature bounds multiple times
    results in the same bounds.
    """
    pcophore2 = Pharmacophore.ExplicitPharmacophore(self.feats, self.radii)
    pcophore2.applyRadiiToBounds()
    self.assertEquals(str(pcophore2), str(self.explicit_pcophore))
    self.explicit_pcophore.applyRadiiToBounds()
    self.assertEquals(str(pcophore2), str(self.explicit_pcophore))

  def testExplicitPharmacophoreFeats(self):
    """
    Ensure that applying radii of 0.0 to the bounds results in the same
    bounds matrix as a regular Pharmacophore.
    """
    explicit_pcophore = self.explicit_pcophore
    self.assertListEqual(self.pcophore.getFeatures(), explicit_pcophore.getFeatures())

    for i in range(len(explicit_pcophore.getFeatures())):
      explicit_pcophore.setRadius(i, 0.0)
    explicit_pcophore.applyRadiiToBounds()
    self.assertEquals(str(explicit_pcophore), str(self.pcophore))


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
