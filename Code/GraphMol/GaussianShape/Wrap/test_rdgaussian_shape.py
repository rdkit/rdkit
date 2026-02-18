import unittest
import numpy as np

from rdkit import Chem
from rdkit.Chem import rdGaussianShape, rdMolTransforms
from rdkit import RDConfig
from rdkit.Geometry import Point3D


datadir = RDConfig.RDBaseDir + '/External/pubchem_shape/test_data'


class TestCase(unittest.TestCase):

  def setUp(self):
    suppl = Chem.SDMolSupplier(datadir + '/test1.sdf')
    self.ref = suppl[0]
    self.probe = suppl[1]

  def test1_Defaults(self):
    tpl = rdGaussianShape.AlignMol(self.ref, self.probe)
    self.assertAlmostEqual(tpl[0], 0.497, places=3)
    self.assertAlmostEqual(tpl[1], 0.759, places=3)
    self.assertAlmostEqual(tpl[2], 0.235, places=3)
    

  def test2_NoColor(self):
    ovOpts = rdGaussianShape.ShapeOverlayOptions()
    ovOpts.optimMode = rdGaussianShape.OptimMode.SHAPE_ONLY
    shpOpts=  rdGaussianShape.ShapeInputOptions()
    shpOpts.useColors = False
    tpl = rdGaussianShape.AlignMol(self.ref, self.probe , shpOpts, shpOpts, ovOpts)
    self.assertAlmostEqual(tpl[0], 0.760, places=3)
    self.assertAlmostEqual(tpl[1], 0.760, places=3)
    self.assertAlmostEqual(tpl[2], 0.0, places=3)

  def test3_FromShape(self):
    ovOpts = rdGaussianShape.ShapeOverlayOptions()
    shpOpts=  rdGaussianShape.ShapeInputOptions()
    shp = rdGaussianShape.ShapeInput(self.ref, -1, shpOpts, ovOpts)
    self.assertAlmostEqual(shp.ShapeVolume, 591.058, places=3)
    self.assertAlmostEqual(shp.ColorVolume, 31.935, places=3)
    self.assertTrue(type(shp) == rdGaussianShape.ShapeInput)
    tpl = rdGaussianShape.AlignMol(shp, self.probe)
    self.assertAlmostEqual(tpl[0], 0.497, places=3)
    self.assertAlmostEqual(tpl[1], 0.759, places=3)
    self.assertAlmostEqual(tpl[2], 0.235, places=3)

  def test4_customFeatures(self):
    m1 = Chem.MolFromSmiles(
      "O=CC=O |(-1.75978,0.148897,0;-0.621382,-0.394324,0;0.624061,0.3656,.1;1.7571,-0.120174,.1)|")
    opts = rdGaussianShape.ShapeInputOptions()
    opts.customFeatures = ((1, Point3D(-1.75978, 0.148897,
                                       0), 1.0), (2, Point3D(1.7571, -0.120174, 0.1), 1.0))
    ovOpts = rdGaussianShape.ShapeOverlayOptions()
    shp = rdGaussianShape.ShapeInput(m1, -1, opts, ovOpts)
    self.assertEqual(shp.NumAtoms, 4)
    self.assertEqual(shp.NumFeatures, 2)
    m2 = Chem.Mol(m1)
    opts2 = rdGaussianShape.ShapeInputOptions()
    opts2.customFeatures = ((2, Point3D(-1.75978, 0.148897,
                                        0), 1.0), (1, Point3D(1.7571, -0.120174, 0.1), 1.0))
    shp2 = rdGaussianShape.ShapeInput(m2, -1, opts2, ovOpts)
    tpl = rdGaussianShape.AlignShapes(shp, shp2, ovOpts)
    self.assertAlmostEqual(tpl[0], 0.999, places=3)
    self.assertAlmostEqual(tpl[1], 1.000, places=3)
    self.assertAlmostEqual(tpl[2], 0.999, places=3)
    tf = tpl[3]
    self.assertGreater(0.0, tf[0])
    self.assertEqual(1.0, tf[15])

    # check the getter:
    cfs = opts2.customFeatures
    self.assertEqual(len(cfs), 2)
    self.assertEqual(cfs[0][0], 2)
    self.assertEqual(cfs[1][0], 1)

  def test5_customFeatures(self):
    m1 = Chem.MolFromSmiles(
      "O=CC=O |(-1.75978,0.148897,0;-0.621382,-0.394324,0;0.624061,0.3656,.1;1.7571,-0.120174,.1)|")
    opts = rdGaussianShape.ShapeInputOptions()
    opts.customFeatures = ((1, Point3D(-1.75978, 0.148897,
                                       0), 1.0), (2, Point3D(1.7571, -0.120174, 0.1), 1.0))
    m2 = Chem.Mol(m1)
    opts2 = rdGaussianShape.ShapeInputOptions()
    opts2.customFeatures = ((2, Point3D(-1.75978, 0.148897,
                                        0), 1.0), (1, Point3D(1.7571, -0.120174, 0.1), 1.0))
    ovOpts = rdGaussianShape.ShapeOverlayOptions()
    tpl = rdGaussianShape.AlignMol(m1, m2, opts, opts2, ovOpts)
    self.assertAlmostEqual(tpl[0], 0.999, places=3)
    self.assertAlmostEqual(tpl[1], 1.000, places=3)
    self.assertAlmostEqual(tpl[2], 0.999, places=3)

  def test6_FixedScore(self):
    ovOpts = rdGaussianShape.ShapeOverlayOptions()
    # Just to make sure it's there and returns a value.
    opts = rdGaussianShape.ShapeInputOptions()
    tpl = rdGaussianShape.ScoreMol(self.ref, self.ref, opts, opts, ovOpts)
    self.assertAlmostEqual(tpl[0], 1.0, places=3)
    self.assertAlmostEqual(tpl[1], 1.0, places=3)
    self.assertAlmostEqual(tpl[2], 1.0, places=3)

    opts = rdGaussianShape.ShapeInputOptions()
    opts.useColors = False
    ovOpts.normalize = False
    shp = rdGaussianShape.ShapeInput(self.ref, -1, opts, ovOpts)
    tpl = rdGaussianShape.ScoreMol(shp, self.probe, opts)
    self.assertAlmostEqual(tpl[0], 0.0, places=3)
    self.assertAlmostEqual(tpl[1], 0.0, places=3)
    self.assertAlmostEqual(tpl[2], 0.0, places=3)
    
    opts.useColors = True
    shp1 = rdGaussianShape.ShapeInput(self.probe, -1, opts, ovOpts)
    shp2 = rdGaussianShape.ShapeInput(self.probe, -1, opts, ovOpts)
    tpl = rdGaussianShape.ScoreShape(shp1, shp2, ovOpts)
    self.assertAlmostEqual(tpl[0], 1.0, places=3)
    self.assertAlmostEqual(tpl[1], 1.0, places=3)
    self.assertAlmostEqual(tpl[2], 1.0, places=3)

    
if __name__ == '__main__':
  unittest.main()
