import unittest
from rdkit import Chem
from rdkit.Chem import rdShapeAlign
from rdkit import RDConfig

datadir = RDConfig.RDBaseDir + '/External/pubchem_shape/test_data'


class TestCase(unittest.TestCase):

  def setUp(self):
    suppl = Chem.SDMolSupplier(datadir + '/test1.sdf')
    self.ref = suppl[0]
    self.probe = suppl[1]

  def test1_Defaults(self):
    tpl = rdShapeAlign.AlignMol(self.ref, self.probe, opt_param=0.5, max_preiters=3,
                                max_postiters=16)
    self.assertAlmostEqual(tpl[0], 0.773, places=3)
    self.assertAlmostEqual(tpl[1], 0.303, places=3)

  def test2_NoColor(self):
    tpl = rdShapeAlign.AlignMol(self.ref, self.probe, useColors=False, opt_param=0.5,
                                max_preiters=3, max_postiters=16)
    self.assertAlmostEqual(tpl[0], 0.773, places=3)
    self.assertAlmostEqual(tpl[1], 0.0, places=3)

  def test3_FromShape(self):
    shp = rdShapeAlign.PrepareConformer(self.ref)
    self.assertTrue(type(shp) == rdShapeAlign.ShapeInput)
    tpl = rdShapeAlign.AlignMol(shp, self.probe, opt_param=0.5, max_preiters=3, max_postiters=16)
    self.assertAlmostEqual(tpl[0], 0.773, places=3)
    self.assertAlmostEqual(tpl[1], 0.303, places=3)

  def test4_ShapeInputOptions(self):
    opts = rdShapeAlign.ShapeInputOptions()
    opts.useColors = False
    shp = rdShapeAlign.PrepareConformer(self.ref, -1, opts)
    tpl = rdShapeAlign.AlignMol(shp, self.probe, opt_param=0.5, max_preiters=3, max_postiters=16)
    self.assertAlmostEqual(tpl[0], 0.773, places=3)
    self.assertAlmostEqual(tpl[1], 0.0, places=3)

    opts.atomSubset = [4, 5, 6, 7, 8, 9]
    shp = rdShapeAlign.PrepareConformer(self.ref, -1, opts)
    self.assertAlmostEqual(shp.sov, 251.946, places=3)
    self.assertAlmostEqual(shp.sof, 0.0, places=3)

    opts.atomRadii = [(4, 1.9)]
    shp = rdShapeAlign.PrepareConformer(self.ref, -1, opts)
    self.assertAlmostEqual(shp.sov, 274.576, places=3)

    with self.assertRaises(AttributeError):
      opts.rhubarb = True

  def test5_ShapeShapeOverlay(self):
    refShp = rdShapeAlign.PrepareConformer(self.ref, -1)
    probeShp = rdShapeAlign.PrepareConformer(self.probe, -1)
    tpl = rdShapeAlign.AlignShapes(refShp, probeShp)
    probeCp = Chem.Mol(self.probe)
    rdShapeAlign.TransformConformer(refShp.shift, tpl[2], probeShp, probeCp.GetConformer(-1))
    # Just show it did something.  The full test is in the C++ layer.
    self.assertNotEqual(self.probe.GetConformer().GetAtomPosition(0),
                        probeCp.GetConformer().GetAtomPosition(0))
    matrix = tpl[2][:10]
    with self.assertRaises(ValueError):
      rdShapeAlign.TransformConformer(refShp.shift, matrix, probeShp, probeCp.GetConformer(-1))

  def test6_notColorAtoms(self):
    m1 = Chem.MolFromSmiles("Nc1ccccc1 |(0.392086,-2.22477,0.190651;"
                            "0.232269,-1.38667,0.118385;-1.06274,-0.918982,0.0342466;"
                            "-1.26098,0.446053,-0.0811879;-0.244035,1.36265,-0.11691;"
                            "1.05134,0.875929,-0.031248;1.28797,-0.499563,0.0864097)"
                            ",atomProp:0.dummyLabel.*|")
    opts = rdShapeAlign.ShapeInputOptions()
    opts.notColorAtoms = [0]
    shp = rdShapeAlign.PrepareConformer(m1, -1, opts)
    self.assertAlmostEqual(shp.sof, 5.074, places=3)

if __name__ == '__main__':
  unittest.main()
