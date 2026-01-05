import unittest
import numpy as np

from rdkit import Chem
from rdkit.Chem import rdShapeAlign2, rdMolTransforms
from rdkit import RDConfig

datadir = RDConfig.RDBaseDir + '/External/pubchem_shape/test_data'


class TestCase(unittest.TestCase):

  def setUp(self):
    suppl = Chem.SDMolSupplier(datadir + '/test1.sdf')
    self.ref = suppl[0]
    self.probe = suppl[1]

  def test1_Defaults(self):
    tpl = rdShapeAlign2.AlignMol(self.ref, self.probe)
    self.assertAlmostEqual(tpl[0], 0.701, places=3)
    self.assertAlmostEqual(tpl[1], 0.556, places=3)
    
  def test2_FromShape(self):
    opts = rdShapeAlign2.ShapeOverlayOptions()
    shp = rdShapeAlign2.ShapeInput(self.ref, -1, opts)
    self.assertAlmostEqual(shp.GetVolume(), 591.1, places=1)
    self.assertAlmostEqual(shp.GetColorVolume(), 16.1, places=1)
    types = shp.GetTypes()
    self.assertEqual(len(types), shp.GetNumAtoms() + shp.GetNumFeatures())

    tpl = rdShapeAlign2.AlignMol(shp, self.probe)
    self.assertAlmostEqual(tpl[0], 0.701, places=3)
    self.assertAlmostEqual(tpl[1], 0.556, places=3)
    
  def test3_OverlayShapes(self):
    opts = rdShapeAlign2.ShapeOverlayOptions()
    refShp = rdShapeAlign2.ShapeInput(self.ref, -1, opts)
    probeShp = rdShapeAlign2.ShapeInput(self.probe, -1, opts)
    tpl = rdShapeAlign2.AlignShape(refShp, probeShp)
    self.assertAlmostEqual(tpl[0], 0.701, places=3)
    self.assertAlmostEqual(tpl[1], 0.556, places=3)
    self.assertEqual(len(tpl[2]), 16)
    expXform = [-0.864046,-0.281027,-0.419039,-15.338341,
                -0.026639,-0.803610,0.595666,-5.658636,
                -0.504465,0.524967,0.686617,2.000176,
                0.000000,0.000000,0.000000,1.000000]
    for t, et in zip(tpl[2], expXform):
      self.assertAlmostEqual(t, et, places=6)

    # Check that the transformation works.  The input
    # structure isn't in canonical order, so its
    # first atom isn't the first one in the CXSmiles.
    ttrans = [tpl[2][0:4], tpl[2][4:8], tpl[2][8:12], tpl[2][12:16]]
    nptrans = np.array(ttrans)
    rdMolTransforms.TransformConformer(self.probe.GetConformer(), nptrans)
    pos1 = self.probe.GetConformer().GetAtomPosition(0)
    self.assertAlmostEqual(pos1.x, -13.57295, places=5)
    self.assertAlmostEqual(pos1.y, -4.87636, places=5)
    self.assertAlmostEqual(pos1.z, 4.55432, places=5)

  def test4_scoreOnly(self):
    # These are 2gu8_lig_796 and 3ama_lig_SKE from Ilenia Giangreco's
    # ligand overlay set.
    mol1 = Chem.MolFromSmiles("CNc1nccc(-c2ccc(C(=O)N[C@H](C[NH3+])Cc3ccc(Cl)cc3Cl)s2)n1 |(-3.234,-9.922,0.988;-2.567,-10.239,2.232;-3.313,-10.306,3.299;-2.759,-10.404,4.518;-3.487,-10.471,5.637;-4.885,-10.444,5.491;-5.456,-10.373,4.213;-6.935,-10.268,4.033;-8.015,-10.254,4.915;-9.243,-10.001,4.289;-9.142,-9.794,2.927;-10.094,-9.522,1.832;-9.69,-9.261,0.684;-11.364,-9.586,2.201;-12.392,-8.987,1.375;-12.597,-7.523,1.75;-13.369,-6.797,0.753;-13.592,-9.862,1.699;-13.475,-11.196,0.951;-13.875,-11.306,-0.371;-13.79,-12.54,-1.033;-13.307,-13.645,-0.365;-13.168,-15.231,-1.141;-12.885,-13.567,0.922;-12.966,-12.343,1.553;-12.406,-12.346,3.24;-7.476,-9.917,2.366;-4.667,-10.294,3.115),wU:14.14|")
    mol2 = Chem.MolFromSmiles("Nc1nc(Nc2ccc(S(N)(=O)=O)cc2)nn1C(=O)c1c(F)cccc1F |(-2.664,-10.133,1.814;-3.368,-9.918,2.943;-2.844,-9.878,4.165;-3.873,-9.66,4.996;-3.747,-9.546,6.343;-4.743,-9.498,7.255;-6.094,-9.327,6.912;-7.084,-9.291,7.904;-6.736,-9.42,9.245;-8.016,-9.385,10.489;-8.886,-8.022,10.364;-7.378,-9.461,11.824;-8.927,-10.558,10.251;-5.396,-9.6,9.596;-4.402,-9.634,8.612;-4.987,-9.539,4.267;-4.691,-9.722,3.056;-5.596,-9.65,2.057;-5.285,-9.738,0.861;-7.013,-9.24,2.504;-7.724,-10.046,3.399;-7.192,-11.212,3.821;-8.993,-9.652,3.839;-9.533,-8.441,3.398;-8.813,-7.631,2.513;-7.551,-8.029,2.069;-6.85,-7.237,1.224)|")
    opts = rdShapeAlign2.ShapeOverlayOptions()
    shape_2gu8 = rdShapeAlign2.ShapeInput(mol1, -1, opts)
    shape_3ama = rdShapeAlign2.ShapeInput(mol2, -1, opts)
    scores = rdShapeAlign2.ScoreShape(shape_2gu8, shape_3ama)
    self.assertAlmostEqual(scores[0], 0.262113, places=6)
    self.assertAlmostEqual(scores[1], 0.355309, places=6)
  
    scores = rdShapeAlign2.ScoreMol(shape_2gu8, mol2)
    self.assertAlmostEqual(scores[0], 0.262113, places=6)
    self.assertAlmostEqual(scores[1], 0.355309, places=6)

    scores = rdShapeAlign2.ScoreMol(mol1, mol2)
    self.assertAlmostEqual(scores[0], 0.262113, places=6)
    self.assertAlmostEqual(scores[1], 0.355309, places=6)

    
if __name__ == '__main__':
  unittest.main()
