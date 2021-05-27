from rdkit import Chem
from rdkit import rdBase
from rdkit import RDConfig
import os

from rdkit.Chem import rdMolDescriptors as rdMD
from rdkit.Chem import AllChem

haveDescrs3D = hasattr(rdMD, 'CalcAUTOCORR3D')

import time, unittest


def _gen3D(m, is3d, calculator):
  if not is3d:
    m = Chem.AddHs(m)
    ps = AllChem.ETKDG()
    ps.randomSeed = 0xf00d
    AllChem.EmbedMolecule(m, ps)
  return calculator(m)


class TestCase(unittest.TestCase):

  def setUp(self):
    self.dataDir = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Descriptors', 'test_data')
    self.suppl = Chem.SDMolSupplier(os.path.join(self.dataDir, 'PBF_egfr.sdf'), removeHs=False)

  @unittest.skipIf(not haveDescrs3D, "3d descriptors not present")
  def test1AUTOCORR2D(self):
    # not really a 3D descriptor, but this was added at the same time
    with open(os.path.join(self.dataDir, 'auto2D.out')) as refFile:
      for i, m in enumerate(self.suppl):
        if i > 10:
          break
        nm = m.GetProp('_Name')
        inl = refFile.readline()
        split = inl.split('\t')
        self.assertEqual(split[0], nm)
        split.pop(0)
        vs = rdMD.CalcAUTOCORR2D(m)
        for rv, nv in zip(split, vs):
          self.assertAlmostEqual(float(rv), nv, delta=0.05)

  @unittest.skipIf(not haveDescrs3D, "3d descriptors not present")
  def test2AUTOCORR3D(self):
    with open(os.path.join(self.dataDir, 'auto3D_dragon.out')) as refFile:
      for i, m in enumerate(self.suppl):
        if i > 10:
          break
        nm = m.GetProp('_Name')
        inl = refFile.readline()
        split = inl.split('\t')
        self.assertEqual(split[0], nm)
        split.pop(0)
        vs = _gen3D(m, True, rdMD.CalcAUTOCORR3D)
        for rv, nv in zip(split, vs):
          self.assertAlmostEqual(float(rv), nv, delta=0.05)

  @unittest.skipIf(not haveDescrs3D, "3d descriptors not present")
  def test3GETAWAY(self):
    with open(os.path.join(self.dataDir, 'GETAWAY.new.out')) as refFile:
      for i, m in enumerate(self.suppl):
        if i > 10:
          break
        nm = m.GetProp('_Name')
        inl = refFile.readline()
        split = inl.split('\t')
        self.assertEqual(split[0], nm)
        split.pop(0)
        vs = _gen3D(m, True, rdMD.CalcGETAWAY)
        for rv, nv in zip(split, vs):
          self.assertAlmostEqual(float(rv), nv, delta=0.05)

  @unittest.skipIf(not haveDescrs3D, "3d descriptors not present")
  def test4MORSE(self):
    with open(os.path.join(self.dataDir, 'MORSE.out')) as refFile:
      for i, m in enumerate(self.suppl):
        if i > 10:
          break
        nm = m.GetProp('_Name')
        inl = refFile.readline()
        split = inl.split('\t')
        self.assertEqual(split[0], nm)
        split.pop(0)
        vs = _gen3D(m, True, rdMD.CalcMORSE)
        for rv, nv in zip(split, vs):
          ref = float(rv)
          self.assertTrue(ref < 1 or abs(ref - nv) / ref < 0.02)

  @unittest.skipIf(not haveDescrs3D, "3d descriptors not present")
  def test5RDF(self):
    with open(os.path.join(self.dataDir, 'RDF.out')) as refFile:
      for i, m in enumerate(self.suppl):
        if i > 10:
          break
        nm = m.GetProp('_Name')
        inl = refFile.readline()
        split = inl.split('\t')
        self.assertEqual(split[0], nm)
        split.pop(0)
        vs = _gen3D(m, True, rdMD.CalcRDF)
        for rv, nv in zip(split, vs):
          ref = float(rv)
          self.assertTrue(ref < 0.5 or abs(ref - nv) / ref < 0.02)

  @unittest.skipIf(not haveDescrs3D, "3d descriptors not present")
  def test6WHIM(self):
    with open(os.path.join(self.dataDir, 'whim.new.out')) as refFile:
      for i, m in enumerate(self.suppl):
        if i > 10:
          break
        nm = m.GetProp('_Name')
        inl = refFile.readline()
        split = inl.split('\t')
        self.assertEqual(split[0], nm)
        split.pop(0)
        vs = _gen3D(m, True, lambda x: rdMD.CalcWHIM(x, thresh=0.01))
        for rv, nv in zip(split, vs):
          self.assertAlmostEqual(float(rv), nv, delta=0.01)

  @unittest.skipIf(not haveDescrs3D, "3d descriptors not present")
  def testGithub2037(self):
    m = Chem.AddHs(Chem.MolFromSmiles("CCCCCCC"))
    cids = AllChem.EmbedMultipleConfs(m, 10)
    # start with defaults (which does not cache results):
    npr1s = []
    npr2s = []
    for cid in cids:
      npr1s.append(rdMD.CalcNPR1(m, confId=cid))
      npr2s.append(rdMD.CalcNPR2(m, confId=cid))
    for i in range(1, len(npr1s)):
      self.assertNotAlmostEqual(npr1s[0], npr1s[i])
      self.assertNotAlmostEqual(npr2s[0], npr2s[i])

    # now ensure that we can cache:
    npr1s = []
    npr2s = []
    for cid in cids:
      npr1s.append(rdMD.CalcNPR1(m, confId=cid, force=False))
      npr2s.append(rdMD.CalcNPR2(m, confId=cid, force=False))
    for i in range(1, len(npr1s)):
      self.assertAlmostEqual(npr1s[0], npr1s[i])
      self.assertAlmostEqual(npr2s[0], npr2s[i])

  @unittest.skipIf(not haveDescrs3D, "3d descriptors not present")
  def testGithub4167(self):
    with Chem.SDMolSupplier(os.path.join(self.dataDir, 'github4167.sdf'), removeHs=False,
                            sanitize=True) as suppl:
      m1 = suppl[0]
      m2 = suppl[1]
    m1.AddConformer(Chem.Conformer(m2.GetConformer()), assignId=True)
    v1_0 = rdMD.CalcSpherocityIndex(m1)
    v1_1 = rdMD.CalcSpherocityIndex(m1, confId=1, force=True)
    v2 = rdMD.CalcSpherocityIndex(m2)
    self.assertNotEqual(v1_0, v1_1)
    self.assertEqual(v1_1, v2)


if (__name__ == '__main__'):
  unittest.main()
