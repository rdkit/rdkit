import os
import unittest

import numpy

from rdkit import Chem, RDConfig
from rdkit.Chem import ChemicalForceFields, rdDistGeom


def feq(v1, v2, tol2=1e-4):
  return abs(v1 - v2) <= tol2


class TestCase(unittest.TestCase):

  def setUp(self):
    self.dirName = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'ForceFieldHelpers', 'UFF',
                                'test_data')

  def test1(self):
    fName = os.path.join(self.dirName, 'benzene.mol')
    m = Chem.MolFromMolFile(fName)
    self.assertFalse(ChemicalForceFields.UFFOptimizeMolecule(m))
    # make sure that keyword arguments work:
    m = Chem.MolFromMolFile(fName)
    self.assertTrue(ChemicalForceFields.UFFOptimizeMolecule(m, maxIters=1))

    m = Chem.MolFromMolFile(fName)
    self.assertFalse(ChemicalForceFields.UFFOptimizeMolecule(m, vdwThresh=2.0))

    m = Chem.MolFromMolFile(fName)
    self.assertFalse(ChemicalForceFields.UFFOptimizeMolecule(m, confId=-1))

    m = Chem.MolFromMolFile(fName)
    self.assertRaises(ValueError, lambda: ChemicalForceFields.UFFOptimizeMolecule(m, confId=1))

  def test2(self):
    fName = os.path.join(self.dirName, 'benzene.mol')
    m = Chem.MolFromMolFile(fName)
    ff = ChemicalForceFields.UFFGetMoleculeForceField(m)
    self.assertTrue(ff)
    e1 = ff.CalcEnergy()
    r = ff.Minimize()
    self.assertTrue(r == 0)
    e2 = ff.CalcEnergy()
    self.assertTrue(e2 < e1)

    # test keyword args:
    r = ff.Minimize(forceTol=1e-8)
    self.assertTrue(r == 0)

    # test keyword args:
    r = ff.Minimize(energyTol=1e-3)
    self.assertTrue(r == 0)

  def test3(self):
    molB = """


  4  4  0  0  0  0  0  0  0  0999 V2000
   -0.8500    0.4512   -0.6671 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3307   -0.9436   -0.3641 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6796   -0.4074    0.5894 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5011    0.8998   -0.1231 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
  1  4  1  0
M  END"""
    m = Chem.MolFromMolBlock(molB)

    ff = ChemicalForceFields.UFFGetMoleculeForceField(m)
    self.assertTrue(ff)
    e1 = ff.CalcEnergy()
    r = ff.Minimize()
    self.assertTrue(r == 0)
    e2 = ff.CalcEnergy()
    self.assertTrue(e2 < e1)

  def test4(self):
    m = Chem.MolFromSmiles('[Cu](C)(C)(C)(C)C')
    self.assertFalse(ChemicalForceFields.UFFHasAllMoleculeParams(m))

    m = Chem.MolFromSmiles('C(C)(C)(C)C')
    self.assertTrue(ChemicalForceFields.UFFHasAllMoleculeParams(m))

  def test5(self):
    fName = os.path.join(self.dirName, 'benzene.mol')
    m = Chem.MolFromMolFile(fName)
    self.assertFalse(ChemicalForceFields.MMFFOptimizeMolecule(m))
    # make sure that keyword arguments work:
    m = Chem.MolFromMolFile(fName)
    self.assertTrue(ChemicalForceFields.MMFFOptimizeMolecule(m, maxIters=1))

    m = Chem.MolFromMolFile(fName)
    self.assertFalse(ChemicalForceFields.MMFFOptimizeMolecule(m, nonBondedThresh=2.0))

    m = Chem.MolFromMolFile(fName)
    self.assertFalse(ChemicalForceFields.MMFFOptimizeMolecule(m, confId=-1))

    m = Chem.MolFromMolFile(fName)
    self.assertRaises(ValueError, lambda: ChemicalForceFields.MMFFOptimizeMolecule(m, confId=1))

  def test6(self):
    fName = os.path.join(self.dirName, 'benzene.mol')
    m = Chem.MolFromMolFile(fName)
    mp = ChemicalForceFields.MMFFGetMoleculeProperties(m)
    ff = ChemicalForceFields.MMFFGetMoleculeForceField(m, mp)
    self.assertTrue(ff)
    e1 = ff.CalcEnergy()
    r = ff.Minimize()
    self.assertTrue(r == 0)
    e2 = ff.CalcEnergy()
    self.assertTrue(e2 < e1)

    # test keyword args:
    r = ff.Minimize(forceTol=1.0e-8)
    self.assertTrue(r == 0)

    # test keyword args:
    r = ff.Minimize(energyTol=1.0e-3)
    self.assertTrue(r == 0)

  def test7(self):
    molB = """


  4  4  0  0  0  0  0  0  0  0999 V2000
   -0.8500    0.4512   -0.6671 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3307   -0.9436   -0.3641 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6796   -0.4074    0.5894 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5011    0.8998   -0.1231 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
  1  4  1  0
M  END"""
    m = Chem.MolFromMolBlock(molB)

    mp = ChemicalForceFields.MMFFGetMoleculeProperties(m)
    ff = ChemicalForceFields.MMFFGetMoleculeForceField(m, mp)
    self.assertTrue(ff)
    e1 = ff.CalcEnergy()
    r = ff.Minimize()
    self.assertTrue(r == 0)
    e2 = ff.CalcEnergy()
    self.assertTrue(e2 < e1)

  def test8(self):
    m = Chem.MolFromSmiles('[Cu](C)(C)(C)(C)C')
    self.assertFalse(ChemicalForceFields.MMFFHasAllMoleculeParams(m))

    m = Chem.MolFromSmiles('C(C)(C)(C)C')
    self.assertTrue(ChemicalForceFields.MMFFHasAllMoleculeParams(m))

  def test9(self):
    m = Chem.MolFromSmiles('c1ccccc1CCNN')
    m = Chem.AddHs(m)
    mp = ChemicalForceFields.MMFFGetMoleculeProperties(m)
    mmffBondStretchParams = mp.GetMMFFBondStretchParams(m, 6, 7)
    self.assertTrue(mmffBondStretchParams)
    self.assertTrue((mmffBondStretchParams[0] == 0)
                    and (int(round(mmffBondStretchParams[1] * 1000) == 4258))
                    and (int(round(mmffBondStretchParams[2] * 1000) == 1508)))
    mmffBondStretchParams = mp.GetMMFFBondStretchParams(m, 0, 7)
    self.assertFalse(mmffBondStretchParams)
    mmffAngleBendParams = mp.GetMMFFAngleBendParams(m, 6, 7, 8)
    self.assertTrue(mmffAngleBendParams)
    self.assertTrue((mmffAngleBendParams[0] == 0)
                    and (int(round(mmffAngleBendParams[1] * 1000) == 777))
                    and (int(round(mmffAngleBendParams[2] * 1000) == 108290)))
    mmffAngleBendParams = mp.GetMMFFAngleBendParams(m, 0, 7, 8)
    self.assertFalse(mmffAngleBendParams)
    mmffStretchBendParams = mp.GetMMFFStretchBendParams(m, 6, 7, 8)
    self.assertTrue(mmffStretchBendParams)
    self.assertTrue((mmffStretchBendParams[0] == 0)
                    and (int(round(mmffStretchBendParams[1] * 1000) == 136))
                    and (int(round(mmffStretchBendParams[2] * 1000) == 282)))
    mmffStretchBendParams = mp.GetMMFFStretchBendParams(m, 0, 7, 8)
    self.assertFalse(mmffStretchBendParams)
    mmffTorsionParams = mp.GetMMFFTorsionParams(m, 6, 7, 8, 9)
    self.assertTrue(mmffTorsionParams)
    self.assertTrue((mmffTorsionParams[0] == 0) and (int(round(mmffTorsionParams[1] * 1000) == 0))
                    and (int(round(mmffTorsionParams[2] * 1000) == -300))
                    and (int(round(mmffTorsionParams[3] * 1000) == 500)))
    mmffTorsionParams = mp.GetMMFFTorsionParams(m, 0, 7, 8, 9)
    self.assertFalse(mmffTorsionParams)
    mmffOopBendParams = mp.GetMMFFOopBendParams(m, 6, 5, 4, 0)
    self.assertTrue(mmffOopBendParams)
    self.assertTrue(int(round(mmffOopBendParams * 1000)) == 40)
    mmffOopBendParams = mp.GetMMFFOopBendParams(m, 6, 5, 4, 1)
    self.assertFalse(mmffOopBendParams)
    sub1 = m.GetSubstructMatch(Chem.MolFromSmarts('NN[H]'))
    self.assertTrue(len(sub1) == 3)
    nIdx = sub1[0]
    hIdx = sub1[2]
    mmffVdWParams = mp.GetMMFFVdWParams(nIdx, hIdx)
    self.assertTrue(mmffVdWParams)
    self.assertTrue((int(round(mmffVdWParams[0] * 1000)) == 3321)
                    and (int(round(mmffVdWParams[1] * 1000)) == 34)
                    and (int(round(mmffVdWParams[2] * 1000)) == 2657)
                    and (int(round(mmffVdWParams[3] * 1000)) == 17))

  def test10(self):
    m = Chem.MolFromSmiles('c1ccccc1CCNN')
    m = Chem.AddHs(m)
    uffBondStretchParams = ChemicalForceFields.GetUFFBondStretchParams(m, 6, 7)
    self.assertTrue(uffBondStretchParams)
    self.assertTrue((int(round(uffBondStretchParams[0] * 1000) == 699592))
                    and (int(round(uffBondStretchParams[1] * 1000) == 1514)))
    uffBondStretchParams = ChemicalForceFields.GetUFFBondStretchParams(m, 0, 7)
    self.assertFalse(uffBondStretchParams)
    uffAngleBendParams = ChemicalForceFields.GetUFFAngleBendParams(m, 6, 7, 8)
    self.assertTrue(uffAngleBendParams)
    self.assertTrue((int(round(uffAngleBendParams[0] * 1000) == 303297))
                    and (int(round(uffAngleBendParams[1] * 1000) == 109470)))
    uffAngleBendParams = ChemicalForceFields.GetUFFAngleBendParams(m, 0, 7, 8)
    self.assertFalse(uffAngleBendParams)
    uffTorsionParams = ChemicalForceFields.GetUFFTorsionParams(m, 6, 7, 8, 9)
    self.assertTrue(uffTorsionParams)
    self.assertTrue((int(round(uffTorsionParams * 1000) == 976)))
    uffTorsionParams = ChemicalForceFields.GetUFFTorsionParams(m, 0, 7, 8, 9)
    self.assertFalse(uffTorsionParams)
    uffInversionParams = ChemicalForceFields.GetUFFInversionParams(m, 6, 5, 4, 0)
    self.assertTrue(uffInversionParams)
    self.assertTrue(int(round(uffInversionParams * 1000)) == 2000)
    uffInversionParams = ChemicalForceFields.GetUFFInversionParams(m, 6, 5, 4, 1)
    self.assertFalse(uffInversionParams)
    uffVdWParams = ChemicalForceFields.GetUFFVdWParams(m, 0, 9)
    self.assertTrue(uffVdWParams)
    self.assertTrue((int(round(uffVdWParams[0] * 1000)) == 3754)
                    and (int(round(uffVdWParams[1] * 1000)) == 85))

  def test11(self):
    query = Chem.MolFromSmarts('c1cccn1')
    for i in [0, 1]:
      m = Chem.MolFromSmiles('Cc1nc(=O)c(C[NH3+])c(-c2c[nH]c3ccccc23)[nH]1')
      aromaticFlagsBefore = []
      for a in m.GetAtoms():
        aromaticFlagsBefore.append(a.GetIsAromatic())
      if (i):
        self.assertTrue(ChemicalForceFields.MMFFGetMoleculeProperties(m))
      else:
        self.assertTrue(ChemicalForceFields.MMFFHasAllMoleculeParams(m))
      aromaticFlagsAfter = []
      for a in m.GetAtoms():
        aromaticFlagsAfter.append(a.GetIsAromatic())
      res = (aromaticFlagsBefore == aromaticFlagsAfter)
      if (i):
        res = not res
      self.assertTrue(res)
      pyrroleNIdx = m.GetSubstructMatch(query)[-1]
      self.assertTrue(m.GetAtomWithIdx(pyrroleNIdx).GetTotalDegree() == 3)

  def test12(self):
    self.dirName = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'ForceFieldHelpers', 'UFF',
                                'test_data')
    fName = os.path.join(self.dirName, 'Issue239.mol')
    m = Chem.MolFromMolFile(fName)
    self.assertIsNotNone(m)
    ff = ChemicalForceFields.UFFGetMoleculeForceField(m)
    ff.Initialize()
    positions = ff.Positions()
    savedPos = list(positions)
    e1 = ff.CalcEnergy()
    ff.Minimize(10000, 1.0e-6, 1.0e-3)
    e2 = ff.CalcEnergy()
    self.assertTrue(e2 < e1)
    e3 = ff.CalcEnergy(savedPos)
    self.assertAlmostEqual(e3, e1, 2)
    savedPos = tuple(positions)
    e3 = ff.CalcEnergy(savedPos)
    self.assertAlmostEqual(e3, e1, 2)
    savedPos = tuple(numpy.array(savedPos))
    e3 = ff.CalcEnergy(savedPos)
    self.assertAlmostEqual(e3, e1, 2)

  def test13(self):
    self.dirName = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'ForceFieldHelpers', 'UFF',
                                'test_data')
    fName = os.path.join(self.dirName, 'Issue239.mol')
    m = Chem.MolFromMolFile(fName)
    self.assertIsNotNone(m)
    mp = ChemicalForceFields.MMFFGetMoleculeProperties(m)
    ff = ChemicalForceFields.MMFFGetMoleculeForceField(m, mp)
    ff.Initialize()
    positions = ff.Positions()
    savedPos = list(positions)
    e1 = ff.CalcEnergy()
    ff.Minimize(10000, 1.0e-6, 1.0e-3)
    e2 = ff.CalcEnergy()
    self.assertTrue(e2 < e1)
    e3 = ff.CalcEnergy(savedPos)
    self.assertAlmostEqual(e3, e1, 2)

  def test14(self):
    self.dirName = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'ForceFieldHelpers', 'UFF',
                                'test_data')
    fName = os.path.join(self.dirName, 'Issue239.mol')
    m = Chem.MolFromMolFile(fName)
    self.assertIsNotNone(m)
    ff = ChemicalForceFields.UFFGetMoleculeForceField(m)
    ff.Initialize()
    positions = ff.Positions()
    savedPos = list(positions)
    grad1 = ff.CalcGrad()
    for v in grad1:
      self.assertNotAlmostEqual(v, 0.0, 3)
    ff.Minimize(10000, 1.0e-6, 1.0e-3)
    grad2 = ff.CalcGrad()
    for v in grad2:
      self.assertAlmostEqual(v, 0.0, 3)
    ff.Initialize()
    grad2 = ff.CalcGrad(savedPos)
    self.assertEqual(len(grad1), len(grad2))
    for i in range(len(grad1)):
      self.assertAlmostEqual(grad1[i], grad2[i], 3)

  def test15(self):
    self.dirName = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'ForceFieldHelpers', 'UFF',
                                'test_data')
    fName = os.path.join(self.dirName, 'Issue239.mol')
    m = Chem.MolFromMolFile(fName)
    self.assertIsNotNone(m)
    mp = ChemicalForceFields.MMFFGetMoleculeProperties(m)
    ff = ChemicalForceFields.MMFFGetMoleculeForceField(m, mp)
    ff.Initialize()
    positions = ff.Positions()
    savedPos = list(positions)
    grad1 = ff.CalcGrad()
    for v in grad1:
      self.assertNotAlmostEqual(v, 0.0, 3)
    ff.Minimize(10000, 1.0e-6, 1.0e-3)
    grad2 = ff.CalcGrad()
    for v in grad2:
      self.assertAlmostEqual(v, 0.0, 3)
    ff.Initialize()
    grad2 = ff.CalcGrad(savedPos)
    self.assertEqual(len(grad1), len(grad2))
    for i in range(len(grad1)):
      self.assertAlmostEqual(grad1[i], grad2[i], 3)

  def testGitHub2820(self):
    m = Chem.MolFromSmiles("[Na]C")
    self.assertIsNotNone(m)
    mp = ChemicalForceFields.MMFFGetMoleculeProperties(m)
    self.assertIsNone(mp)
    rdDistGeom.EmbedMultipleConfs(m, 2)
    res = ChemicalForceFields.MMFFOptimizeMoleculeConfs(m)
    self.assertEqual(len(res), 2)
    self.assertEqual(res[0], res[1])
    self.assertEqual(res[0], (-1, -1.0))

  def testOptimizeMolecule(self):
    m = Chem.AddHs(Chem.MolFromSmiles("CCCO"))
    self.assertIsNotNone(m)
    self.assertEqual(rdDistGeom.EmbedMolecule(m), 0)
    mp = ChemicalForceFields.MMFFGetMoleculeProperties(m)
    ff = ChemicalForceFields.MMFFGetMoleculeForceField(m, mp)
    before = ff.CalcEnergy()
    self.assertEqual(ChemicalForceFields.OptimizeMolecule(ff, maxIters=200), 0)
    after = ff.CalcEnergy()
    self.assertLess(after, before)

  def testOptimizeMoleculeConfs(self):
    m = Chem.AddHs(Chem.MolFromSmiles("CCCO"))
    self.assertIsNotNone(m)
    cids = rdDistGeom.EmbedMultipleConfs(m, numConfs=10)
    self.assertEqual(len(cids), 10)
    mp = ChemicalForceFields.MMFFGetMoleculeProperties(m)
    ff = ChemicalForceFields.MMFFGetMoleculeForceField(m, mp)
    before = [
      ChemicalForceFields.MMFFGetMoleculeForceField(m, mp, confId=cid).CalcEnergy() for cid in cids
    ]
    res, after = tuple(zip(*ChemicalForceFields.OptimizeMoleculeConfs(m, ff, maxIters=200)))
    self.assertEqual(len(res), 10)
    self.assertEqual(len(before), len(after))
    self.assertTrue(all(map(lambda i: i == 0, res)))
    self.assertTrue(all(after[i] < b for i, b in enumerate(before)))

  def testEmptyFF(self) -> None:
    m = Chem.MolFromSmiles(
      'CCCO |(-1.28533,-0.0567758,0.434662;-0.175447,0.695786,-0.299881;0.918409,-0.342619,-0.555572;1.30936,-0.801512,0.71705)|'
    )
    self.assertIsNotNone(m)
    ff = ChemicalForceFields.CreateEmptyForceFieldForMol(m)
    posa = m.GetConformer().GetAtomPosition(0)
    posb = m.GetConformer().GetAtomPosition(1)
    dist = (posa - posb).Length()
    self.assertTrue(ff)
    self.assertFalse(ff.Initialize())
    self.assertEqual(ff.CalcEnergy(), 0.0)
    self.assertTrue(all(v == 0.0 for v in ff.CalcGrad()))
    self.assertEqual(ff.NumPoints(), m.GetNumAtoms())
    self.assertEqual(len(ff.Positions()) / 3, m.GetNumAtoms())
    self.assertFalse(ff.Minimize())
    posa = m.GetConformer().GetAtomPosition(0)
    posb = m.GetConformer().GetAtomPosition(1)
    self.assertEqual((posa - posb).Length(), dist)

    ff.MMFFAddDistanceConstraint(0, 1, False, 100, 100, 100)
    self.assertFalse(ff.Minimize())
    pos = ff.Positions()
    dist = ((pos[0] - pos[3])**2 + (pos[1] - pos[4])**2 + (pos[2] - pos[5])**2)**0.5
    self.assertAlmostEqual(dist, 100, delta=10e-5)
    posa = m.GetConformer().GetAtomPosition(0)
    posb = m.GetConformer().GetAtomPosition(1)
    self.assertAlmostEqual((posa - posb).Length(), 100, delta=10e-5)


if __name__ == "__main__":
  unittest.main()
