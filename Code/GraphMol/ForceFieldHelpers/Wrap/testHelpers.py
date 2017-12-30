from rdkit import Chem
from rdkit.Chem import ChemicalForceFields
from rdkit import RDConfig
import unittest
import os
import numpy


def feq(v1, v2, tol2=1e-4):
  return abs(v1 - v2) <= tol2


class TestCase(unittest.TestCase):

  def setUp(self):
    self.dirName = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'ForceFieldHelpers', 'UFF',
                                'test_data')

  def test1(self):
    fName = os.path.join(self.dirName, 'benzene.mol')
    m = Chem.MolFromMolFile(fName)
    self.failIf(ChemicalForceFields.UFFOptimizeMolecule(m))
    # make sure that keyword arguments work:
    m = Chem.MolFromMolFile(fName)
    self.failUnless(ChemicalForceFields.UFFOptimizeMolecule(m, maxIters=1))

    m = Chem.MolFromMolFile(fName)
    self.failIf(ChemicalForceFields.UFFOptimizeMolecule(m, vdwThresh=2.0))

    m = Chem.MolFromMolFile(fName)
    self.failIf(ChemicalForceFields.UFFOptimizeMolecule(m, confId=-1))

    m = Chem.MolFromMolFile(fName)
    self.failUnlessRaises(ValueError, lambda: ChemicalForceFields.UFFOptimizeMolecule(m, confId=1))

  def test2(self):
    fName = os.path.join(self.dirName, 'benzene.mol')
    m = Chem.MolFromMolFile(fName)
    ff = ChemicalForceFields.UFFGetMoleculeForceField(m)
    self.failUnless(ff)
    e1 = ff.CalcEnergy()
    r = ff.Minimize()
    self.failUnless(r == 0)
    e2 = ff.CalcEnergy()
    self.failUnless(e2 < e1)

    # test keyword args:
    r = ff.Minimize(forceTol=1e-8)
    self.failUnless(r == 0)

    # test keyword args:
    r = ff.Minimize(energyTol=1e-3)
    self.failUnless(r == 0)

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
    self.failUnless(ff)
    e1 = ff.CalcEnergy()
    r = ff.Minimize()
    self.failUnless(r == 0)
    e2 = ff.CalcEnergy()
    self.failUnless(e2 < e1)

  def test4(self):
    m = Chem.MolFromSmiles('[Cu](C)(C)(C)(C)C')
    self.failIf(ChemicalForceFields.UFFHasAllMoleculeParams(m))

    m = Chem.MolFromSmiles('C(C)(C)(C)C')
    self.failUnless(ChemicalForceFields.UFFHasAllMoleculeParams(m))

  def test5(self):
    fName = os.path.join(self.dirName, 'benzene.mol')
    m = Chem.MolFromMolFile(fName)
    self.failIf(ChemicalForceFields.MMFFOptimizeMolecule(m))
    # make sure that keyword arguments work:
    m = Chem.MolFromMolFile(fName)
    self.failUnless(ChemicalForceFields.MMFFOptimizeMolecule(m, maxIters=1))

    m = Chem.MolFromMolFile(fName)
    self.failIf(ChemicalForceFields.MMFFOptimizeMolecule(m, nonBondedThresh=2.0))

    m = Chem.MolFromMolFile(fName)
    self.failIf(ChemicalForceFields.MMFFOptimizeMolecule(m, confId=-1))

    m = Chem.MolFromMolFile(fName)
    self.failUnlessRaises(ValueError, lambda: ChemicalForceFields.MMFFOptimizeMolecule(m, confId=1))

  def test6(self):
    fName = os.path.join(self.dirName, 'benzene.mol')
    m = Chem.MolFromMolFile(fName)
    mp = ChemicalForceFields.MMFFGetMoleculeProperties(m)
    ff = ChemicalForceFields.MMFFGetMoleculeForceField(m, mp)
    self.failUnless(ff)
    e1 = ff.CalcEnergy()
    r = ff.Minimize()
    self.failUnless(r == 0)
    e2 = ff.CalcEnergy()
    self.failUnless(e2 < e1)

    # test keyword args:
    r = ff.Minimize(forceTol=1.0e-8)
    self.failUnless(r == 0)

    # test keyword args:
    r = ff.Minimize(energyTol=1.0e-3)
    self.failUnless(r == 0)

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
    self.failUnless(ff)
    e1 = ff.CalcEnergy()
    r = ff.Minimize()
    self.failUnless(r == 0)
    e2 = ff.CalcEnergy()
    self.failUnless(e2 < e1)

  def test8(self):
    m = Chem.MolFromSmiles('[Cu](C)(C)(C)(C)C')
    self.failIf(ChemicalForceFields.MMFFHasAllMoleculeParams(m))

    m = Chem.MolFromSmiles('C(C)(C)(C)C')
    self.failUnless(ChemicalForceFields.MMFFHasAllMoleculeParams(m))

  def test9(self):
    m = Chem.MolFromSmiles('c1ccccc1CCNN')
    m = Chem.AddHs(m)
    mp = ChemicalForceFields.MMFFGetMoleculeProperties(m)
    mmffBondStretchParams = mp.GetMMFFBondStretchParams(m, 6, 7)
    self.failUnless(mmffBondStretchParams)
    self.failUnless((mmffBondStretchParams[0] == 0) and
                    (int(round(mmffBondStretchParams[1] * 1000) == 4258)) and
                    (int(round(mmffBondStretchParams[2] * 1000) == 1508)))
    mmffBondStretchParams = mp.GetMMFFBondStretchParams(m, 0, 7)
    self.failIf(mmffBondStretchParams)
    mmffAngleBendParams = mp.GetMMFFAngleBendParams(m, 6, 7, 8)
    self.failUnless(mmffAngleBendParams)
    self.failUnless((mmffAngleBendParams[0] == 0) and
                    (int(round(mmffAngleBendParams[1] * 1000) == 777)) and
                    (int(round(mmffAngleBendParams[2] * 1000) == 108290)))
    mmffAngleBendParams = mp.GetMMFFAngleBendParams(m, 0, 7, 8)
    self.failIf(mmffAngleBendParams)
    mmffStretchBendParams = mp.GetMMFFStretchBendParams(m, 6, 7, 8)
    self.failUnless(mmffStretchBendParams)
    self.failUnless((mmffStretchBendParams[0] == 0) and
                    (int(round(mmffStretchBendParams[1] * 1000) == 136)) and
                    (int(round(mmffStretchBendParams[2] * 1000) == 282)))
    mmffStretchBendParams = mp.GetMMFFStretchBendParams(m, 0, 7, 8)
    self.failIf(mmffStretchBendParams)
    mmffTorsionParams = mp.GetMMFFTorsionParams(m, 6, 7, 8, 9)
    self.failUnless(mmffTorsionParams)
    self.failUnless((mmffTorsionParams[0] == 0) and
                    (int(round(mmffTorsionParams[1] * 1000) == 0)) and
                    (int(round(mmffTorsionParams[2] * 1000) == -300)) and
                    (int(round(mmffTorsionParams[3] * 1000) == 500)))
    mmffTorsionParams = mp.GetMMFFTorsionParams(m, 0, 7, 8, 9)
    self.failIf(mmffTorsionParams)
    mmffOopBendParams = mp.GetMMFFOopBendParams(m, 6, 5, 4, 0)
    self.failUnless(mmffOopBendParams)
    self.failUnless(int(round(mmffOopBendParams * 1000)) == 40)
    mmffOopBendParams = mp.GetMMFFOopBendParams(m, 6, 5, 4, 1)
    self.failIf(mmffOopBendParams)
    sub1 = m.GetSubstructMatch(Chem.MolFromSmarts('NN[H]'))
    self.failUnless(len(sub1) == 3)
    nIdx = sub1[0]
    hIdx = sub1[2]
    mmffVdWParams = mp.GetMMFFVdWParams(nIdx, hIdx)
    self.failUnless(mmffVdWParams)
    self.failUnless((int(round(mmffVdWParams[0] * 1000)) == 3321) and
                    (int(round(mmffVdWParams[1] * 1000)) == 34) and
                    (int(round(mmffVdWParams[2] * 1000)) == 2657) and
                    (int(round(mmffVdWParams[3] * 1000)) == 17))

  def test10(self):
    m = Chem.MolFromSmiles('c1ccccc1CCNN')
    m = Chem.AddHs(m)
    uffBondStretchParams = ChemicalForceFields.GetUFFBondStretchParams(m, 6, 7)
    self.failUnless(uffBondStretchParams)
    self.failUnless((int(round(uffBondStretchParams[0] * 1000) == 699592)) and
                    (int(round(uffBondStretchParams[1] * 1000) == 1514)))
    uffBondStretchParams = ChemicalForceFields.GetUFFBondStretchParams(m, 0, 7)
    self.failIf(uffBondStretchParams)
    uffAngleBendParams = ChemicalForceFields.GetUFFAngleBendParams(m, 6, 7, 8)
    self.failUnless(uffAngleBendParams)
    self.failUnless((int(round(uffAngleBendParams[0] * 1000) == 303297)) and
                    (int(round(uffAngleBendParams[1] * 1000) == 109470)))
    uffAngleBendParams = ChemicalForceFields.GetUFFAngleBendParams(m, 0, 7, 8)
    self.failIf(uffAngleBendParams)
    uffTorsionParams = ChemicalForceFields.GetUFFTorsionParams(m, 6, 7, 8, 9)
    self.failUnless(uffTorsionParams)
    self.failUnless((int(round(uffTorsionParams * 1000) == 976)))
    uffTorsionParams = ChemicalForceFields.GetUFFTorsionParams(m, 0, 7, 8, 9)
    self.failIf(uffTorsionParams)
    uffInversionParams = ChemicalForceFields.GetUFFInversionParams(m, 6, 5, 4, 0)
    self.failUnless(uffInversionParams)
    self.failUnless(int(round(uffInversionParams * 1000)) == 2000)
    uffInversionParams = ChemicalForceFields.GetUFFInversionParams(m, 6, 5, 4, 1)
    self.failIf(uffInversionParams)
    uffVdWParams = ChemicalForceFields.GetUFFVdWParams(m, 0, 9)
    self.failUnless(uffVdWParams)
    self.failUnless((int(round(uffVdWParams[0] * 1000)) == 3754) and
                    (int(round(uffVdWParams[1] * 1000)) == 85))

  def test11(self):
    query = Chem.MolFromSmarts('c1cccn1')
    for i in [0, 1]:
      m = Chem.MolFromSmiles('Cc1nc(=O)c(C[NH3+])c(-c2c[nH]c3ccccc23)[nH]1')
      aromaticFlagsBefore = []
      for a in m.GetAtoms():
        aromaticFlagsBefore.append(a.GetIsAromatic())
      if (i):
        self.failUnless(ChemicalForceFields.MMFFGetMoleculeProperties(m))
      else:
        self.failUnless(ChemicalForceFields.MMFFHasAllMoleculeParams(m))
      aromaticFlagsAfter = []
      for a in m.GetAtoms():
        aromaticFlagsAfter.append(a.GetIsAromatic())
      res = (aromaticFlagsBefore == aromaticFlagsAfter)
      if (i):
        res = not res
      self.failUnless(res)
      pyrroleNIdx = m.GetSubstructMatch(query)[-1]
      self.failUnless(m.GetAtomWithIdx(pyrroleNIdx).GetTotalDegree() == 3)

  def test12(self):
    self.dirName = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'ForceFieldHelpers', 'UFF',
                                'test_data')
    fName = os.path.join(self.dirName, 'Issue239.mol')
    m = Chem.MolFromMolFile(fName)
    self.assertIsNotNone(m);
    ff = ChemicalForceFields.UFFGetMoleculeForceField(m)
    ff.Initialize()
    positions = ff.Positions()
    savedPos = list(positions)
    e1 = ff.CalcEnergy()
    ff.Minimize(10000, 1.0e-6, 1.0e-3)
    e2 = ff.CalcEnergy()
    self.failUnless(e2 < e1)
    e3 = ff.CalcEnergy(savedPos)
    self.assertAlmostEqual(e3, e1, 2);
    savedPos = tuple(positions)
    e3 = ff.CalcEnergy(savedPos)
    self.assertAlmostEqual(e3, e1, 2);
    savedPos = tuple(numpy.array(savedPos))
    e3 = ff.CalcEnergy(savedPos)
    self.assertAlmostEqual(e3, e1, 2);

  def test13(self):
    self.dirName = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'ForceFieldHelpers', 'UFF',
                                'test_data')
    fName = os.path.join(self.dirName, 'Issue239.mol')
    m = Chem.MolFromMolFile(fName)
    self.assertIsNotNone(m);
    mp = ChemicalForceFields.MMFFGetMoleculeProperties(m)
    ff = ChemicalForceFields.MMFFGetMoleculeForceField(m, mp)
    ff.Initialize()
    positions = ff.Positions()
    savedPos = list(positions)
    e1 = ff.CalcEnergy()
    ff.Minimize(10000, 1.0e-6, 1.0e-3)
    e2 = ff.CalcEnergy()
    self.failUnless(e2 < e1)
    e3 = ff.CalcEnergy(savedPos)
    self.assertAlmostEqual(e3, e1, 2);

  def test14(self):
    self.dirName = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'ForceFieldHelpers', 'UFF',
                                'test_data')
    fName = os.path.join(self.dirName, 'Issue239.mol')
    m = Chem.MolFromMolFile(fName)
    self.assertIsNotNone(m);
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
    self.assertIsNotNone(m);
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


if __name__ == '__main__':
  unittest.main()
