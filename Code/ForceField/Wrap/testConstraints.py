from rdkit import RDConfig
import os
import unittest
from rdkit import Chem
from rdkit.Chem import ChemicalForceFields
from rdkit.Chem import rdMolTransforms

class TestCase(unittest.TestCase):
  def setUp(self) :
    self.molB = """butane
     RDKit          3D
butane
 17 16  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4280    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7913   -0.2660    0.9927 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.9040    1.3004   -0.3485 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5407    2.0271    0.3782 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.5407    1.5664   -1.3411 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.3320    1.3004   -0.3485 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6953    1.5162   -1.3532 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.8080    0.0192    0.0649 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4447   -0.7431   -0.6243 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.4447   -0.1966    1.0697 H   0  0  0  0  0  0  0  0  0  0  0  0
    4.8980    0.0192    0.0649 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.6954    2.0627    0.3408 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.7913   -0.7267   -0.7267 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3633    0.7267    0.7267 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3633   -0.9926    0.2660 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3633    0.2660   -0.9926 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1 15  1  0  0  0  0
  1 16  1  0  0  0  0
  1 17  1  0  0  0  0
  2  3  1  0  0  0  0
  2  4  1  0  0  0  0
  2 14  1  0  0  0  0
  4  5  1  0  0  0  0
  4  6  1  0  0  0  0
  4  7  1  0  0  0  0
  7  8  1  0  0  0  0
  7  9  1  0  0  0  0
  7 13  1  0  0  0  0
  9 10  1  0  0  0  0
  9 11  1  0  0  0  0
  9 12  1  0  0  0  0
M  END"""

  def testUFFDistanceConstraints(self) :
    m = Chem.MolFromMolBlock(self.molB, True, False)
    ff = ChemicalForceFields.UFFGetMoleculeForceField(m)
    self.assertTrue(ff)
    ff.UFFAddDistanceConstraint(1, 3, False, 2.0, 2.0, 1.0e5)
    r = ff.Minimize()
    self.assertTrue(r == 0)
    conf = m.GetConformer()
    dist = rdMolTransforms.GetBondLength(conf, 1, 3)
    self.assertTrue(dist > 1.99)
    ff = ChemicalForceFields.UFFGetMoleculeForceField(m)
    self.assertTrue(ff)
    ff.UFFAddDistanceConstraint(1, 3, True, -0.2, 0.2, 1.0e5)
    r = ff.Minimize()
    self.assertTrue(r == 0)
    conf = m.GetConformer()
    dist = rdMolTransforms.GetBondLength(conf, 1, 3)
    self.assertTrue(dist > 1.79)

  def testUFFAngleConstraints(self) :
    m = Chem.MolFromMolBlock(self.molB, True, False)
    ff = ChemicalForceFields.UFFGetMoleculeForceField(m)
    self.assertTrue(ff)
    ff.UFFAddAngleConstraint(1, 3, 6, False, 90.0, 90.0, 1.0e5)
    r = ff.Minimize()
    self.assertTrue(r == 0)
    conf = m.GetConformer()
    angle = rdMolTransforms.GetAngleDeg(conf, 1, 3, 6)
    self.assertTrue(int(angle) == 90)
    ff = ChemicalForceFields.UFFGetMoleculeForceField(m)
    self.assertTrue(ff)
    ff.UFFAddAngleConstraint(1, 3, 6, True, -10.0, 10.0, 1.0e5)
    r = ff.Minimize()
    self.assertTrue(r == 0)
    conf = m.GetConformer()
    angle = rdMolTransforms.GetAngleDeg(conf, 1, 3, 6)
    self.assertTrue(int(angle) == 100)
    ff = ChemicalForceFields.UFFGetMoleculeForceField(m)
    self.assertTrue(ff)
    ff.UFFAddAngleConstraint(1, 3, 6, False, -10.0, 10.0, 1.0e5)
    r = ff.Minimize()
    self.assertTrue(r == 0)
    conf = m.GetConformer()
    angle = rdMolTransforms.GetAngleDeg(conf, 1, 3, 6)
    self.assertTrue(int(angle) == 10)

  def testUFFTorsionConstraints(self) :
    m = Chem.MolFromMolBlock(self.molB, True, False)
    conf = m.GetConformer()
    rdMolTransforms.SetDihedralDeg(conf, 1, 3, 6, 8, 50.0);
    ff = ChemicalForceFields.UFFGetMoleculeForceField(m)
    self.assertTrue(ff)
    ff.UFFAddTorsionConstraint(1, 3, 6, 8, False, 10.0, 20.0, 1.0e5)
    r = ff.Minimize()
    self.assertTrue(r == 0)
    dihedral = rdMolTransforms.GetDihedralDeg(conf, 1, 3, 6, 8)
    self.assertTrue(int(dihedral) == 20)
    rdMolTransforms.SetDihedralDeg(conf, 1, 3, 6, 8, -30.0);
    ff = ChemicalForceFields.UFFGetMoleculeForceField(m)
    self.assertTrue(ff)
    ff.UFFAddTorsionConstraint(1, 3, 6, 8, True, -10.0, -8.0, 1.0e5)
    r = ff.Minimize()
    self.assertTrue(r == 0)
    conf = m.GetConformer()
    dihedral = rdMolTransforms.GetDihedralDeg(conf, 1, 3, 6, 8)
    self.assertTrue(int(dihedral) == -40)
    rdMolTransforms.SetDihedralDeg(conf, 1, 3, 6, 8, 30.0);
    ff = ChemicalForceFields.UFFGetMoleculeForceField(m)
    self.assertTrue(ff)
    ff.UFFAddTorsionConstraint(1, 3, 6, 8, False, -10.0, -8.0, 1.0e5)
    r = ff.Minimize(1000)
    self.assertTrue(r == 0)
    conf = m.GetConformer()
    dihedral = rdMolTransforms.GetDihedralDeg(conf, 1, 3, 6, 8)
    self.assertTrue(int(dihedral) == -10)

  def testUFFPositionConstraints(self) :
    m = Chem.MolFromMolBlock(self.molB, True, False)
    ff = ChemicalForceFields.UFFGetMoleculeForceField(m)
    self.assertTrue(ff)
    conf = m.GetConformer()
    p = conf.GetAtomPosition(1)
    ff.UFFAddPositionConstraint(1, 0.3, 1.0e5)
    r = ff.Minimize()
    self.assertTrue(r == 0)
    q = conf.GetAtomPosition(1)
    self.assertTrue((p - q).Length() < 0.3)

  def testUFFFixedAtoms(self) :
    m = Chem.MolFromMolBlock(self.molB, True, False)
    ff = ChemicalForceFields.UFFGetMoleculeForceField(m)
    self.assertTrue(ff)
    conf = m.GetConformer()
    fp = conf.GetAtomPosition(1)
    ff.AddFixedPoint(1)
    r = ff.Minimize()
    self.assertTrue(r == 0)
    fq = conf.GetAtomPosition(1)
    self.assertTrue((fp - fq).Length() < 0.01)

  def testMMFFDistanceConstraints(self) :
    m = Chem.MolFromMolBlock(self.molB, True, False)
    mp = ChemicalForceFields.MMFFGetMoleculeProperties(m)
    ff = ChemicalForceFields.MMFFGetMoleculeForceField(m, mp)
    self.assertTrue(ff)
    ff.MMFFAddDistanceConstraint(1, 3, False, 2.0, 2.0, 1.0e5)
    r = ff.Minimize()
    self.assertTrue(r == 0)
    conf = m.GetConformer()
    dist = rdMolTransforms.GetBondLength(conf, 1, 3)
    self.assertTrue(dist > 1.99)
    ff = ChemicalForceFields.MMFFGetMoleculeForceField(m, mp)
    self.assertTrue(ff)
    ff.MMFFAddDistanceConstraint(1, 3, True, -0.2, 0.2, 1.0e5)
    r = ff.Minimize()
    self.assertTrue(r == 0)
    conf = m.GetConformer()
    dist = rdMolTransforms.GetBondLength(conf, 1, 3)
    self.assertTrue(dist > 1.79)

  def testMMFFAngleConstraints(self) :
    m = Chem.MolFromMolBlock(self.molB, True, False)
    mp = ChemicalForceFields.MMFFGetMoleculeProperties(m)
    ff = ChemicalForceFields.MMFFGetMoleculeForceField(m, mp)
    self.assertTrue(ff)
    ff.MMFFAddAngleConstraint(1, 3, 6, False, 90.0, 90.0, 1.0e5)
    r = ff.Minimize()
    self.assertTrue(r == 0)
    conf = m.GetConformer()
    angle = rdMolTransforms.GetAngleDeg(conf, 1, 3, 6)
    self.assertTrue(int(angle) == 90)
    ff = ChemicalForceFields.MMFFGetMoleculeForceField(m, mp)
    self.assertTrue(ff)
    ff.MMFFAddAngleConstraint(1, 3, 6, True, -10.0, 10.0, 1.0e5)
    r = ff.Minimize()
    self.assertTrue(r == 0)
    conf = m.GetConformer()
    angle = rdMolTransforms.GetAngleDeg(conf, 1, 3, 6)
    self.assertTrue(int(angle) == 100)
    ff = ChemicalForceFields.MMFFGetMoleculeForceField(m, mp)
    self.assertTrue(ff)
    ff.MMFFAddAngleConstraint(1, 3, 6, False, -10.0, 10.0, 1.0e5)
    r = ff.Minimize()
    self.assertTrue(r == 0)
    conf = m.GetConformer()
    angle = rdMolTransforms.GetAngleDeg(conf, 1, 3, 6)
    self.assertTrue(int(angle) == 10)

  def testMMFFTorsionConstraints(self) :
    m = Chem.MolFromMolBlock(self.molB, True, False)
    conf = m.GetConformer()
    rdMolTransforms.SetDihedralDeg(conf, 1, 3, 6, 8, 50.0);
    mp = ChemicalForceFields.MMFFGetMoleculeProperties(m)
    ff = ChemicalForceFields.MMFFGetMoleculeForceField(m, mp)
    self.assertTrue(ff)
    ff.MMFFAddTorsionConstraint(1, 3, 6, 8, False, 10.0, 20.0, 1.0e5)
    r = ff.Minimize()
    self.assertTrue(r == 0)
    dihedral = rdMolTransforms.GetDihedralDeg(conf, 1, 3, 6, 8)
    self.assertTrue(int(dihedral) == 20)
    rdMolTransforms.SetDihedralDeg(conf, 1, 3, 6, 8, -30.0);
    ff = ChemicalForceFields.MMFFGetMoleculeForceField(m, mp)
    self.assertTrue(ff)
    ff.MMFFAddTorsionConstraint(1, 3, 6, 8, True, -10.0, -8.0, 1.0e5)
    r = ff.Minimize()
    self.assertTrue(r == 0)
    conf = m.GetConformer()
    dihedral = rdMolTransforms.GetDihedralDeg(conf, 1, 3, 6, 8)
    self.assertTrue(int(dihedral) == -40)
    rdMolTransforms.SetDihedralDeg(conf, 1, 3, 6, 8, 30.0);
    ff = ChemicalForceFields.MMFFGetMoleculeForceField(m, mp)
    self.assertTrue(ff)
    ff.MMFFAddTorsionConstraint(1, 3, 6, 8, False, -10.0, -8.0, 1.0e5)
    r = ff.Minimize(1000)
    self.assertTrue(r == 0)
    conf = m.GetConformer()
    dihedral = rdMolTransforms.GetDihedralDeg(conf, 1, 3, 6, 8)
    self.assertTrue(int(dihedral) == -10)

  def testMMFFPositionConstraints(self) :
    m = Chem.MolFromMolBlock(self.molB, True, False)
    mp = ChemicalForceFields.MMFFGetMoleculeProperties(m)
    ff = ChemicalForceFields.MMFFGetMoleculeForceField(m, mp)
    self.assertTrue(ff)
    conf = m.GetConformer()
    p = conf.GetAtomPosition(1)
    ff.MMFFAddPositionConstraint(1, 0.3, 1.0e5)
    r = ff.Minimize()
    self.assertTrue(r == 0)
    q = conf.GetAtomPosition(1)
    self.assertTrue((p - q).Length() < 0.3)

  def testMMFFFixedAtoms(self) :
    m = Chem.MolFromMolBlock(self.molB, True, False)
    mp = ChemicalForceFields.MMFFGetMoleculeProperties(m)
    ff = ChemicalForceFields.MMFFGetMoleculeForceField(m, mp)
    self.assertTrue(ff)
    conf = m.GetConformer()
    fp = conf.GetAtomPosition(1)
    ff.AddFixedPoint(1)
    r = ff.Minimize()
    self.assertTrue(r == 0)
    fq = conf.GetAtomPosition(1)
    self.assertTrue((fp - fq).Length() < 0.01)





if __name__== '__main__':
    unittest.main()
