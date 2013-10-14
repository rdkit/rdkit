# $Id$
#
#  Copyright (C) 2004  Rational Discovery LLC
#         All Rights Reserved
#
from rdkit import RDConfig
import os,sys,math
import unittest
from rdkit import Chem
from rdkit.Geometry import Point3D

def feq(v1,v2,tol2=1e-4):
  return abs(v1-v2)<=tol2

def ptEq(pt1, pt2, tol=1e-4):
  return feq(pt1.x,pt2.x,tol) and feq(pt1.y,pt2.y,tol) and feq(pt1.z,pt2.z,tol)

def addConf(mol):
  conf = Chem.Conformer(mol.GetNumAtoms())
  for i in range(mol.GetNumAtoms()):
    conf.SetAtomPosition(i,(0.,0.,0.))
  mol.AddConformer(conf)
  mb = Chem.MolToMolBlock(mol)
  mb = Chem.MolToMolBlock(mol)
  
class TestCase(unittest.TestCase):
  def setUp(self):
    pass

  def test0Conformers(self) :
    """Test the conformer data structure"""
    mol = Chem.MolFromSmiles("CC")
    conf = Chem.Conformer(2)
    conf.SetAtomPosition(0, (-0.5, 0.0, 0.0))
    conf.SetAtomPosition(1, (1.0, 0.0, 0.0))
    conf.SetId(0)
    cid = mol.AddConformer(conf)

    self.failUnless(cid == 0)
    
    conf2 = mol.GetConformer(0)
    self.failUnless(conf2.GetId() == cid)
    pt1 = conf2.GetAtomPosition(0)
    self.failUnless(ptEq(pt1, Point3D(-0.5, 0.0, 0.0)))
    
    pt2 = conf2.GetAtomPosition(1)
    self.failUnless(ptEq(pt2, Point3D(1.0, 0.0, 0.0)))
    
    #changing conf should not change conf2 - related to issue 217
    conf.SetAtomPosition(1, Point3D(2.0, 0.0, 0.0))
    pt2 = conf2.GetAtomPosition(1)
    self.failUnless(feq(pt2[0], 1.0))

    conf = Chem.Conformer(2)
    conf.SetAtomPosition(0, Point3D(-0.5, 0.0, 0.0))
    conf.SetAtomPosition(1, Point3D(1.0, 0.0, 0.0))
    conf.SetId(2)

    cid = mol.AddConformer(conf, 0)
    self.failUnless(cid == 2)
    conf3 = mol.GetConformer(2)
    
  def test0AddHds(self) :
    mol = Chem.MolFromSmiles("CC")
    conf = Chem.Conformer(1)
    conf.SetAtomPosition(0, Point3D(-0.5, 0.0, 0.0))
    conf.SetAtomPosition(1, Point3D(1.0, 0.0, 0.0))
    cid = mol.AddConformer(conf)

    conf2 = mol.GetConformer()
    self.failUnless(conf2.GetNumAtoms() == 2)

    nmol = Chem.AddHs(mol, 0,1)
    conf3 = nmol.GetConformer()
    self.failUnless(conf3.GetNumAtoms() == 8)
    self.failUnless(conf2.GetNumAtoms() == 2)

    targetCoords = [[-0.5, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [-0.8667, 0.0, 1.03709],
                    [-0.8667, 0.8981, -0.5185],
                    [-0.8667, -0.8981, -0.5185],
                    [1.3667, 0.0, -1.0371],
                    [1.36667, 0.8981, 0.5185],
                    [1.36667, -0.8981, 0.5185]]

    for i in range(8) :
      pt = conf3.GetAtomPosition(i)
      self.failUnless(ptEq(pt, apply(Point3D,tuple(targetCoords[i]))))

  def test2Issue217(self) :
    smi = 'c1ccccc1'
    m = Chem.MolFromSmiles(smi)
    addConf(m)
    self.failUnless(m.GetNumConformers()==1);
    mb2 = Chem.MolToMolBlock(m)

  def test3Exceptions(self) :
    smi = 'c1ccccc1'
    m = Chem.MolFromSmiles(smi)
    addConf(m)
    self.failUnless(m.GetNumConformers()==1)
    self.failUnlessRaises(ValueError,lambda:m.GetConformer(2))

  def test4ConfTuple(self):
    smi = 'c1ccccc1'
    m = Chem.MolFromSmiles(smi)
    for i in range(10):
      addConf(m)

    confs = m.GetConformers()
    self.failUnless(len(confs) == 10)

    for conf in confs:
      for i in range(6):
        pt = conf.GetAtomPosition(i)
        self.failUnless(ptEq(pt, Point3D(0.0, 0.0, 0.0)))

    m.RemoveAllConformers()
    self.failUnless(m.GetNumConformers() == 0)
    confs = m.GetConformers()
    self.failUnless(confs == ())

  def testGetSetBondLength(self):
    file = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
                         'test_data', '3-cyclohexylpyridine.mol')

    m = Chem.MolFromMolFile(file, True, False)
    conf = m.GetConformer()
    dist = conf.GetBondLength(0, 19)
    self.failUnlessAlmostEqual(dist, 1.36, 2)
    conf.SetBondLength(0, 19, 2.5)
    dist = conf.GetBondLength(0, 19)
    self.failUnlessAlmostEqual(dist, 2.5, 1)
    conf.SetBondLength(19, 0, 3.0)
    dist = conf.GetBondLength(0, 19)
    self.failUnlessAlmostEqual(dist, 3.0, 1)

  def testGetSetAngle(self):
    file = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
                         'test_data', '3-cyclohexylpyridine.mol')

    m = Chem.MolFromMolFile(file, True, False)
    conf = m.GetConformer()
    angle = conf.GetAngleDeg(0, 19, 21)
    self.failUnlessAlmostEqual(angle, 109.7, 1)
    conf.SetAngleDeg(0, 19, 21, 125.0)
    angle = conf.GetAngleDeg(0, 19, 21);
    self.failUnlessAlmostEqual(angle, 125.0, 1)
    conf.SetAngleRad(21, 19, 0, math.pi / 2.)
    angle = conf.GetAngleRad(0, 19, 21)
    self.failUnlessAlmostEqual(angle, math.pi / 2., 1)
    angle = conf.GetAngleDeg(0, 19, 21)
    self.failUnlessAlmostEqual(angle, 90.0, 1)

  def testGetSetDihedral(self):
    file = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
                         'test_data', '3-cyclohexylpyridine.mol')

    m = Chem.MolFromMolFile(file, True, False)
    conf = m.GetConformer()
    dihedral = conf.GetDihedralDeg(0, 19, 21, 24)
    self.failUnlessAlmostEqual(dihedral, 176.05, 2)
    conf.SetDihedralDeg(8, 0, 19, 21, 65.0)
    dihedral = conf.GetDihedralDeg(8, 0, 19, 21)
    self.failUnlessAlmostEqual(dihedral, 65.0, 1)
    conf.SetDihedralDeg(8, 0, 19, 21, -130.0)
    dihedral = conf.GetDihedralDeg(8, 0, 19, 21)
    self.failUnlessAlmostEqual(dihedral, -130.0, 1)
    conf.SetDihedralRad(21, 19, 0, 8, -2. / 3. * math.pi)
    dihedral = conf.GetDihedralRad(8, 0, 19, 21)
    self.failUnlessAlmostEqual(dihedral, -2. / 3. * math.pi, 1)
    dihedral = conf.GetDihedralDeg(8, 0, 19, 21)
    self.failUnlessAlmostEqual(dihedral, -120.0, 1)


if __name__ == '__main__':
  unittest.main()


