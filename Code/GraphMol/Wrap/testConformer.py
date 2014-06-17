# $Id$
#
#  Copyright (C) 2004  Rational Discovery LLC
#         All Rights Reserved
#
from rdkit import RDConfig
import os,sys
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

    self.assertTrue(cid == 0)
    
    conf2 = mol.GetConformer(0)
    self.assertTrue(conf2.GetId() == cid)
    pt1 = conf2.GetAtomPosition(0)
    self.assertTrue(ptEq(pt1, Point3D(-0.5, 0.0, 0.0)))
    
    pt2 = conf2.GetAtomPosition(1)
    self.assertTrue(ptEq(pt2, Point3D(1.0, 0.0, 0.0)))
    
    #changing conf should not change conf2 - related to issue 217
    conf.SetAtomPosition(1, Point3D(2.0, 0.0, 0.0))
    pt2 = conf2.GetAtomPosition(1)
    self.assertTrue(feq(pt2[0], 1.0))

    conf = Chem.Conformer(2)
    conf.SetAtomPosition(0, Point3D(-0.5, 0.0, 0.0))
    conf.SetAtomPosition(1, Point3D(1.0, 0.0, 0.0))
    conf.SetId(2)

    cid = mol.AddConformer(conf, 0)
    self.assertTrue(cid == 2)
    conf3 = mol.GetConformer(2)
    
  def test0AddHds(self) :
    mol = Chem.MolFromSmiles("CC")
    conf = Chem.Conformer(1)
    conf.SetAtomPosition(0, Point3D(-0.5, 0.0, 0.0))
    conf.SetAtomPosition(1, Point3D(1.0, 0.0, 0.0))
    cid = mol.AddConformer(conf)

    conf2 = mol.GetConformer()
    self.assertTrue(conf2.GetNumAtoms() == 2)

    nmol = Chem.AddHs(mol, 0,1)
    conf3 = nmol.GetConformer()
    self.assertTrue(conf3.GetNumAtoms() == 8)
    self.assertTrue(conf2.GetNumAtoms() == 2)

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
      self.assertTrue(ptEq(pt, Point3D(*tuple(targetCoords[i]))))

  def test2Issue217(self) :
    smi = 'c1ccccc1'
    m = Chem.MolFromSmiles(smi)
    addConf(m)
    self.assertTrue(m.GetNumConformers()==1);
    mb2 = Chem.MolToMolBlock(m)

  def test3Exceptions(self) :
    smi = 'c1ccccc1'
    m = Chem.MolFromSmiles(smi)
    addConf(m)
    self.assertTrue(m.GetNumConformers()==1)
    self.assertRaises(ValueError,lambda:m.GetConformer(2))

  def test4ConfTuple(self):
    smi = 'c1ccccc1'
    m = Chem.MolFromSmiles(smi)
    for i in range(10):
      addConf(m)

    confs = m.GetConformers()
    self.assertTrue(len(confs) == 10)

    for conf in confs:
      for i in range(6):
        pt = conf.GetAtomPosition(i)
        self.assertTrue(ptEq(pt, Point3D(0.0, 0.0, 0.0)))

    m.RemoveAllConformers()
    self.assertTrue(m.GetNumConformers() == 0)
    confs = m.GetConformers()
    self.assertTrue(confs == ())
    
if __name__ == '__main__':
  unittest.main()


