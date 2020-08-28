#  Copyright (c) 2017 Greg Landrum
#  All rights reserved.
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.


import unittest
import os,sys, copy

from rdkit.Chem import rdCoordGen, rdMolAlign
from rdkit import Chem, Geometry

def compareConfs(c1,c2,match, tol=1e-2, alignIt=False):
  for i,j in enumerate(match):
    pi = c2.GetAtomPosition(i)
    pj = c1.GetAtomPosition(j)
    if (pj-pi).Length()>=tol:
      return False
  return True

class TestCase(unittest.TestCase) :
  def test_basics(self):
    mol = Chem.MolFromSmiles('CCOC')
    self.assertEqual(mol.GetNumConformers(),0)
    rdCoordGen.AddCoords(mol)
    self.assertEqual(mol.GetNumConformers(),1)
    rdCoordGen.AddCoords(mol)
    self.assertEqual(mol.GetNumConformers(),1)
    rwmol = Chem.RWMol(mol)
    rdCoordGen.AddCoords(rwmol)
    self.assertEqual(rwmol.GetNumConformers(),1)
  def test_template1(self):
    template = Chem.MolFromSmiles('C1OCOCOCOCCCNCCC1')
    template.SetProp("_Name",'template')
    mol = Chem.MolFromSmiles('C1OCOCOCOCCCNC(OC(=O)C2CC2)CC1')
    mol.SetProp("_Name","mol")
    rdCoordGen.AddCoords(template)
    rdCoordGen.AddCoords(mol)
    self.assertFalse(compareConfs(mol.GetConformer(),template.GetConformer(),mol.GetSubstructMatch(template)))
    match = mol.GetSubstructMatch(template)
    mapd = dict()
    for i,aid in enumerate(match):
      p = template.GetConformer().GetAtomPosition(i)
      mapd[aid] = Geometry.Point2D(p.x,p.y)
    ps = rdCoordGen.CoordGenParams()
    ps.SetCoordMap(mapd)
    ps.dbg_useFixed = True
    rdCoordGen.AddCoords(mol,ps)
    self.assertTrue(compareConfs(mol.GetConformer(),template.GetConformer(),mol.GetSubstructMatch(template)))

  def test_template2(self):
    # the easier way...
    template = Chem.MolFromSmiles('C1OCOCOCOCCCNCCC1')
    template.SetProp("_Name",'template')
    mol = Chem.MolFromSmiles('C1OCOCOCOCCCNC(OC(=O)C2CC2)CC1')
    mol.SetProp("_Name","mol")
    rdCoordGen.AddCoords(template)
    ps = rdCoordGen.CoordGenParams()
    ps.SetTemplateMol(template)
    ps.dbg_useFixed = True
    rdCoordGen.AddCoords(mol,ps)
    self.assertTrue(compareConfs(mol.GetConformer(),template.GetConformer(),mol.GetSubstructMatch(template)))

  def test_template3(self):
    # the easier way, test lifetime...
    template = Chem.MolFromSmiles('C1OCOCOCOCCCNCCC1')
    mol = Chem.MolFromSmiles('C1OCOCOCOCCCNC(OC(=O)C2CC2)CC1')
    rdCoordGen.AddCoords(template)
    template2 = Chem.Mol(template)
    ps = rdCoordGen.CoordGenParams()
    ps.SetTemplateMol(template2)
    template2 = None
    ps.dbg_useFixed = True
    rdCoordGen.AddCoords(mol,ps)
    self.assertTrue(compareConfs(mol.GetConformer(),template.GetConformer(),mol.GetSubstructMatch(template)))

  def testMinimizeOnly(self):
    m = Chem.MolFromMolBlock('''
  Mrv2014 08052005392D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 14 15 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.4287 -1.4523 0 0
M  V30 2 C -2.9638 -1.5752 0 0
M  V30 3 C -3.8377 -0.3072 0 0
M  V30 4 C -3.1766 1.0837 0 0
M  V30 5 C -1.6416 1.2066 0 0
M  V30 6 C -0.7675 -0.0614 0 0
M  V30 7 C 0.7675 0.0614 0 0
M  V30 8 C 1.6416 -1.2066 0 0
M  V30 9 C 0.9804 -2.5975 0 0
M  V30 10 N 3.1766 -1.0837 0 0
M  V30 11 C 3.8377 0.3072 0 0
M  V30 12 C 2.9638 1.5752 0 0
M  V30 13 C 1.4287 1.4523 0 0
M  V30 14 F -0.5548 -2.7203 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 7
M  V30 7 2 7 8
M  V30 8 1 8 9
M  V30 9 1 8 10
M  V30 10 2 10 11
M  V30 11 1 11 12
M  V30 12 2 12 13
M  V30 13 1 6 1
M  V30 14 1 13 7
M  V30 15 1 1 14
M  V30 END BOND
M  V30 END CTAB
M  END
''')
    ref = Chem.MolFromMolBlock('''
     RDKit          2D

 14 15  0  0  0  0  0  0  0  0999 V2000
   -1.5379   -1.4859    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1218   -1.5958    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9595   -0.2554    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2641    1.1663    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6778    1.2217    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7941   -0.0886    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7983    0.0383    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7524   -1.2209    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1246   -2.6443    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.3306   -1.0787    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    3.9439    0.3747    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0337    1.6656    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4617    1.4692    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6924   -2.7917    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  2  0
  8  9  1  0
  8 10  1  0
 10 11  2  0
 11 12  1  0
 12 13  2  0
  6  1  1  0
 13  7  1  0
  1 14  1  0
M  END''')
    ps = rdCoordGen.CoordGenParams()
    ps.minimizeOnly = True
    m2 = Chem.Mol(m)
    rdCoordGen.AddCoords(m2,ps)
    self.assertGreater(rdMolAlign.AlignMol(m, ref), 0.1)
    self.assertLess(rdMolAlign.AlignMol(m2, ref), 0.1)



if __name__ == '__main__':
  unittest.main()
