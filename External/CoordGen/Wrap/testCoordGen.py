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
  Mrv2014 03252113022D          

 14 15  0  0  0  0            999 V2000
   -0.7654   -0.7780    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5877   -0.8439    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0559   -0.1646    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7017    0.5805    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8794    0.6464    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4112   -0.0329    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4112    0.0329    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8794   -0.6464    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7771   -1.4237    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7017   -0.5805    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.0559    0.1646    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5877    0.8439    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7654    0.7780    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7654   -1.6519    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  7  1  0  0  0  0
  7  8  2  0  0  0  0
  8  9  1  0  0  0  0
  8 10  1  0  0  0  0
 10 11  2  0  0  0  0
 11 12  1  0  0  0  0
 12 13  2  0  0  0  0
  6  1  1  0  0  0  0
 13  7  1  0  0  0  0
  1 14  1  0  0  0  0
M  END
''')
    ref = Chem.MolFromMolBlock('''
     RDKit          2D

 14 15  0  0  0  0  0  0  0  0999 V2000
   -0.5957   -0.8221    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4185   -0.8441    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8520   -0.1434    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4595    0.5814    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6368    0.6072    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2063   -0.0955    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2087   -0.0625    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7102   -0.7158    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4040   -1.4774    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5259   -0.6099    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.8401    0.1516    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3405    0.8062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5252    0.6981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1735   -1.5240    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
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
M  END
''')
    ps = rdCoordGen.CoordGenParams()
    ps.minimizeOnly = True
    m2 = Chem.Mol(m)
    rdCoordGen.AddCoords(m2,ps)
    self.assertGreater(rdMolAlign.AlignMol(m, ref), 0.1)
    self.assertLess(rdMolAlign.AlignMol(m2, ref), 0.1)



if __name__ == '__main__':
  unittest.main()
