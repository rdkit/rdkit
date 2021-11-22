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

if __name__ == '__main__':
  unittest.main()
