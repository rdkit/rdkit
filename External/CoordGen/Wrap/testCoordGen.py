#  Copyright (c) 2017 Greg Landrum
#  All rights reserved.
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
from __future__ import print_function

import unittest
import os,sys, copy

from rdkit.Chem import rdCoordGen
from rdkit import Chem, Geometry

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
    def test_template(self):
        template = Chem.MolFromSmiles('C1OCOCOCOCCCNCCC1')
        template.SetProp("_Name",'template')
        mol = Chem.MolFromSmiles('C1OCOCOCOCCCNC(OC(=O)C2CC2)CC1')
        mol.SetProp("_Name","mol")
        rdCoordGen.AddCoords(template)
        rdCoordGen.AddCoords(mol)
        print(Chem.MolToMolBlock(template))
        print(Chem.MolToMolBlock(mol))
        match = mol.GetSubstructMatch(template)
        mapd = dict()
        for i,aid in enumerate(match):
            p = template.GetConformer().GetAtomPosition(i)
            mapd[aid] = Geometry.Point2D(p.x,p.y)
        ps = rdCoordGen.CoordGenParams()
        ps.SetCoordMap(mapd)
        rdCoordGen.AddCoords(mol,ps)
        mol.SetProp("_Name","templated")
        print(Chem.MolToMolBlock(mol))



if __name__ == '__main__':
  unittest.main()
