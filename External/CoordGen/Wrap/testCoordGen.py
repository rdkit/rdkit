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
from rdkit import Chem

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


if __name__ == '__main__':
  unittest.main()
