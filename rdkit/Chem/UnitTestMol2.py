# $Id$
#
#  Copyright (C) 2001-2006  greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""basic unit testing code for the molecule boost wrapper

"""
from rdkit import RDConfig
import unittest,cPickle,os
from rdkit import Chem
import tempfile, os

class TestCase(unittest.TestCase):
  def setUp(self):
    self._files=[]
  def tearDown(self):
    for fileN in self._files:
      try:
        os.unlink(fileN)
      except OSError:
        pass

  def testMol2(self, removeHs=False):
    self.testMol2String(True)
    self.testMol2String(False)
    self.testMol2File(True)
    self.testMol2File(False)
    
  def testMol2String(self, removeHs=False): 
    " testing 5k molecule pickles "
    mol = Chem.MolFromMol2File('{}/CSAR_2014_01_FXA_gtc101/Decoys.mol2'.format(RDConfig.RDDataDir), removeHs=removeHs)
    mol_block = Chem.MolToMol2Block(mol)
    mol2 = Chem.MolFromMol2Block(mol_block, removeHs=removeHs)
    self.assertEqual(mol.GetNumAtoms(), mol2.GetNumAtoms())
    self.assertEqual(mol.GetNumBonds(), mol2.GetNumBonds())

  def testMol2File(self, removeHs=False): 
    " testing 5k molecule pickles "
    import tempfile
    fileN, filename = tempfile.mkstemp(suffix='.mol2')
    
    mol = Chem.MolFromMol2File('{}/CSAR_2014_01_FXA_gtc101/Decoys.mol2'.format(RDConfig.RDDataDir), removeHs=removeHs)
    Chem.MolToMol2File(mol, filename)
    mol2 = Chem.MolFromMol2File(filename, removeHs=removeHs)
    self.assertEqual(mol.GetNumAtoms(), mol2.GetNumAtoms())
    self.assertEqual(mol.GetNumBonds(), mol2.GetNumBonds())    


if __name__ == '__main__':
  unittest.main()


