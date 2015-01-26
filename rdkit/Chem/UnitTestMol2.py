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
import gzip

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
    self._testMol2String(True)
    self._testMol2String(False)
    self._testMol2File(True)
    self._testMol2File(False)
    pass
    
  def _testMol2String(self, removeHs=False): 
    " testing 5k molecule pickles "
    for mol1 in Chem.ForwardSDMolSupplier(gzip.open('{}/Regress/Data/mols.1000.sdf.gz'.format(RDConfig.RDBaseDir))):
        if mol1:
            mol2 = Chem.MolFromMol2Block(Chem.MolToMol2Block(mol1))
            if mol2:
                self.assertEqual(mol1.GetNumAtoms(), mol2.GetNumAtoms())
                self.assertEqual(Chem.MolToSmiles(mol1), Chem.MolToSmiles(mol2))
            else:
                print 'Could not read molecule back from mol2'

  def _testMol2File(self, removeHs=False): 
    " testing 5k molecule pickles "
    import tempfile
    fileN, filename = tempfile.mkstemp(suffix='.mol2')
    
    for mol1 in Chem.ForwardSDMolSupplier(gzip.open('{}/Regress/Data/mols.1000.sdf.gz'.format(RDConfig.RDBaseDir))):
        if mol1:
            Chem.MolToMol2File(mol1, filename)
            mol2 = Chem.MolFromMol2File(filename)
            if mol2:
                self.assertEqual(mol1.GetNumAtoms(), mol2.GetNumAtoms())
                self.assertEqual(Chem.MolToSmiles(mol1), Chem.MolToSmiles(mol2))
            else:
                print 'Could not read molecule back from mol2'

if __name__ == '__main__':
  unittest.main()


