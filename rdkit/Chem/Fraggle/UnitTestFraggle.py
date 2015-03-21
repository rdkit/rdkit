# $Id$
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
from rdkit import RDConfig
import unittest,os
from rdkit.six.moves import cPickle
from rdkit import Chem
from rdkit.Chem.Fraggle import FraggleSim

class TestCase(unittest.TestCase):
  def testFragmentation(self):
    """ 
    
    """
    mol = Chem.MolFromSmiles('COc1cc(CN2CCC(CC2)NC(=O)c2cncc(C)c2)c(OC)c2ccccc12')
    frags = FraggleSim.generate_fraggle_fragmentation(mol)
    self.assertEqual(len(frags),16)

    expected=('[*]C(=O)NC1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1',
 '[*]C(=O)c1cncc(C)c1.[*]C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1',
 '[*]C(=O)c1cncc(C)c1.[*]c1cc(OC)c2ccccc2c1OC',
 '[*]C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1',
 '[*]C(=O)c1cncc(C)c1.[*]Cc1cc(OC)c2ccccc2c1OC',
 '[*]Cc1cc(OC)c2ccccc2c1OC.[*]NC(=O)c1cncc(C)c1',
 '[*]Cc1cc(OC)c2ccccc2c1OC.[*]c1cncc(C)c1',
 '[*]NC(=O)c1cncc(C)c1.[*]c1cc(OC)c2ccccc2c1OC',
 '[*]NC1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1',
 '[*]NC1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1.[*]c1cncc(C)c1',
 '[*]c1c(CN2CCC(NC(=O)c3cncc(C)c3)CC2)cc(OC)c2ccccc12',
 '[*]c1c(OC)cc(CN2CCC(NC(=O)c3cncc(C)c3)CC2)c(OC)c1[*]',
 '[*]c1cc(CN2CCC(NC(=O)c3cncc(C)c3)CC2)c(OC)c2ccccc12',
 '[*]N1CCC(NC(=O)c2cncc(C)c2)CC1.[*]c1cc(OC)c2ccccc2c1OC',
 '[*]C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1.[*]c1cncc(C)c1',
 '[*]c1cc(OC)c2ccccc2c1OC.[*]c1cncc(C)c1')
    for smi in frags:
        self.assertTrue(smi in expected)

if __name__ == '__main__':
  unittest.main()

