#
#  Copyright (C) 2024 Tad Hurst
#         All Rights Reserved
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

#
import os
import sys
import unittest

from rdkit import Chem
from rdkit.Chem import RDConfig


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def testScsr(self):
    """Test the Scsr system"""

    ofile = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data','macromols',
                         'Triplet.mol')
    with open(ofile) as inf:
      scsrBlock = inf.read()

   
    scsr = Chem.ScsrFromScsrBlock(scsrBlock,False,False)    
    self.assertTrue(scsr.GetNumTemplates() == 3)
    t1 = scsr.GetTemplate(0)
    self.assertTrue(t1.GetNumAtoms() == 11)
    self.assertTrue(scsr.GetMol().GetNumAtoms() == 3)

    molFromScsrParams = Chem.MolFromScsrParams()
    molFromScsrParams.includeLeavingGroups = True
    
    mol = Chem.ScsrToMol(scsr, molFromScsrParams)

    self.assertTrue(mol.GetNumAtoms() == 30)
    sgs = Chem.GetMolSubstanceGroups(mol)
    self.assertTrue(len(sgs) == 3)

if __name__ == '__main__':
  unittest.main()
