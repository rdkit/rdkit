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

  def testSCSR(self):
    """Test the SCSR system"""

    ofile = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'macromols', 'Triplet.mol')
    with open(ofile) as inf:
      scsrBlock = inf.read()

    scsr = Chem.SCSRFromSCSRBlock(scsrBlock, False, False)
    self.assertTrue(scsr.GetNumTemplates() == 3)
    t1 = scsr.GetTemplate(0)
    self.assertTrue(t1.GetNumAtoms() == 11)
    self.assertTrue(scsr.GetMol().GetNumAtoms() == 3)

    molFromSCSRParams = Chem.MolFromSCSRParams()
    molFromSCSRParams.includeLeavingGroups = True
    molFromSCSRParams.scsrTemplateNames = Chem.SCSRTemplateNamesVal.SCSRTemplateNamesAsEntered
    molFromSCSRParams.scsrBaseHbondOptions = Chem.SCSRBaseHbondOptionsVal.SCSRBaseHbondOptionsIgnore

    mol = Chem.SCSRToMol(scsr, molFromSCSRParams)

    self.assertTrue(mol.GetNumAtoms() == 30)
    sgs = Chem.GetMolSubstanceGroups(mol)
    self.assertTrue(len(sgs) == 6)

  def testSCSRRna(self):
    """Test the SCSR system with and RNA double strand"""

    ofile = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'macromols', 'DnaTest.mol')
    with open(ofile) as inf:
      scsrBlock = inf.read()

    scsr = Chem.SCSRFromSCSRBlock(scsrBlock, False, False)
    self.assertTrue(scsr.GetNumTemplates() == 7)
    t1 = scsr.GetTemplate(0)
    self.assertTrue(t1.GetNumAtoms() == 12)
    self.assertTrue(scsr.GetMol().GetNumAtoms() == 34)
    self.assertTrue(scsr.GetMol().GetNumBonds() == 38)

    molFromSCSRParams = Chem.MolFromSCSRParams()
    molFromSCSRParams.includeLeavingGroups = True
    molFromSCSRParams.scsrTemplateNames = Chem.SCSRTemplateNamesVal.SCSRTemplateNamesAsEntered
    molFromSCSRParams.scsrBaseHbondOptions = Chem.SCSRBaseHbondOptionsVal.SCSRBaseHbondOptionsAuto

    mol = Chem.SCSRToMol(scsr, molFromSCSRParams)

    self.assertTrue(mol.GetNumAtoms() == 254)
    self.assertTrue(mol.GetNumBonds() == 300)
    sgs = Chem.GetMolSubstanceGroups(mol)
    self.assertTrue(len(sgs) == 38)

    molFromSCSRParams.scsrBaseHbondOptions = Chem.SCSRBaseHbondOptionsVal.SCSRBaseHbondOptionsIgnore
    mol = Chem.SCSRToMol(scsr, molFromSCSRParams)

    self.assertTrue(mol.GetNumAtoms() == 254)
    self.assertTrue(mol.GetNumBonds() == 282)
    sgs = Chem.GetMolSubstanceGroups(mol)
    self.assertTrue(len(sgs) == 38)

  def testThreeLetterCodes(self):
    """Test the SCSR system with three letter codes"""

    ofile = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'macromols', 'PepTla.mol')
    with open(ofile) as inf:
      scsrBlock = inf.read()

    scsr = Chem.SCSRFromSCSRBlock(scsrBlock, False, False)
    self.assertTrue(scsr.GetNumTemplates() == 4)
    t1 = scsr.GetTemplate(0)
    self.assertEqual(t1.GetNumAtoms(), 7)
    self.assertEqual(scsr.GetMol().GetNumAtoms(), 4)
    self.assertEqual(scsr.GetMol().GetNumBonds(), 3)

    molFromSCSRParams = Chem.MolFromSCSRParams()
    molFromSCSRParams.includeLeavingGroups = True
    molFromSCSRParams.scsrTemplateNames = Chem.SCSRTemplateNamesVal.SCSRTemplateNamesAsEntered

    mol = Chem.SCSRToMol(scsr, molFromSCSRParams)

    self.assertEqual(mol.GetNumAtoms(), 26)
    self.assertEqual(mol.GetNumBonds(), 25)
    sgs = Chem.GetMolSubstanceGroups(mol)
    self.assertTrue(len(sgs) == 7)
    sgs[0].GetProp('LABEL')
    self.assertEqual(sgs[0].GetProp('LABEL'), 'Ala_4')

    molFromSCSRParams.scsrTemplateNames = Chem.SCSRTemplateNamesVal.SCSRTemplateNamesUseFirstName

    mol = Chem.SCSRToMol(scsr, molFromSCSRParams)

    self.assertEqual(mol.GetNumAtoms(), 26)
    sgs = Chem.GetMolSubstanceGroups(mol)
    self.assertEqual(len(sgs), 7)

    self.assertEqual(sgs[0].GetProp('LABEL'), 'Ala_4')

    molFromSCSRParams.scsrTemplateNames = Chem.SCSRTemplateNamesVal.SCSRTemplateNamesUseLastName

    mol = Chem.SCSRToMol(scsr, molFromSCSRParams)

    self.assertEqual(mol.GetNumAtoms(), 26)
    sgs = Chem.GetMolSubstanceGroups(mol)
    self.assertEqual(len(sgs), 7)

    self.assertEqual(sgs[0].GetProp('LABEL'), 'A_4')


if __name__ == '__main__':
  unittest.main()
