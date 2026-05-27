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

    molFromMACROMolParams = Chem.MolFromMACROMolParams()
    molFromMACROMolParams.includeLeavingGroups = True
    molFromMACROMolParams.macroTemplateNames = Chem.MACROTemplateNames.AsEntered
    scsrBaseHbondOptions = Chem.SCSRBaseHbondOptions()
    scsrBaseHbondOptions = Chem.SCSRBaseHbondOptions.Ignore

    for mol in (Chem.MolFromSCSRBlock(scsrBlock, False, False, molFromMACROMolParams, scsrBaseHbondOptions),
                Chem.MolFromSCSRFile(ofile, False, False, molFromMACROMolParams, scsrBaseHbondOptions)):

      self.assertTrue(mol.GetNumAtoms() == 30)
      sgs = Chem.GetMolSubstanceGroups(mol)
      self.assertEqual(len(sgs), 6)

    # check defaults:
    molFromMACROMolParamsDefault = Chem.MolFromMACROMolParams()
    for mol in (Chem.MolFromSCSRBlock(scsrBlock,True,True, molFromMACROMolParamsDefault, scsrBaseHbondOptions), Chem.MolFromSCSRFile(ofile,True,True, molFromMACROMolParamsDefault,scsrBaseHbondOptions)):
      self.assertTrue(mol.GetNumAtoms() == 30)
      sgs = Chem.GetMolSubstanceGroups(mol)
      self.assertEqual(len(sgs), 6)

  def testSCSRRna(self):
    """Test the SCSR system with and RNA double strand"""

    ofile = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'macromols', 'DnaTest.mol')
    with open(ofile) as inf:
      scsrBlock = inf.read()

    molFromMACROMolParams = Chem.MolFromMACROMolParams()
    molFromMACROMolParams.includeLeavingGroups = True
    molFromMACROMolParams.macroTemplateNames = Chem.MACROTemplateNames.AsEntered
    scsrBaseHbondOptions = Chem.SCSRBaseHbondOptions()
    scsrBaseHbondOptions = Chem.SCSRBaseHbondOptions.Auto

    mol = Chem.MolFromSCSRBlock(scsrBlock, False, False, molFromMACROMolParams, scsrBaseHbondOptions)

    self.assertTrue(mol.GetNumAtoms() == 254)
    self.assertTrue(mol.GetNumBonds() == 300)
    sgs = Chem.GetMolSubstanceGroups(mol)
    self.assertTrue(len(sgs) == 38)

    scsrBaseHbondOptions = Chem.SCSRBaseHbondOptions()
    scsrBaseHbondOptions = Chem.SCSRBaseHbondOptions.Ignore
    mol = Chem.MolFromSCSRBlock(scsrBlock, False, False, molFromMACROMolParams, scsrBaseHbondOptions)

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

    molFromMACROMolParams = Chem.MolFromMACROMolParams()
    molFromMACROMolParams.includeLeavingGroups = True
    molFromMACROMolParams.macroTemplateNames = Chem.MACROTemplateNames.AsEntered
    scsrBaseHbondOptions = Chem.SCSRBaseHbondOptions()

    mol = Chem.MolFromSCSRBlock(scsrBlock, False, False, molFromMACROMolParams,scsrBaseHbondOptions)

    self.assertEqual(mol.GetNumAtoms(), 26)
    self.assertEqual(mol.GetNumBonds(), 25)
    sgs = Chem.GetMolSubstanceGroups(mol)
    self.assertTrue(len(sgs) == 7)
    sgs[0].GetProp('LABEL')
    self.assertEqual(sgs[0].GetProp('LABEL'), 'AA_4_Ala')

    molFromMACROMolParams.macroTemplateNames = Chem.MACROTemplateNames.UseFirstName
    scsrBaseHbondOptions = Chem.SCSRBaseHbondOptions()

    mol = Chem.MolFromSCSRBlock(scsrBlock, False, False, molFromMACROMolParams, scsrBaseHbondOptions)

    self.assertEqual(mol.GetNumAtoms(), 26)
    sgs = Chem.GetMolSubstanceGroups(mol)
    self.assertEqual(len(sgs), 7)

    self.assertEqual(sgs[0].GetProp('LABEL'), 'AA_4_Ala')

    molFromMACROMolParams.macroTemplateNames = Chem.MACROTemplateNames.UseSecondName
    mol = Chem.MolFromSCSRBlock(scsrBlock, False, False, molFromMACROMolParams, scsrBaseHbondOptions)

    self.assertEqual(mol.GetNumAtoms(), 26)
    sgs = Chem.GetMolSubstanceGroups(mol)
    self.assertEqual(len(sgs), 7)

    self.assertEqual(sgs[0].GetProp('LABEL'), 'AA_4_A')


if __name__ == '__main__':
  unittest.main()
