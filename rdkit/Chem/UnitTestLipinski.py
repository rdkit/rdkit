# $Id$
#
#  Copyright (C) 2003-2006  Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" unit testing code for Lipinski parameter calculation

  This provides a workout for the SMARTS matcher

"""
from __future__ import print_function
from rdkit import RDConfig
import unittest,os
from rdkit.six.moves import cPickle
from rdkit import Chem
from rdkit.Chem import Lipinski

class TestCase(unittest.TestCase):
  def setUp(self):
    print('\n%s: '%self.shortDescription(),end='')
    self.inFileName = '%s/NCI/first_200.props.sdf'%(RDConfig.RDDataDir)
  def test1(self):
    " testing first 200 mols from NCI "
    suppl = Chem.SDMolSupplier(self.inFileName)
    idx = 1
    oldDonorSmarts = Chem.MolFromSmarts('[NH1,NH2,OH1]')
    OldDonorCount = lambda x,y=oldDonorSmarts:Lipinski._NumMatches(x,y)
    oldAcceptorSmarts = Chem.MolFromSmarts('[N,O]')
    OldAcceptorCount = lambda x,y=oldAcceptorSmarts:Lipinski._NumMatches(x,y)
    for m in suppl:
      if m:
        calc = Lipinski.NHOHCount(m)
        orig = int(m.GetProp('NUM_LIPINSKIHDONORS'))
        assert calc==orig,'bad num h donors for mol %d (%s): %d != %d'%(idx,m.GetProp('SMILES'),calc,orig)

        calc = Lipinski.NOCount(m)
        orig = int(m.GetProp('NUM_LIPINSKIHACCEPTORS'))
        assert calc==orig,'bad num h acceptors for mol %d (%s): %d != %d'%(idx,m.GetProp('SMILES'),calc,orig)

        calc = Lipinski.NumHDonors(m)
        orig = int(m.GetProp('NUM_HDONORS'))
        assert calc==orig,'bad num h donors for mol %d (%s): %d != %d'%(idx,m.GetProp('SMILES'),calc,orig)

        calc = Lipinski.NumHAcceptors(m)
        orig = int(m.GetProp('NUM_HACCEPTORS'))
        assert calc==orig,'bad num h acceptors for mol %d (%s): %d != %d'%(idx,m.GetProp('SMILES'),calc,orig)

        calc = Lipinski.NumHeteroatoms(m)
        orig = int(m.GetProp('NUM_HETEROATOMS'))
        assert calc==orig,'bad num heteroatoms for mol %d (%s): %d != %d'%(idx,m.GetProp('SMILES'),calc,orig)

        calc = Lipinski.NumRotatableBonds(m)
        orig = int(m.GetProp('NUM_ROTATABLEBONDS'))
        assert calc==orig,'bad num rotors for mol %d (%s): %d != %d'%(idx,m.GetProp('SMILES'),calc,orig)
      idx += 1

  def testIssue2183420(self):
    " testing a problem with the acceptor definition "
    self.assertTrue(Lipinski.NumHAcceptors(Chem.MolFromSmiles('NC'))==1)
    self.assertTrue(Lipinski.NumHAcceptors(Chem.MolFromSmiles('CNC'))==1)
    self.assertTrue(Lipinski.NumHAcceptors(Chem.MolFromSmiles('CN(C)C'))==1)
    self.assertTrue(Lipinski.NumHAcceptors(Chem.MolFromSmiles('NC(=O)'))==1)
    self.assertTrue(Lipinski.NumHAcceptors(Chem.MolFromSmiles('NC(=O)C'))==1)
    self.assertTrue(Lipinski.NumHAcceptors(Chem.MolFromSmiles('CNC(=O)'))==1)
    self.assertTrue(Lipinski.NumHAcceptors(Chem.MolFromSmiles('CNC(=O)C'))==1)
    self.assertTrue(Lipinski.NumHAcceptors(Chem.MolFromSmiles('O=CNC(=O)C'))==2)
    self.assertTrue(Lipinski.NumHAcceptors(Chem.MolFromSmiles('O=C(C)NC(=O)C'))==2)
      

if __name__ == '__main__':
  unittest.main()

