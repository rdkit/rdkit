# $Id: UnitTestLipinski.py 5007 2006-02-22 15:14:41Z glandrum $
#
#  Copyright (C) 2003-2006  Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" unit testing code for Lipinski parameter calculation

  This provides a workout for the SMARTS matcher

"""
import RDConfig
import unittest,cPickle,os
import Chem
from Chem import Lipinski

class TestCase(unittest.TestCase):
  def setUp(self):
    print '\n%s: '%self.shortDescription(),
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
        calc = OldDonorCount(m)
        orig = int(m.GetProp('NUM_HDONORS'))
        assert calc==orig,'bad num h donors for mol %d (%s): %d != %d'%(idx,m.GetProp('SMILES'),calc,orig)

        calc = OldAcceptorCount(m)
        orig = int(m.GetProp('NUM_HACCEPTORS'))
        assert calc==orig,'bad num h acceptors for mol %d (%s): %d != %d'%(idx,m.GetProp('SMILES'),calc,orig)

        calc = Lipinski.NumHeteroatoms(m)
        orig = int(m.GetProp('NUM_HETEROATOMS'))
        assert calc==orig,'bad num heteroatoms for mol %d (%s): %d != %d'%(idx,m.GetProp('SMILES'),calc,orig)

        calc = Lipinski.NumRotatableBonds(m)
        orig = int(m.GetProp('NUM_ROTATABLEBONDS'))
        assert calc==orig,'bad num rotors for mol %d (%s): %d != %d'%(idx,m.GetProp('SMILES'),calc,orig)
      idx += 1
      

if __name__ == '__main__':
  unittest.main()

