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
"""basic unit testing code for SMILES canonicalization

"""
import unittest,os
from rdkit.six.moves import cPickle
from rdkit import Chem

class TestCase(unittest.TestCase):
  def setUp(self):
    pass

  def _testSpellings(self,smiList):
    for smi,spellings in smiList:
      m = Chem.MolFromSmiles(smi)
      canSmi = Chem.MolToSmiles(m)
      for spelling in spellings:
        m = Chem.MolFromSmiles(spelling)
        trySmi = Chem.MolToSmiles(m)
        assert canSmi==trySmi,'Non-canonical: mol %s gave %s (should be %s)'%(spelling,trySmi,canSmi)
        m2 = Chem.MolFromSmiles(trySmi.strip())
        try2 = Chem.MolToSmiles(m2)
        assert canSmi==try2,'Non-canonical: mol %s gave %s (should be %s) on second pass'%(spelling,try2,canSmi)
        
    

  def testLinear1(self):
    " testing first batch of linear mols "
    smiList = [('O=CCO',('OCC=O',
                         'C(O)C=O',
                         'C(C=O)O',
                         'C(CO)=O',
                         )),
               ('OCC(C=C)CCC(C#N)CC',('C=CC(CO)CCC(C#N)CC',
                                      'C(CO)(C=C)CCC(CC)C#N',
                                      'C(CO)(C=C)CCC(C#N)CC',
                                      'C(C=C)(CO)CCC(C#N)CC',
                                      'C(C=C)(CO)CCC(CC)C#N',
                                      )),
               ('[Se]=CCO',('OCC=[Se]',
                         'C(O)C=[Se]',
                         'C(C=[Se])O',
                         'C(CO)=[Se]',
                         )),
               ('O=C(OC)C',('COC(=O)C','O=C(C)OC','O(C)C(=O)C')),
               ]
    self._testSpellings(smiList)
  def testRings1(self):
    " testing first batch of rings "
    smiList = [('C1OCCCC1',('O1CCCCC1',
                            'C1COCCC1',
                            'C1CCOCC1',
                            'C1CCCOC1',
                            'C1CCCCO1',)),
               ('CC1=CCCCC1',('C1=C(C)CCCC1',
                              'C1CC=C(C)CC1',)),
               ('CC1C=CCCC1',('C1=CC(C)CCC1',
                              'C1CC=CC(C)C1',)),
               ]
    self._testSpellings(smiList)

  def testRings2(self):
    " testing second batch of rings "
    smiList = [('c1c(cc2nc3cc(ccc3cc2c1))',('c1ccc2cc3ccccc3nc2c1',
                                            'c1ccc2nc3ccccc3cc2c1',
                                            'c1c2nc3ccccc3cc2ccc1')),
               ('Cc1ccc2nc3ccccc3cc2c1',('c1ccc2nc3ccc(C)cc3cc2c1',
                                          )),
               ('c1c(C)cc2nc3ccccc3cc2c1',('c1ccc2nc3cc(C)ccc3cc2c1',
                                           )),

               ]
    self._testSpellings(smiList)
      
  def testProblems(self):
    " testing molecules which have been problematic "
    smiList = [ ('[Al+3]CCC',('CCC[Al+3]','C(C)(C[Al+3])' ) ),
                ('C(=O)(Cl)CC(=O)Cl',('ClC(CC(Cl)=O)=O','C(Cl)(=O)CC(=O)Cl','C(Cl)(=O)CC(Cl)=O')),
                ('C(=O)(Cl)c1ccc(C(=O)Cl)cc1',('O=C(Cl)c1ccc(cc1)C(Cl)=O','C(Cl)(=O)C1=CC=C(C=C1)C(Cl)=O',
                                               'ClC(=O)c1ccc(cc1)C(=O)Cl')),
                ('[N+](=O)([O-])c1ccc([N+](=O)[O-])cc1',('[N+]([O-])(=O)C1=CC=C(C=C1)[N+](=O)[O-]',
                                                         'O=[N+1]([O-1])c1ccc(cc1)[N+1]([O-1])=O',
                                                         '[O-1][N+1](=O)c1ccc(cc1)[N+1]([O-1])=O')),
                ('Oc1c3c(cc(c1)S(=O)(=O)O)cc(NC(=O)c2ccccc2)cc3',
                 ('C1=C(C2=C(C=C1S(O)(=O)=O)C=C(C=C2)NC(C3=CC=CC=C3)=O)O',
                  'O=S(=O)(O)c1cc(O)c2ccc(NC(=O)c3ccccc3)cc2c1',
                  'OS(=O)(=O)c1cc(O)c2ccc(NC(=O)c3ccccc3)cc2c1')),
                #('C(C(C)(C)O)C(C)O',('C([C@@](C)(O)C)C(C)O','CC(O)CC(C)(O)C','CC(O)CC(O)(C)C')),
                ('C',('C')),
                ('C(Cl)(Br)(F)CC(Cl)(Br)(F)',('C(Cl)(F)(Br)CC(F)(Cl)(Br)','C(Cl)(Br)(F)CC(Cl)(F)(Br)',
                                              'C(F)(Br)(Cl)CC(Br)(Cl)(F)','C(C(Cl)(Br)(F))C(F)(Cl)Br')),
                  
               ]
    self._testSpellings(smiList)


  def testHighSymmetry(self):
    " testing tricky (high-symmetry) molecules "
    smiList = [('CC(C)CC',('CCC(C)C',)),
               ('C1CCCC1CCC',('CCCC1CCCC1',)),
               ('C1(C)CC(C)CCC1',('CC1CCCC(C)C1',)),
               ('CCC1CCCCC1CC',('CCC1CCCCC1CC',)),
               ('CCC1CC(CC)CCC1',('CCC1CCCC(CC)C1',)),
               ('C1CCCCC1CC(CC)CC',('CCC(CC)CC1CCCCC1',)),
               ('C1CCCC2C1CC(CC)CC2',('CCC1CCC2CCCCC2C1',)),
               ('CC1CCCC2C1C(C)CCC2',('CC1CCCC2CCCC(C)C12',)),
               ('C2CCC1CCC(C)C12',('CC1CCC2CCCC12',)),
               ('CC(C)CCCC(C)C',('CC(CCCC(C)C)C',)),
               ]
    self._testSpellings(smiList)


  def testFailures(self):
    " EXPECT FAILURES -> testing molecules which are known to fail "
    smiList = [('C13C6C1C2C4C2C3C5C4C56',
                ('C45C1C6C3C6C5C4C2C3C12',
                 'C45C2C6C3C6C5C4C1C3C12')),
                ]
    self._testSpellings(smiList)


if __name__ == '__main__':
  unittest.main()

