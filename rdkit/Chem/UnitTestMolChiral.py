"""unit testing code for chirality and enantiomer stuff

"""
import unittest, os
import sys
import random
from rdkit import Chem
from rdkit.Chem import MolChiral

class TestCase(unittest.TestCase):

  def testCorrectEnantiomer(self):
    m1 = Chem.MolFromSmiles('CC(=CCC[C@](C)(C=C)O)C') # (-)-linalool
    m2 = Chem.MolFromSmiles('CC(=CCC[C@@](C)(C=C)O)C') # (+)-linalool

    # test that the correct enantiomer is generated
    em1 = MolChiral.GetEnantiomer(m1)

    self.assertEqual(Chem.MolToSmiles(m2), Chem.MolToSmiles(em1))

  def testNoEnantiomer(self):
    m1 = Chem.MolFromSmiles('CCC1(C(=O)NC(=O)NC1=O)CC') # barbital

    # test that None is returned for no enantiomer
    em1 = MolChiral.GetEnantiomer(m1)

    self.assertIsNone(em1)

  def testIsChiral(self):
    m1 = Chem.MolFromSmiles('CCC(=O)OC(CC1=CC=CC=C1)(C2=CC=CC=C2)C(C)CN(C)C') # (+)-propoxyphene

    # test correct recognition of chiral compound
    self.assertTrue(MolChiral.IsMolChiral(m1))
  
  def testIsNotChiral(self):
    m1 = Chem.MolFromSmiles('C1=CC(=O)C=CC1=O') # 1,4-benzoquinone

    # test correct recognition of achiral compound
    self.assertFalse(MolChiral.IsMolChiral(m1))

if __name__ == '__main__':
  unittest.main()
