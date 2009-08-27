# $Id$
#
#  Copyright (C) 2001-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
"""testing code inspired by old bugs

The bugs were in the OELib code, so these are maybe no longer
relevant... but tests are tests

"""
from rdkit import RDConfig
import unittest,cPickle,os
from rdkit import Chem

def feq(n1,n2,tol=1e-4):
  return abs(n1-n2)<=tol

class TestCase(unittest.TestCase):
  def testBug14(self):
    """ 
    
    """
    smi = '[O-][N+](=O)C1=CNC(=N)S1'
    mol = Chem.MolFromSmiles(smi)
    at = mol.GetAtomWithIdx(5)
    assert at.GetHybridization()==Chem.HybridizationType.SP2,'bad hyb'
    assert at.GetTotalNumHs()==1,'bad H count'

    mol = Chem.MolFromSmiles(smi)
    at = mol.GetAtomWithIdx(5)
    assert at.GetTotalNumHs()==1,'bad H count'
    assert at.GetHybridization()==Chem.HybridizationType.SP2,'bad hyb'


if __name__ == '__main__':
  unittest.main()

