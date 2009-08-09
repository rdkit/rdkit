# $Id$
#
#  Copyright (C) 2003-2008  Greg Landrum and Rational Discovery LLC
#         All Rights Reserved
#
""" This is a rough coverage test of the python wrapper

it's intended to be shallow, but broad

"""
from rdkit import RDConfig
import unittest
from rdkit import Chem
from rdkit.Chem import rdchemtransforms as ChemTransforms

class TestCase(unittest.TestCase):
  def setUp(self):
    pass

  def test22DeleteSubstruct(self) :
    query = Chem.MolFromSmarts('C(=O)O')
    mol = Chem.MolFromSmiles('CCC(=O)O')
    nmol = ChemTransforms.DeleteSubstructs(mol, query)
    
    self.failUnless(Chem.MolToSmiles(nmol) == 'CC')

    mol = Chem.MolFromSmiles('CCC(=O)O.O=CO')
    # now delete only fragments
    nmol = ChemTransforms.DeleteSubstructs(mol, query, 1)
    self.failUnless(Chem.MolToSmiles(nmol) == 'CCC(=O)O',Chem.MolToSmiles(nmol))
    
    mol = Chem.MolFromSmiles('CCC(=O)O.O=CO')
    nmol = ChemTransforms.DeleteSubstructs(mol, query, 0)
    self.failUnless(Chem.MolToSmiles(nmol) == 'CC')
    
    mol = Chem.MolFromSmiles('CCCO')
    nmol = ChemTransforms.DeleteSubstructs(mol, query, 0)
    self.failUnless(Chem.MolToSmiles(nmol) == 'CCCO')

    # Issue 96 prevented this from working:
    mol = Chem.MolFromSmiles('CCC(=O)O.O=CO')
    nmol = ChemTransforms.DeleteSubstructs(mol, query, 1)
    self.failUnless(Chem.MolToSmiles(nmol) == 'CCC(=O)O')
    nmol = ChemTransforms.DeleteSubstructs(nmol, query, 1)
    self.failUnless(Chem.MolToSmiles(nmol) == 'CCC(=O)O')
    nmol = ChemTransforms.DeleteSubstructs(nmol, query, 0)
    self.failUnless(Chem.MolToSmiles(nmol) == 'CC')
    
  def test46ReplaceCore(self):
    """ test the ReplaceCore functionality

    """

    core = Chem.MolFromSmiles('C=O')

    smi = 'CCC=O'
    m = Chem.MolFromSmiles(smi)
    r = ChemTransforms.ReplaceCore(m,core)
    self.failUnless(r)
    self.failUnless(Chem.MolToSmiles(r,True)=='[1*]CC')

    smi = 'C1CC(=O)CC1'
    m = Chem.MolFromSmiles(smi)
    r = ChemTransforms.ReplaceCore(m,core)
    self.failUnless(r)
    self.failUnless(Chem.MolToSmiles(r,True) =='[1*]CCCC[2*]')

    smi = 'C1CC(=N)CC1'
    m = Chem.MolFromSmiles(smi)
    r = ChemTransforms.ReplaceCore(m,core)
    self.failIf(r)

if __name__ == '__main__':
  unittest.main()

