# $Id$
#
#  Copyright (C) 2002-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""unit testing code for the EState indices

validation values are from the paper (JCICS _31_ 76-81 (1991))

"""
from __future__ import print_function
import unittest
import numpy
from rdkit import Chem
from rdkit.Chem import EState

class TestCase(unittest.TestCase):
  def setUp(self):
    #print '\n%s: '%self.shortDescription(),
    pass

  def _validate(self,vals,tol=1e-2,show=0):
    for smi,ans in vals:
      mol = Chem.MolFromSmiles(smi)
      ans = numpy.array(ans)
      inds = EState.EStateIndices(mol)

      maxV = max(abs(ans-inds))
      if show: print(inds)
      assert maxV<tol,'bad EStates for smiles: %s'%(smi)
      
  
  def test1(self):
    """ simple molecules
    """
    data = [
      ('CCCC',[2.18,1.32,1.32,2.18]),
      ('CCCCC',[2.21,1.34,1.39,1.34,2.21]),
      ('CCCCCCC',[2.24,1.36,1.42,1.44,1.42,1.36,2.24]),
      ('CCCCCCCCCC',[2.27,1.37,1.44,1.46,1.47,1.47,1.46,1.44,1.37,2.27]),
      ]
    self._validate(data)
    
  def test2(self):
    """ isomers
    """
    data = [
      ('CCCCCC',[2.23,1.36,1.41,1.41,1.36,2.23]),
      ('CCC(C)CC',[2.23,1.33,0.94,2.28,1.33,2.23]),
      ('CC(C)CCC',[2.25,0.90,2.25,1.38,1.33,2.22]),
      ('CC(C)(C)CC',[2.24,0.54,2.24,2.24,1.27,2.20]),
      ]
    self._validate(data)
    
  def test3(self):
    """ heteroatoms
    """
    data = [
      ('CCCCOCCCC',[2.18,1.24,1.21,0.95,5.31,0.95,1.21,1.24,2.18]),
      ('CCC(C)OC(C)CC',[2.15,1.12,0.43,2.12,5.54,0.43,2.12,1.12,2.15]),
      ('CC(C)(C)OC(C)(C)C',[2.07,-0.02,2.07,2.07,5.63,-0.02,2.07,2.07,2.07]),
      ('CC(C)CC',[2.22,0.88,2.22,1.31,2.20]),
      ('CC(C)CN',[2.10,0.66,2.10,0.81,5.17]),
      ('CC(C)CO',[1.97,0.44,1.97,0.31,8.14]),
      ('CC(C)CF',[1.85,0.22,1.85,-0.19,11.11]),
      ('CC(C)CCl',[2.09,0.65,2.09,0.78,5.34]),
      ('CC(C)CBr',[2.17,0.80,2.17,1.11,3.31]),
      ('CC(C)CI',[2.21,0.87,2.21,1.28,2.38]),
      ]
    self._validate(data,show=0)
    
  def test4(self):
    """ more heteroatoms
    """
    data = [
      ('CC(N)C(=O)O',[1.42,-0.73,4.84,-0.96,9.57,7.86]),
      ('CCOCC',[1.99,0.84,4.83,0.84,1.99]),
      ('CCSCC',[2.17,1.26,1.96,1.26,2.17]),  #NOTE: this does not match the values in the paper
      ('CC(=O)OC',[1.36,-0.24,9.59,4.11,1.35]),
      ('CC(=S)OC',[1.73,0.59,4.47,4.48,1.56]),
      ]
    self._validate(data,show=0)
    
  def test5(self):
    """ aromatics with heteroatoms
    """
    data = [
      ('Fc1ccc(C)cc1',[12.09,-0.17,1.45,1.75,1.09,1.93,1.75,1.45]),
      ('Clc1ccc(C)cc1',[5.61,0.80,1.89,1.99,1.24,2.04,1.99,1.89]),
      ('Brc1ccc(C)cc1',[3.35,1.14,2.04,2.07,1.30,2.08,2.07,2.04]),
      ('Ic1ccc(C)cc1',[2.30,1.30,2.10,2.11,1.32,2.09,2.11,2.10]),
      ]
    self._validate(data,show=0)
    
if __name__ == '__main__':
  unittest.main()

