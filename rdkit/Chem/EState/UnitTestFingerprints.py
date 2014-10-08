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
"""unit testing code for the EState fingerprinting

validation values are from the paper (JCICS _35_ 1039-1045 (1995))

"""
from __future__ import print_function
import unittest
import numpy
from rdkit import Chem
from rdkit.Chem import EState
from rdkit.Chem.EState import Fingerprinter

class TestCase(unittest.TestCase):
  def setUp(self):
    pass

  def _validate(self,vals,tol=1e-2,show=0):
    for smi,c,v in vals:
      mol = Chem.MolFromSmiles(smi)
      counts,vals = Fingerprinter.FingerprintMol(mol)
      counts = counts[numpy.nonzero(counts)]
      vals = vals[numpy.nonzero(vals)]
      if show:
        print(counts)
        print(vals)
      assert len(c)==len(counts),'bad count len for smiles: %s'%(smi)
      assert len(v)==len(vals),'bad val len for smiles: %s'%(smi)
      c = numpy.array(c)
      assert max(abs(c-counts))<tol,'bad count for SMILES: %s'%(smi)
      v = numpy.array(v)
      assert max(abs(v-vals))<tol,'bad val for SMILES: %s'%(smi)
      
  def test1(self):
    """ molecules
    """
    data = [
      ('c1[nH]cnc1CC(N)C(O)=O',[1,2,1,1,1,1,1,1,1,1],
       [0.26,3.12,-0.86,-1.01,0.67,5.25,2.71,3.84,8.42,10.26]),
      ('NCCc1ccc(O)c(O)c1',[2,3,3,1,2],
       [1.26,4.71,0.75,5.30,17.97]),
      ]
    self._validate(data,show=0)
    
    
if __name__ == '__main__':
  unittest.main()

