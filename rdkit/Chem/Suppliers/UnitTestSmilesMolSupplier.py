# $Id$
#
#  Copyright (C) 2003-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""unit testing code for the Smiles file handling stuff

"""
import unittest,sys,os
from rdkit import RDConfig
from rdkit import Chem
from rdkit.six import next

class TestCase(unittest.TestCase):
  def setUp(self):
    self.smis = ['CC','CCC','CCCCC','CCCCCC','CCCCCCC','CC','CCCCOC']

  def test1LazyReader(self):
    " tests lazy reads """
    supp = Chem.SmilesMolSupplierFromText('\n'.join(self.smis),',',0,-1,0)
    for i in range(4):
      m = next(supp)
      assert m,'read %d failed'%i
      assert m.GetNumAtoms(),'no atoms in mol %d'%i
    i = len(supp)-1
    m = supp[i]
    assert m,'read %d failed'%i
    assert m.GetNumAtoms(),'no atoms in mol %d'%i

    ms = [x for x in supp]
    for i in range(len(supp)):
      m = ms[i]
      if m:
        ms[i] = Chem.MolToSmiles(m)

    
    l = len(supp)
    assert l == len(self.smis),'bad supplier length: %d'%(l)

    i = len(self.smis)-3
    m = supp[i-1]
    assert m,'back index %d failed'%i
    assert m.GetNumAtoms(),'no atoms in mol %d'%i
    
    try:
      m = supp[len(self.smis)]
    except:
      fail = 1
    else:
      fail = 0
    assert fail,'out of bound read did not fail'

  def test2LazyIter(self):
    " tests lazy reads using the iterator interface "
    supp = Chem.SmilesMolSupplierFromText('\n'.join(self.smis),',',0,-1,0)

    nDone = 0
    for mol in supp:
      assert mol,'read %d failed'%i
      assert mol.GetNumAtoms(),'no atoms in mol %d'%i
      nDone += 1
    assert nDone==len(self.smis),'bad number of molecules'   

    l = len(supp)
    assert l == len(self.smis),'bad supplier length: %d'%(l)

    i = len(self.smis)-3
    m = supp[i-1]
    assert m,'back index %d failed'%i
    assert m.GetNumAtoms(),'no atoms in mol %d'%i


    try:
      m = supp[len(self.smis)]
    except:
      fail = 1
    else:
      fail = 0
    assert fail,'out of bound read did not fail'


  def test3BoundaryConditions(self):
    smis = ['CC','CCOC','fail','CCO']
    supp = Chem.SmilesMolSupplierFromText('\n'.join(smis),',',0,-1,0)
    assert len(supp)==4
    assert supp[2] is None
    assert supp[3]

    supp = Chem.SmilesMolSupplierFromText('\n'.join(smis),',',0,-1,0)
    assert supp[2] is None
    assert supp[3]
    assert len(supp)==4
    try:
      supp[4]
    except:
      ok=1
    else:
      ok=0
    assert ok

    supp = Chem.SmilesMolSupplierFromText('\n'.join(smis),',',0,-1,0)
    assert len(supp)==4
    assert supp[3]
    try:
      supp[4]
    except:
      ok=1
    else:
      ok=0
    assert ok

    supp = Chem.SmilesMolSupplierFromText('\n'.join(smis),',',0,-1,0)
    try:
      supp[4]
    except:
      ok=1
    else:
      ok=0
    assert ok
    assert len(supp)==4
    assert supp[3]



if __name__ == '__main__':
  unittest.main()

