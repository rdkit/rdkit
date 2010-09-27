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
"""unit testing code for the lazy signature generator

"""
import unittest
from rdkit import Chem
from rdkit.Chem.Pharm2D import SigFactory,LazyGenerator

class TestCase(unittest.TestCase):
  def setUp(self):
    self.factory = SigFactory.SigFactory()
    self.factory.SetPatternsFromSmarts(['O','N'])
    self.factory.SetBins([(0,2),(2,5),(5,8)])
    self.factory.SetMinCount(2)
    self.factory.SetMaxCount(3)

  def test1(self):
    """ simple tests

    """
    mol = Chem.MolFromSmiles('OCC(=O)CCCN')
    sig = self.factory.GetSignature()
    assert sig.GetSize()==105,'bad signature size: %d'%(sig.GetSize())
    sig.SetIncludeBondOrder(0)
    gen = LazyGenerator.Generator(sig,mol)
    assert len(gen) == sig.GetSize(),'length mismatch %d!=%d'%(len(gen),sig.GetSize())
    
    tgt = (1,5,48)
    for bit in tgt:
      assert gen[bit],'bit %d not properly set'%(bit)
      assert gen.GetBit(bit),'bit %d not properly set'%(bit)
      assert not gen[bit+50],'bit %d improperly set'%(bit+100)
    
    sig = self.factory.GetSignature()
    assert sig.GetSize()==105,'bad signature size: %d'%(sig.GetSize())
    sig.SetIncludeBondOrder(1)
    gen = LazyGenerator.Generator(sig,mol)
    assert len(gen) == sig.GetSize(),'length mismatch %d!=%d'%(len(gen),sig.GetSize())

    tgt = (1,4,5,45)
    for bit in tgt:
      assert gen[bit],'bit %d not properly set'%(bit)
      assert gen.GetBit(bit),'bit %d not properly set'%(bit)
      assert not gen[bit+50],'bit %d improperly set'%(bit+100)

    try:
      gen[sig.GetSize()+1]
    except IndexError:
      ok = 1
    else:
      ok = 0
    assert ok,'accessing bogus bit did not fail'
    try:
      gen[-1]
    except IndexError:
      ok = 1
    else:
      ok = 0
    assert ok,'accessing bogus bit did not fail'

if __name__ == '__main__':
  unittest.main()

