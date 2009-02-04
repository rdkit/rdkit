# $Id$
#
#  Copyright (C) 2002-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
"""unit testing code for the signatures

"""
import unittest
import Chem
from Chem.Pharm2D import Generate,SigFactory,Matcher,Gobbi_Pharm2D

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
    Generate.Gen2DFingerprint(mol,sig)
    tgt = (1,5,48)
    onBits = sig.GetOnBits()
    assert len(onBits)==len(tgt),'bad on-bit length (%d!=%d)'%(len(onBits),len(tgt))
    for i in range(len(onBits)):
      self.failUnless(onBits[i]==tgt[i],'bad on-bits (%s != %s)'%(list(onBits),tgt))

    bitMatches = ([((0,),(3,))],
                  [((0,),(7,)),((3,),(7,))],
                  [((0,),(3,),(7,))],
                  )
    for i in range(len(onBits)):
      bit = onBits[i]
      tgt = bitMatches[i]
      matches = Matcher.GetAtomsMatchingBit(sig,bit,mol)
      assert len(matches)==len(tgt),\
             'bad match length for bit %d (%d != %d)'%(bit,len(matches),len(tgt))
      assert matches==tgt,\
             'bad match for bit %d (%s != %s)'%(bit,matches,tgt)

    
    sig = self.factory.GetSignature()
    assert sig.GetSize()==105,'bad signature size: %d'%(sig.GetSize())
    sig.SetIncludeBondOrder(1)
    Generate.Gen2DFingerprint(mol,sig)
    tgt = (1,4,5,45)
    onBits = sig.GetOnBits()
    assert len(onBits)==len(tgt),'bad on-bit length (%d!=%d)'%(len(onBits),len(tgt))
    for i in range(len(onBits)):
      self.failUnless(onBits[i]==tgt[i],'bad on-bits (%s != %s)'%(list(onBits),tgt))

    bitMatches = ([((0,),(3,))],
                  [((3,),(7,))],
                  [((0,),(7,))],
                  [((0,),(3,),(7,))],
                  )
    for i in range(len(onBits)):
      bit = onBits[i]
      tgt = bitMatches[i]
      matches = Matcher.GetAtomsMatchingBit(sig,bit,mol)
      assert len(matches)==len(tgt),\
             'bad match length for bit %d (%d != %d)'%(bit,len(matches),len(tgt))
      assert matches==tgt,\
             'bad match for bit %d (%s != %s)'%(bit,matches,tgt)


  def testBug28(self):
    smi = 'Cc([s]1)nnc1SCC(\CS2)=C(/C([O-])=O)N3C(=O)[C@H]([C@@H]23)NC(=O)C[n]4cnnn4'
    mol = Chem.MolFromSmiles(smi)
    factory = Gobbi_Pharm2D.factory
    factory.SetBins([(2,3),(3,4),(4,5),(5,8),(8,100)])
    sig = Generate.Gen2DFingerprint(mol,factory)
    onBits = sig.GetOnBits()
    for bit in onBits:
      as = Matcher.GetAtomsMatchingBit(sig,bit,mol,justOne=1)
      assert len(as),'bit %d failed to match'%(bit)

  def testRoundtrip(self):
    """ longer-running Bug 28 test
    """
    import RDConfig,os
    nToDo=20
    inD = open(os.path.join(RDConfig.RDDataDir,'NCI','first_5K.smi'),'r').readlines()[:nToDo]
    factory = Gobbi_Pharm2D.factory
    factory.SetBins([(2,3),(3,4),(4,5),(5,8),(8,100)])
    for line in inD:
      smi = line.split('\t')[0]
      mol = Chem.MolFromSmiles(smi)
      sig = Generate.Gen2DFingerprint(mol,factory)
      onBits = sig.GetOnBits()
      for bit in onBits:
        as = Matcher.GetAtomsMatchingBit(sig,bit,mol,justOne=1)
        assert len(as),'bit %d failed to match for smi %s'%(bit,smi)

    

if __name__ == '__main__':
  unittest.main()

