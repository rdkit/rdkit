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
from Chem.Pharm2D import Signature,Generate,SigFactory,Utils
import cPickle

class TestCase(unittest.TestCase):
  def setUp(self):
    self.factory = SigFactory.SigFactory()
    self.factory.SetPatternsFromSmarts(['O','N'])
    self.factory.SetBins([(0,2),(2,5),(5,8)])
    self.factory.SetMinCount(2)
    self.factory.SetMaxCount(3)
    
  def testSizes(self):
    sig = Signature.Pharm2DSig()
    sig.SetPatternsFromSmarts(['O','N'])
    sig.SetBins([(1,2),(2,5),(5,6)])
    sig.SetMinCount(2)

    sig.SetMaxCount(2)
    sig.Init()
    assert sig.GetSize()==9,'bad 2 point size %d'%(sig.GetSize())
    
    sig.SetMaxCount(3)
    sig.Init()
    assert sig.GetSize()==105,'bad 3 point size %d'%(sig.GetSize())
    
    sig.SetMaxCount(4)
    sig.Init()
    assert sig.GetSize()==1075,'bad 4 point size %d'%(sig.GetSize())
    
  def testBitIdx(self,sig=None):
    data = [
      ( (0,0),[2],1 ),
      ( (0,0),[5],2 ),
      ( (0,1),[5],5 ),
      ( (1,1),[4],7 ),
      ( (1,1),[7],8 ),
      ( (0,0,0),[1,1,1],9),
      ( (0,0,1),[1,1,1],33),
      ( (0,0,1),[1,1,3],34),
      ( (0,0,1),[3,1,1],40),
      ( (0,0,1),[3,3,1],43),
      ]
    if sig is None:
      sig = self.factory.GetSignature()
    for tpl in data:
      patts,dists,bit = tpl
      try:
        idx = sig.GetBitIdx(patts,dists)
      except:
        assert 0,'GetBitIdx failed for probe %s'%(str(tpl))
      else:
        assert bit==idx,'bad idx (%d) for probe %s'%(idx,str(tpl))
            
  def testSimpleSig(self,sig=None):
    if sig is None:
      sig = Signature.Pharm2DSig()
      sig.SetPatternsFromSmarts(['O'])
      sig.SetBins([(1,3),(3,4),(4,8)])
      sig.SetMinCount(2)
      sig.SetMaxCount(3)
      sig.Init()

    mol = Chem.MolFromSmiles('OCCC1COCCO1')
    Generate.Gen2DFingerprint(mol,sig)
    assert sig.GetSize()==30,'bad sig size: %d'%(sig.GetSize())
    bs = tuple(sig.GetOnBits())
    assert bs==(1,2,20),'bad bit list: %s'%(str(bs))

  def testSimpleSig2(self,sig=None):
    if sig is None:
      sig = Signature.Pharm2DSig()
      sig.SetPatternsFromSmarts(['O'])
      sig.SetBins([(1,3),(3,4),(4,8)])
      sig.SetMinCount(2)
      sig.SetMaxCount(3)
      sig.Init()

    mol = Chem.MolFromSmiles('OCCC1COCCO1')
    Generate.Gen2DFingerprint(mol,sig)
    assert sig.GetSize()==30,'bad sig size: %d'%(sig.GetSize())
    bs = tuple(sig.GetOnBits())
    assert bs==(1,2,20),'bad bit list: %s'%(str(bs))

  def testBitIDs1(self):
    """ test 3 point p'cophore ids,
    you can never have too much of this stuff

    """
    sig = self.factory.GetSignature()
    sig.SetBins(((0,2),(2,4),(4,8)))
    sig.Init()
    assert sig.GetSize()==117,'bad signature size: %d'%(sig.GetSize())
    probes = [((0,0,0),(1,3,1),12),
              ((0,0,0),(1,3,3),13),
              ((0,0,1),(1,3,1),39),
              ]
    for patts,bins,ans in probes:
      idx = sig.GetBitIdx(patts,bins)
      assert idx==ans,'bad idx: %d != %d'%(idx,ans)
    patts,bins = (1,0,0),(1,3,1)
    try:
      sig.GetBitIdx(patts,bins,checkPatts=1)
    except ValueError:
      pass
    except:
      assert 0,'bad exception type'
    else:
      assert 0,'expected exception was not raised'
  
    # we don't bother checking the return value here because it's bogus
    try:
      sig.GetBitIdx(patts,bins,checkPatts=0)
    except:
      assert 0,'should not have thrown an exception here'


    

  def testBitIDs2(self):
    """ test 3 point p'cophore ids where the triangle
      inequality has been used to remove some bits
    """
    sig = self.factory.GetSignature()
    assert sig.GetSize()==105,'bad signature size: %d'%(sig.GetSize())
    probes = [((0,0,0),(1,3,1),11),
              ((0,0,0),(1,3,3),12),
              ((0,0,1),(1,3,1),35),
              ]
    for patts,bins,ans in probes:
      idx = sig.GetBitIdx(patts,bins)
      assert idx==ans,'bad idx: %d != %d'%(idx,ans)
    patts,bins = (1,0,0),(1,3,1)
    try:
      sig.GetBitIdx(patts,bins,checkPatts=1)
    except ValueError:
      pass
    except:
      assert 0,'bad exception type'
    else:
      assert 0,'expected exception was not raised'
    # we don't bother checking the return value here because it's bogus
    try:
      sig.GetBitIdx(patts,bins,checkPatts=0)
    except:
      assert 0,'should not have thrown an exception here'

  

  # FIX: add test for perms argument to Gen2DFingerprint

  def testBondOrderSigs(self):
    """ test sigs where bond order is used

    """
    mol = Chem.MolFromSmiles('OCCC(=O)CCCN')
    sig = self.factory.GetSignature()
    assert sig.GetSize()==105,'bad signature size: %d'%(sig.GetSize())
    sig.SetIncludeBondOrder(0)

    Generate.Gen2DFingerprint(mol,sig)
    tgt = (1,5,48)
    onBits = tuple(sig.GetOnBits())
    assert len(onBits)==len(tgt),'bad on bit length (%d!=%d)'%(len(onBits),len(tgt))
    assert onBits==tgt,'bad on bits (%s != %s)'%(onBits,tgt)
    
    mol = Chem.MolFromSmiles('OCCC(=O)CCCN')
    sig = self.factory.GetSignature()
    sig.SetIncludeBondOrder(1)
    Generate.Gen2DFingerprint(mol,sig)
    tgt = (1,4,5,45)
    onBits = tuple(sig.GetOnBits())
    assert len(onBits)==len(tgt),'bad on bit length (%d!=%d)'%(len(onBits),len(tgt))

    assert onBits==tgt,'bad on bits (%s != %s)'%(onBits,tgt)
    
  def testPickle(self):
    onBits = (12,25)
    sig = self.factory.GetSignature()
    self.testBitIdx(sig=sig)
    for bit in onBits:
      sig._bv.SetBit(bit)
    
    text = cPickle.dumps(sig)
    sig = cPickle.loads(text)
    assert tuple(sig._bv.GetOnBits())==onBits,'bad onbits (%s != %s)'%(str(sig._bv.GetOnBits()),
                                                                str(onBits))
    self.testBitIdx(sig=sig)


if __name__ == '__main__':
  unittest.main()

