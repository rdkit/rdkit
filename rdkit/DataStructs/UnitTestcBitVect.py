# $Id$
#
# Copyright (C) 2001-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" unit testing code for the C++ BitVects
"""
from __future__ import print_function
import unittest,os,sys
from rdkit.six.moves import cPickle
from rdkit.DataStructs import cDataStructs
klass = cDataStructs.SparseBitVect

def feq(n1,n2,tol=1e-4):
  return abs(n1-n2)<=tol
def ieq(n1,n2):
  return abs(n1-n2)==0

class TestCase(unittest.TestCase):
  def setUp(self):
    print('\n%s: '%self.shortDescription(),end='')
    sys.stdout.flush()
  def testSparseIdx(self):
    """ test indexing into SparseBitVects
    """
    v = klass(10)
    ok = 1
    v[0] = 1
    v[2] = 1
    v[9] = 1
    
    try:
      v[10] = 1
    except IndexError:
      ok = 1
    except:
      assert 0, 'setting high bit should have failed with an IndexError'
    else:
      assert 0, 'setting high bit should have failed'

    assert v[0] == 1, 'bad bit'
    assert v[1] == 0, 'bad bit'
    assert v[2] == 1, 'bad bit'
    assert v[9] == 1, 'bad bit'
    assert v[-1] == 1, 'bad bit'
    assert v[-2] == 0, 'bad bit'

    try:
      foo = v[10]
    except IndexError:
      ok = 1
    except:
      assert 0, 'getting high bit should have failed with an IndexError'
    else:
      assert 0, 'getting high bit should have failed'
    

  def testSparseBitGet(self):
    """ test operations to get sparse bits
    """
    v = klass(10)
    v[0] = 1
    v[2] = 1
    v[6] = 1
    assert len(v)==10,'len(SparseBitVect) failed'
    assert v.GetNumOnBits()==3,'NumOnBits failed'
    assert tuple(v.GetOnBits())==(0,2,6), 'GetOnBits failed'
    
  def testSparseBitOps(self):
    """ test bit operations on SparseBitVects
    """
    v1 = klass(10)
    v1[0] = 1
    v1[2] = 1
    v1[6] = 1
    v2 = klass(10)
    v2[0] = 1
    v2[3] = 1
    v2[6] = 1
    assert tuple((v1&v2).GetOnBits()) == (0,6),'binary & failed'
    assert tuple((v1&v2).GetOnBits()) == (0,6),'binary & failed'
    assert tuple((v1|v2).GetOnBits()) == (0,2,3,6),'binary | failed'
    assert tuple((v1^v2).GetOnBits()) == (2,3),'binary ^ failed'
    
  def testTanimotoSim(self):
    """ test Tanimoto Similarity measure
    """
    v1 = klass(10)
    v1[0] = 1
    v1[2] = 1
    v1[6] = 1
    v2 = klass(10)
    v2[0] = 1
    v2[3] = 1
    v2[6] = 1
    v3 = klass(10)
    v3[1] = 1
    v3[4] = 1
    v3[8] = 1
    assert feq(cDataStructs.TanimotoSimilarity(v1,v1),1.0),'bad v1,v1 TanimotoSimilarity'
    assert feq(cDataStructs.TanimotoSimilarity(v2,v2),1.0),'bad v2,v2 TanimotoSimilarity'
    assert feq(cDataStructs.TanimotoSimilarity(v1,v2),0.5),'bad v1,v2 TanimotoSimilarity'
    assert feq(cDataStructs.TanimotoSimilarity(v2,v1),0.5),'bad v2,v1 TanimotoSimilarity'
    assert feq(cDataStructs.TanimotoSimilarity(v1,v3),0.0),'bad v1,v3 TanimotoSimilarity'
    assert feq(cDataStructs.TanimotoSimilarity(v2,v3),0.0),'bad v2,v3 TanimotoSimilarity'

  def testOnBitSim(self):
    """ test On Bit Similarity measure
    """
    v1 = klass(10)
    v1[0] = 1
    v1[2] = 1
    v1[6] = 1
    v2 = klass(10)
    v2[0] = 1
    v2[3] = 1
    v2[6] = 1
    v3 = klass(10)
    v3[1] = 1
    v3[4] = 1
    v3[8] = 1
    assert feq(cDataStructs.OnBitSimilarity(v1,v1),1.0),'bad v1,v1 OnBitSimilarity'
    assert feq(cDataStructs.OnBitSimilarity(v2,v2),1.0),'bad v2,v2 OnBitSimilarity'
    assert feq(cDataStructs.OnBitSimilarity(v1,v2),0.5),'bad v1,v2 OnBitSimilarity'
    assert feq(cDataStructs.OnBitSimilarity(v2,v1),0.5),'bad v2,v1 OnBitSimilarity'
    assert feq(cDataStructs.OnBitSimilarity(v1,v3),0.0),'bad v1,v3 OnBitSimilarity'
    assert feq(cDataStructs.OnBitSimilarity(v2,v3),0.0),'bad v2,v3 OnBitSimilarity'

  def testNumBitsInCommon(self):
    """ test calculation of Number of Bits in Common
    """
    v1 = klass(10)
    v1[0] = 1
    v1[2] = 1
    v1[6] = 1
    v2 = klass(10)
    v2[0] = 1
    v2[3] = 1
    v2[6] = 1
    v3 = klass(10)
    v3[1] = 1
    v3[4] = 1
    v3[8] = 1
    assert ieq(cDataStructs.NumBitsInCommon(v1,v1),10),'bad v1,v1 NumBitsInCommon'
    assert ieq(cDataStructs.NumBitsInCommon(v2,v2),10),'bad v2,v2 NumBitsInCommon'
    assert ieq(cDataStructs.NumBitsInCommon(v1,v2),8),'bad v1,v2 NumBitsInCommon'
    assert ieq(cDataStructs.NumBitsInCommon(v2,v1),8),'bad v2,v1 NumBitsInCommon'
    assert ieq(cDataStructs.NumBitsInCommon(v1,v3),4),'bad v1,v3 NumBitsInCommon'
    assert ieq(cDataStructs.NumBitsInCommon(v2,v3),4),'bad v2,v3 NumBitsInCommon'

  def testAllBitSim(self):
    """ test All Bit Similarity measure
    """
    v1 = klass(10)
    v1[0] = 1
    v1[2] = 1
    v1[6] = 1
    v2 = klass(10)
    v2[0] = 1
    v2[3] = 1
    v2[6] = 1
    v3 = klass(10)
    v3[1] = 1
    v3[4] = 1
    v3[8] = 1
    assert feq(cDataStructs.AllBitSimilarity(v1,v1),1.0),'bad v1,v1 AllBitSimilarity'
    assert feq(cDataStructs.AllBitSimilarity(v2,v2),1.0),'bad v2,v2 AllBitSimilarity'
    assert feq(cDataStructs.AllBitSimilarity(v1,v2),0.8),'bad v1,v2 AllBitSimilarity'
    assert feq(cDataStructs.AllBitSimilarity(v2,v1),0.8),'bad v2,v1 AllBitSimilarity'
    assert feq(cDataStructs.AllBitSimilarity(v1,v3),0.4),'bad v1,v3 AllBitSimilarity'
    assert feq(cDataStructs.AllBitSimilarity(v2,v3),0.4),'bad v2,v3 AllBitSimilarity'

  def testStringOps(self):
    """ test serialization operations
    """
    v1 = klass(10)
    v1[0] = 1
    v1[2] = 1
    v1[6] = 1
    s = v1.ToBinary()
    v2 = klass(s)
    assert tuple(v2.GetOnBits())==tuple(v1.GetOnBits()),'To/From string failed'

  def testOnBitsInCommon(self):
    """ test OnBitsInCommon
    """
    v1 = klass(10)
    v1[0] = 1
    v1[2] = 1
    v1[6] = 1
    v2 = klass(10)
    v2[0] = 1
    v2[3] = 1
    v2[6] = 1
    v3 = cDataStructs.OnBitsInCommon(v1,v2)
    assert tuple(v3)==(0,6),'bad on bits in common'
  def testOffBitsInCommon(self):
    """ test OffBitsInCommon
    """
    v1 = klass(10)
    v1[0] = 1
    v1[2] = 1
    v1[6] = 1
    v2 = klass(10)
    v2[0] = 1
    v2[3] = 1
    v2[6] = 1
    v3 = cDataStructs.OffBitsInCommon(v1,v2)
    assert tuple(v3)==(1,4,5,7,8,9),'bad off bits in common'

  def testOnBitProjSimilarity(self):
    """ test OnBitProjSimilarity
    """
    v1 = klass(10)
    v1[1] = 1
    v1[2] = 1
    v1[3] = 1
    v2 = klass(10)
    v2[2] = 1
    v2[3] = 1
    res = cDataStructs.OnBitProjSimilarity(v1,v2)
    assert feq(res[0],0.666667),'bad 1st OnBitsProjSimilarity'
    assert feq(res[1],1.0),'bad 2nd OnBitsProjSimilarity'
    res = cDataStructs.OnBitProjSimilarity(v2,v1)
    assert feq(res[1],0.666667),'bad 1st OnBitsProjSimilarity'
    assert feq(res[0],1.0),'bad 2nd OnBitsProjSimilarity'
    
  def testOffBitProjSimilarity(self):
    """ test OffBitProjSimilarity
    """
    v1 = klass(10)
    v1[1] = 1
    v1[2] = 1
    v1[3] = 1
    v2 = klass(10)
    v2[2] = 1
    v2[3] = 1
    res = cDataStructs.OffBitProjSimilarity(v1,v2)
    assert feq(res[0],1.0),'bad 1st OffBitsProjSimilarity'
    assert feq(res[1],0.875),'bad 2nd OffBitsProjSimilarity'
    res = cDataStructs.OffBitProjSimilarity(v2,v1)
    assert feq(res[1],1.0),'bad 1st OffBitsProjSimilarity'
    assert feq(res[0],0.875),'bad 2nd OffBitsProjSimilarity'
    
  def testPkl(self):
    """ test pickling 
    """
    v1 = klass(10)
    v1[1] = 1
    v1[2] = 1
    v1[3] = 1
    pklName = 'foo.pkl'
    outF = open(pklName,'wb+')
    cPickle.dump(v1,outF)
    outF.close()
    inF = open(pklName,'rb')
    v2 = cPickle.load(inF)
    inF.close()
    os.unlink(pklName)
    assert tuple(v1.GetOnBits())==tuple(v2.GetOnBits()),'pkl failed'

  def testFingerprints(self):
    " test the parsing of daylight fingerprints "
    #actual daylight output:
    rawD="""
0,Cc1n[nH]c(=O)nc1N,.b+HHa.EgU6+ibEIr89.CpX0g8FZiXH+R0+Ps.mr6tg.2
1,Cc1n[nH]c(=O)[nH]c1=O,.b7HEa..ccc+gWEIr89.8lV8gOF3aXFFR.+Ps.mZ6lg.2
2,Cc1nnc(NN)nc1O,.H+nHq2EcY09y5EIr9e.8p50h0NgiWGNx4+Hm+Gbslw.2
3,Cc1nnc(N)nc1C,.1.HHa..cUI6i5E2rO8.Op10d0NoiWGVx.+Hm.Gb6lo.2
"""
    dists="""0,0,1.000000
0,1,0.788991
0,2,0.677165
0,3,0.686957
1,1,1.000000
1,2,0.578125
1,3,0.591304
2,2,1.000000
2,3,0.732759
3,3,1.000000
"""

    fps = []
    for line in rawD.split('\n'):
      if line:
        sbv = klass(256)
        id,smi,fp=line.split(',')
        cDataStructs.InitFromDaylightString(sbv,fp)
        fps.append(sbv)

    ds = dists.split('\n')
    whichd=0
    for i in range(len(fps)):
      for j in range(i,len(fps)):
        idx1,idx2,tgt = ds[whichd].split(',')
        whichd += 1
        tgt = float(tgt)
        dist = cDataStructs.TanimotoSimilarity(fps[i],fps[j])
        assert feq(tgt,dist),'tanimoto between fps %d and %d failed'%(int(idx1),int(idx2))

  def testFold(self):
    """ test folding fingerprints
    """
    v1 = klass(16)
    v1[1] = 1
    v1[12] = 1
    v1[9] = 1
    try:
      v2 = cDataStructs.FoldFingerprint(v1)
    except:
      assert 0,'Fold with no args failed'
    assert v1.GetNumBits()/2==v2.GetNumBits(),'bad num bits post folding'
    try:
      v2 = cDataStructs.FoldFingerprint(v1,2)
    except:
      assert 0,'Fold with arg failed'
      
    assert v1.GetNumBits()/2==v2.GetNumBits(),'bad num bits post folding'


    v2 = cDataStructs.FoldFingerprint(v1,4)
    assert v1.GetNumBits()/4==v2.GetNumBits(),'bad num bits post folding'

  def testOtherSims(self):
    """ test other similarity measures
    """
    v1 = klass(10)
    v1[0] = 1
    v1[2] = 1
    v1[6] = 1
    v2 = klass(10)
    v2[0] = 1
    v2[3] = 1
    v2[6] = 1
    assert feq(cDataStructs.CosineSimilarity(v1,v2),.6667)
    assert feq(cDataStructs.KulczynskiSimilarity(v1,v2),.6667)
    assert feq(cDataStructs.DiceSimilarity(v1,v2),.6667)
    assert feq(cDataStructs.SokalSimilarity(v1,v2),.3333)
    assert feq(cDataStructs.McConnaugheySimilarity(v1,v2),.3333)
    assert feq(cDataStructs.AsymmetricSimilarity(v1,v2),.6667)
    assert feq(cDataStructs.BraunBlanquetSimilarity(v1,v2),.6667)
    assert feq(cDataStructs.RusselSimilarity(v1,v2),.2000)
    assert feq(cDataStructs.RogotGoldbergSimilarity(v1,v2),.7619)


  def testQuickSims(self):
    """ the asymmetric similarity stuff (bv,pkl)
    """
    v1 = klass(10)
    v1[0] = 1
    v1[2] = 1
    v1[6] = 1
    v2 = klass(10)
    v2[0] = 1
    v2[3] = 1
    v2[6] = 1
    pkl = v2.ToBinary()
    v2 = pkl
    assert feq(cDataStructs.CosineSimilarity(v1,v2),.6667)
    assert feq(cDataStructs.KulczynskiSimilarity(v1,v2),.6667)
    assert feq(cDataStructs.DiceSimilarity(v1,v2),.6667)
    assert feq(cDataStructs.SokalSimilarity(v1,v2),.3333)
    assert feq(cDataStructs.McConnaugheySimilarity(v1,v2),.3333)
    assert feq(cDataStructs.AsymmetricSimilarity(v1,v2),.6667)
    assert feq(cDataStructs.BraunBlanquetSimilarity(v1,v2),.6667)
    assert feq(cDataStructs.RusselSimilarity(v1,v2),.2000)
    assert feq(cDataStructs.RogotGoldbergSimilarity(v1,v2),.7619)


if __name__ == '__main__':
  unittest.main()
