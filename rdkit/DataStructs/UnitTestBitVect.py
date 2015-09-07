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
""" unit testing code for the python BitVects
"""

import unittest
from rdkit.six.moves import cPickle
from rdkit.DataStructs import BitVect

def feq(v1,v2,tol=1e-4):
  return abs(v1-v2)<tol

class TestCase(unittest.TestCase):
  def _test_set_valid_bits(self, v):
    v[0] = 1
    v[2] = 1
    self.assertEquals(v[0], 1, "bad bit")
    self.assertEquals(v[1], 0, "bad bit")
    self.assertEquals(v[2], 1, "bad bit")
    self.assertEquals(v[3], 0, "bad bit")

  def _test_out_of_range_bits(self, v):
    with self.assertRaisesRegexp(Exception, "bad index"):
      v[10] = 1
    with self.assertRaisesRegexp(Exception, "bad index"):
      v[-1] = 1

    with self.assertRaisesRegexp(Exception, "bad index"):
      foo = v[10]
    with self.assertRaisesRegexp(Exception, "bad index"):
      foo = v[-1]

  def testBitVect(self):
    v = BitVect.BitVect(10)
    self._test_set_valid_bits(v)
    self._test_out_of_range_bits(v)
    
  def testSparseIdx(self):
    v = BitVect.SparseBitVect(10)
    self._test_set_valid_bits(v)
    self._test_out_of_range_bits(v)

  def _test_bit_get(self, v):
    v[0] = 1
    v[2] = 1
    v[6] = 1
    self.assertEqual(len(v), 10)
    self.assertEqual(v.NumOnBits(), 3)
    self.assertEqual(v.GetOnBits(), [0,2,6])
    self.assertEqual(v.GetOnBits(reverse=1), [6,2,0])
    
  def testVectBitGet(self):
    v = BitVect.BitVect(10)
    self._test_bit_get(v)
    
  def testSparseBitGet(self):
    v = BitVect.SparseBitVect(10)
    self._test_bit_get(v)

  def _test_bit_ops(self, v1, v2):
    v1[0] = 1
    v1[2] = 1
    v1[6] = 1
    v2[0] = 1
    v2[3] = 1
    v2[6] = 1
    self.assertEqual(v1&v2, [0,6],'binary & failed')
    self.assertEqual(v1|v2, [0,2,3,6],'binary | failed')
    self.assertEqual(v1^v2 , [2,3],'binary ^ failed')
    
  def testVectBitOps(self):
    v1 = BitVect.BitVect(10)
    v2 = BitVect.BitVect(10)
    self._test_bit_ops(v1, v2)

  def testCrossBitOps(self):
    v1 = BitVect.SparseBitVect(10)
    v2 = BitVect.BitVect(10)
    self._test_bit_ops(v1, v2)
    
    v1 = BitVect.BitVect(10)
    v2 = BitVect.SparseBitVect(10)
    self._test_bit_ops(v1, v2)
    
  def testSparseBitOps(self):
    v1 = BitVect.SparseBitVect(10)
    v2 = BitVect.SparseBitVect(10)
    self._test_bit_ops(v1, v2)

  def _test_tanimoto(self, v1, v2):
    v1[0] = 1
    v1[2] = 1
    v1[6] = 1
    v2[0] = 1
    v2[3] = 1
    v2[6] = 1
    self.assertEqual(v1.TanimotoSimilarity(v2), 0.5, 'TanimotoSimilarity failed')
    self.assertEqual(v2.TanimotoSimilarity(v1), 0.5, 'TanimotoSimilarity failed')
    self.assertEqual(v1.TanimotoSimilarity(v1), 1.0, 'TanimotoSimilarity failed')
    self.assertEqual(v2.TanimotoSimilarity(v2), 1.0, 'TanimotoSimilarity failed')
    
  def testVectTanimoto(self):
    v1 = BitVect.BitVect(10)
    v2 = BitVect.BitVect(10)
    self._test_tanimoto(v1, v2)
    
  def testSparseTanimoto(self):
    v1 = BitVect.SparseBitVect(10)
    v2 = BitVect.SparseBitVect(10)
    self._test_tanimoto(v1, v2)
    
  def testCrossTanimoto(self):
    v1 = BitVect.BitVect(10)
    v2 = BitVect.SparseBitVect(10)
    self._test_tanimoto(v1, v2)
    
    v1 = BitVect.SparseBitVect(10)
    v2 = BitVect.BitVect(10)
    self._test_tanimoto(v1, v2)


  def _test_Euclid(self, v1, v2):
    v1[0] = 1
    v1[2] = 1
    v1[6] = 1
    v2[0] = 1
    v2[3] = 1
    v2[6] = 1
    self.assertLess(abs(v1.EuclideanDistance(v2)-.14142), .0001, 'EuclideanDistance failed')
    self.assertLess(abs(v2.EuclideanDistance(v1)-.14142), .0001, 'EuclideanDistance failed')
    self.assertEqual(v1.EuclideanDistance(v1), 0.0, 'EuclideanDistance failed')
    self.assertEqual(v2.EuclideanDistance(v2), 0.0, 'EuclideanDistance failed')
    
  def testVectEuclid(self):
    v1 = BitVect.BitVect(10)
    v2 = BitVect.BitVect(10)
    self._test_Euclid(v1, v2)

  def testSparseEuclid(self):
    v1 = BitVect.SparseBitVect(10)
    v2 = BitVect.SparseBitVect(10)
    self._test_Euclid(v1, v2)

  def testCrossEuclid(self):
    v1 = BitVect.BitVect(10)
    v2 = BitVect.SparseBitVect(10)
    self._test_Euclid(v1, v2)

    v1 = BitVect.SparseBitVect(10)
    v2 = BitVect.BitVect(10)
    self._test_Euclid(v1, v2)
    
  def test90BulkDistances(self):
    """
      verify that the base similarity (tanimoto) works using an 18 fp regression
      panel of pubchem compounds (chosen to have different lengths: 5x2048,
      5x1024, 5x512, 3x256)
    """  
    from rdkit import DataStructs
    import os
    from rdkit import RDConfig
    fps = cPickle.load(file(os.path.join(RDConfig.RDCodeDir,'DataStructs','test_data',
                                         'pubchem_fps.pkl'),'rb'))
    dm = cPickle.load(file(os.path.join(RDConfig.RDCodeDir,'DataStructs','test_data',
                                         'pubchem_fps.dm.pkl'),'rb'))
    dmIdx=0
    for i in range(len(fps)):
      nmi,fpi = fps[i]
      for j in range(i+1,len(fps)):
        nmj,fpj = fps[j]
        sim = DataStructs.FingerprintSimilarity(fpi,fpj)
        self.failUnless(feq(sim,dm[dmIdx]))
        dmIdx+=1
        
  def test91BoundDistances(self):
    """
      verify that the bounded similarity (tanimoto) works
    """  
    from rdkit import DataStructs
    import os
    from rdkit import RDConfig
    fps = cPickle.load(file(os.path.join(RDConfig.RDCodeDir,'DataStructs','test_data',
                                         'pubchem_fps.pkl'),'rb'))
    dm = cPickle.load(file(os.path.join(RDConfig.RDCodeDir,'DataStructs','test_data',
                                         'pubchem_fps.dm.pkl'),'rb'))
    dmIdx=0
    for i in range(len(fps)):
      nmi,fpi = fps[i]
      for j in range(i+1,len(fps)):
        nmj,fpj = fps[j]
        sim = DataStructs.FingerprintSimilarity(fpi,fpj)
        self.failUnless(feq(sim,dm[dmIdx]))
        dmIdx+=1
        
if __name__ == '__main__':
  unittest.main()
