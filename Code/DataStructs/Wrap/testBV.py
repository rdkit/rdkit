from rdkit import DataStructs
from rdkit import RDConfig
import unittest
import cPickle as pickle
import random

def feq(a,b,tol=1e-4):
  return abs(a-b)<tol

class TestCase(unittest.TestCase):
   def setUp(self) :
      pass

   def test0FromList(self) :
      bv1 = DataStructs.SparseBitVect(1000)
      bv2 = DataStructs.SparseBitVect(1000)
      obits = range(0,1000, 3)
      
      for bit in obits :
         bv1.SetBit(bit)
         
      bv2.SetBitsFromList(obits)
         
      for i in range(1000) :
         assert bv1.GetBit(i) == bv2.GetBit(i)

      self.failUnless(bv1==bv2)
      bv2.SetBit(1)
      self.failUnless(bv1!=bv2)
      bv2.UnSetBit(1)
      self.failUnless(bv1==bv2)

      bv1 = DataStructs.ExplicitBitVect(1000)
      bv2 = DataStructs.ExplicitBitVect(1000)
      obits = range(0,1000, 3)

      for bit in obits :
         bv1.SetBit(bit)

      bv2.SetBitsFromList(obits)

      for i in range(1000) :
         assert bv1.GetBit(i) == bv2.GetBit(i)

   def test1SparsePickle(self) :
      nbits = 10000
      bv1 = DataStructs.SparseBitVect(nbits)
      for i in range(1000) :
         x = random.randrange(0,nbits)
         bv1.SetBit(x)
         
      pkl = pickle.dumps(bv1,1)
      bv2 = pickle.loads(pkl)
      for i in range(nbits) :
         assert bv1[i] == bv2[i]

   def test2ExplicitPickle(self):
      nbits = 10000
      bv1 = DataStructs.ExplicitBitVect(nbits)
      for i in range(1000) :
         x = random.randrange(0,nbits)
         bv1.SetBit(x)
      
      pkl = pickle.dumps(bv1,1)
      bv2 = pickle.loads(pkl)
      for i in range(nbits) :
         assert bv1[i] == bv2[i]

   def test3Bounds(self) :
      nbits = 10
      bv1 = DataStructs.ExplicitBitVect(nbits)
      bv1[0]
      try:
         bv1[11]
      except IndexError:
         ok = 1
      except:
         ok = -1
      else:
         ok = 0
      assert ok>0,ok  
      
   def test4OnBitsInCommon(self) :
      sz=100
      bv1 = DataStructs.ExplicitBitVect(sz)
      bv2 = DataStructs.ExplicitBitVect(sz)
      for i in range(0,sz,2):
         bv1.SetBit(i)
         if i < 3*sz/4:
            bv2.SetBit(i)
      self.failUnless(DataStructs.AllProbeBitsMatch(bv1,bv1.ToBinary()))
      self.failUnless(DataStructs.AllProbeBitsMatch(bv2,bv1.ToBinary()))
      self.failIf(DataStructs.AllProbeBitsMatch(bv1,bv2.ToBinary()))
      self.failUnless(DataStructs.AllProbeBitsMatch(bv2,bv2.ToBinary()))

   def test5FromBitString(self):
      s1 = '1010'
      bv = DataStructs.CreateFromBitString(s1)
      self.failUnless(len(bv)==4)
      self.failUnless(list(bv.GetOnBits())==[0,2])

   def test6BulkOps(self):
      nbits = 10000
      bvs = []
      for bvi in range(10):
        bv = DataStructs.ExplicitBitVect(nbits)
        for j in range(nbits) :
           x = random.randrange(0,nbits)
           bv.SetBit(x)
        bvs.append(bv)
      sims = DataStructs.BulkTanimotoSimilarity(bvs[0],bvs)
      for i in range(len(bvs)):
        sim = DataStructs.TanimotoSimilarity(bvs[0],bvs[i])
        self.failUnless(feq(sim,sims[i]))

      sims = DataStructs.BulkDiceSimilarity(bvs[0],bvs)
      for i in range(len(bvs)):
        sim = DataStructs.DiceSimilarity(bvs[0],bvs[i])
        self.failUnless(feq(sim,sims[i]))

      sims = DataStructs.BulkAllBitSimilarity(bvs[0],bvs)
      for i in range(len(bvs)):
        sim = DataStructs.AllBitSimilarity(bvs[0],bvs[i])
        self.failUnless(feq(sim,sims[i]))

      sims = DataStructs.BulkOnBitSimilarity(bvs[0],bvs)
      for i in range(len(bvs)):
        sim = DataStructs.OnBitSimilarity(bvs[0],bvs[i])
        self.failUnless(feq(sim,sims[i]))

      sims = DataStructs.BulkRogotGoldbergSimilarity(bvs[0],bvs)
      for i in range(len(bvs)):
        sim = DataStructs.RogotGoldbergSimilarity(bvs[0],bvs[i])
        self.failUnless(feq(sim,sims[i]))

      sims = DataStructs.BulkTverskySimilarity(bvs[0],bvs,1,1)
      for i in range(len(bvs)):
        sim = DataStructs.TverskySimilarity(bvs[0],bvs[i],1,1)
        self.failUnless(feq(sim,sims[i]))
        sim = DataStructs.TanimotoSimilarity(bvs[0],bvs[i])
        self.failUnless(feq(sim,sims[i]))

      sims = DataStructs.BulkTverskySimilarity(bvs[0],bvs,.5,.5)
      for i in range(len(bvs)):
        sim = DataStructs.TverskySimilarity(bvs[0],bvs[i],.5,.5)
        self.failUnless(feq(sim,sims[i]))
        sim = DataStructs.DiceSimilarity(bvs[0],bvs[i])
        self.failUnless(feq(sim,sims[i]))

   def test7FPS(self):
      bv = DataStructs.ExplicitBitVect(32)
      bv.SetBit(0)
      bv.SetBit(1)
      bv.SetBit(17)
      bv.SetBit(23)
      bv.SetBit(31)

      self.failUnlessEqual(DataStructs.BitVectToFPSText(bv),"03008280")
      bv2 = DataStructs.CreateFromFPSText("03008280")
      self.failUnlessEqual(bv,bv2)

      self.failUnlessRaises(ValueError,lambda : DataStructs.CreateFromFPSText("030082801"))
      
      bv2 = DataStructs.CreateFromFPSText("")
      self.failUnlessEqual(bv2.GetNumBits(),0)


   def test8BinText(self):
      bv = DataStructs.ExplicitBitVect(32)
      bv.SetBit(0)
      bv.SetBit(1)
      bv.SetBit(17)
      bv.SetBit(23)
      bv.SetBit(31)

      bv2 = DataStructs.CreateFromBinaryText(DataStructs.BitVectToBinaryText(bv))
      self.failUnlessEqual(bv,bv2)

      bv2 = DataStructs.CreateFromBinaryText("")
      self.failUnlessEqual(bv2.GetNumBits(),0)

   def test9ToNumpy(self):
      import numpy
      bv = DataStructs.ExplicitBitVect(32)
      bv.SetBit(0)
      bv.SetBit(1)
      bv.SetBit(17)
      bv.SetBit(23)
      bv.SetBit(31)
      arr = numpy.zeros((3,),'i')
      DataStructs.ConvertToNumpyArray(bv,arr)
      for i in range(bv.GetNumBits()):
        self.failUnlessEqual(bv[i],arr[i])


      
if __name__ == '__main__':
   unittest.main()
