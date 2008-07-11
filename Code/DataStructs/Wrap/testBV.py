import DataStructs
import RDConfig
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

if __name__ == '__main__':
   unittest.main()
