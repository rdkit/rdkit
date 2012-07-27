from rdkit.ML import FeatureSelect as FS
from rdkit import DataStructs as DS
from rdkit import RDConfig
import unittest

class TestCase(unittest.TestCase):
   def setUp(self) :
      pass

   def test0FromList(self) :
     examples = []

     bv = DS.ExplicitBitVect(5)
     bv.SetBitsFromList([0,2,4])
     examples.append([0,bv,0])

     bv = DS.ExplicitBitVect(5)
     bv.SetBitsFromList([0,2,4])
     examples.append([0,bv,0])

     bv = DS.ExplicitBitVect(5)
     bv.SetBitsFromList([0,3,4])
     examples.append([0,bv,1])

     bv = DS.ExplicitBitVect(5)
     bv.SetBitsFromList([0,2,4])
     examples.append([0,bv,0])
     
     bv = DS.ExplicitBitVect(5)
     bv.SetBitsFromList([0,2])
     examples.append([0,bv,1])

     r = FS.selectCMIM(examples,2)
     self.failUnlessEqual(r,(2,4))

     r = FS.selectCMIM(examples,1)
     self.failUnlessEqual(r,(2,))

     r = FS.selectCMIM(examples,3)
     self.failUnlessEqual(r,(2,4,-1))

   def test1FromList(self) :
     examples = []

     bv = DS.SparseBitVect(5)
     bv.SetBitsFromList([0,2,4])
     examples.append([0,bv,0])

     bv = DS.SparseBitVect(5)
     bv.SetBitsFromList([0,2,4])
     examples.append([0,bv,0])

     bv = DS.SparseBitVect(5)
     bv.SetBitsFromList([0,3,4])
     examples.append([0,bv,1])

     bv = DS.SparseBitVect(5)
     bv.SetBitsFromList([0,2,4])
     examples.append([0,bv,0])
     
     bv = DS.SparseBitVect(5)
     bv.SetBitsFromList([0,2])
     examples.append([0,bv,1])

     r = FS.selectCMIM(examples,2)
     self.failUnlessEqual(r,(2,4))

     r = FS.selectCMIM(examples,1)
     self.failUnlessEqual(r,(2,))

     r = FS.selectCMIM(examples,3)
     self.failUnlessEqual(r,(2,4,-1))


if __name__ == '__main__':
   unittest.main()
