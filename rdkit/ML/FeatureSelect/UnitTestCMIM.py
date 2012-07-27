try:
  from rdkit.ML.FeatureSelect import CMIM
except ImportError:
  CMIM = None
from rdkit import DataStructs as DS
from rdkit import RDConfig
import unittest

class TestCase(unittest.TestCase):
   def setUp(self) :
      pass

   def test0FromList(self) :
     if CMIM is None:
       return
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

     r = CMIM.SelectFeatures(examples,2)
     self.failUnless(r==(2,4))

     # here we ask for three features, but there are only two
     # that are non-redundant:
     r = CMIM.SelectFeatures(examples,3)
     self.failUnless(r==(2,4))

if __name__ == '__main__':
   unittest.main()
