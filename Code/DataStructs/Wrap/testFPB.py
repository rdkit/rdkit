from rdkit import DataStructs
from rdkit import RDConfig
import unittest,os

def feq(a,b,tol=1e-4):
  return abs(a-b)<tol

class TestCase(unittest.TestCase):
   def setUp(self) :
      self.filename = os.path.join(RDConfig.RDBaseDir,'Code','DataStructs','testData','zim.head100.fpb')
      self.fpbr = DataStructs.FPBReader(self.filename)
      self.fpbr.Init()
   def test1Basics(self) :
      self.assertEqual(len(self.fpbr),100)
      self.assertEqual(self.fpbr.GetId(0),"ZINC00902219")
      self.assertEqual(self.fpbr.GetId(3),"ZINC04803506")

      fp = self.fpbr.GetFP(0)
      self.assertEqual(fp.GetNumBits(),2048)
      self.assertEqual(fp.GetNumOnBits(),17)
      obs = (1,   80,  183, 222,  227,  231,  482,  650, 807,
                        811, 831, 888, 1335, 1411, 1664, 1820, 1917)
      obl = tuple(fp.GetOnBits())
      self.assertEqual(obs,obl)


if __name__ == '__main__':
   unittest.main()
