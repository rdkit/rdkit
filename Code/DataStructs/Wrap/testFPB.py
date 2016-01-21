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
      self.assertEqual(self.fpbr.GetNumBits(),2048)
      self.assertEqual(self.fpbr.GetId(0),"ZINC00902219")
      self.assertEqual(self.fpbr.GetId(3),"ZINC04803506")

      fp = self.fpbr.GetFP(0)
      self.assertEqual(fp.GetNumBits(),2048)
      self.assertEqual(fp.GetNumOnBits(),17)
      obs = (1,   80,  183, 222,  227,  231,  482,  650, 807,
                        811, 831, 888, 1335, 1411, 1664, 1820, 1917)
      obl = tuple(fp.GetOnBits())
      self.assertEqual(obs,obl)

      # test operator[]
      fp,nm = self.fpbr[0]
      self.assertEqual(nm,"ZINC00902219")
      self.assertEqual(fp.GetNumOnBits(),17)


   def test2Tanimoto(self) :
      bv = self.fpbr.GetBytes(0)
      self.assertAlmostEqual(self.fpbr.GetTanimoto(0,bv),1.0,4)
      self.assertAlmostEqual(self.fpbr.GetTanimoto(1,bv),0.3704,4)
      tpl = self.fpbr.GetTanimotoNeighbors(bv)
      self.assertEqual(len(tpl),1)
      self.assertEqual(tpl[0][1],0)
      self.assertAlmostEqual(tpl[0][0],1.,4)
      tpl = self.fpbr.GetTanimotoNeighbors(bv,threshold=0.3)
      self.assertEqual(len(tpl),5)
      self.assertEqual(tpl[0][1],0)
      self.assertAlmostEqual(tpl[0][0],1.,4)
      self.assertEqual(tpl[1][1],1)
      self.assertAlmostEqual(tpl[1][0],0.3704,4)

   def test2Tversky(self) :
      bv = self.fpbr.GetBytes(0)
      self.assertAlmostEqual(self.fpbr.GetTversky(0,bv,1,1),1.0,4)
      self.assertAlmostEqual(self.fpbr.GetTversky(1,bv,1,1),0.3704,4)
      tpl = self.fpbr.GetTverskyNeighbors(bv,1,1)
      self.assertEqual(len(tpl),1)
      self.assertEqual(tpl[0][1],0)
      self.assertAlmostEqual(tpl[0][0],1.,4)
      tpl = self.fpbr.GetTverskyNeighbors(bv,1,1,threshold=0.3)
      self.assertEqual(len(tpl),5)
      self.assertEqual(tpl[0][1],0)
      self.assertAlmostEqual(tpl[0][0],1.,4)
      self.assertEqual(tpl[1][1],1)
      self.assertAlmostEqual(tpl[1][0],0.3704,4)


if __name__ == '__main__':
   unittest.main()
