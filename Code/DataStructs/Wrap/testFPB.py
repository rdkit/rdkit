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

   def test3Tversky(self) :
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

   def test4Contains(self):
       bv = self.fpbr.GetBytes(0)
       nbrs = self.fpbr.GetContainingNeighbors(bv)
       self.assertEqual(len(nbrs),1)
       self.assertEqual(nbrs[0],0)

       bv = self.fpbr.GetBytes(1)
       nbrs = self.fpbr.GetContainingNeighbors(bv)
       self.assertEqual(len(nbrs),4)
       self.assertEqual(nbrs,(1,2,3,4))

   def test5Contains(self):
       " an example based on substructure screening "
       filename = os.path.join(RDConfig.RDBaseDir,'Code','DataStructs','testData','zinc_all_clean.100.patt1k.fpb')
       fpbr = DataStructs.FPBReader(filename)
       fpbr.Init()
       # these are the pattern bytes for "c1cncnc1", generated like this:
       # DataStructs.BitVectToFPSText(Chem.PatternFingerprint(Chem.MolFromSmiles('c1cncnc1'),1024))
       bytes = b'\x00\x00\x00\x00\x00\x00@\x00\x00\x00\x00\x00\x00\x00\x00\x00\x000\x00@\x00 \x00\x00 \x00\x00\x02@\x00\x00\x00\x00\x00\x00\x80\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00`\x07\x00\x04\x00"\x14\x02\x00\x00"\x00\x00\x00\x00\x08\x00\x80\x00\x00@\x00@\x00\x80\x00\x00\x00\x00B\x00\x00\x80\x00\x80\x08\x00\x04\x00@\x00\x00\x00\x00\x10\x00\x00\x00\x00\x00  \x00\x00\x00\x00\x00\x00\x00\x08\x00\x00\x00\x80\x04\x00\x00\x0c\x00\x00\x00@\x88\x10\x10\x00\x00\x88\x00@'
       nbrs = fpbr.GetContainingNeighbors(bytes)
       self.assertEqual(len(nbrs),9)
       ids = sorted(fpbr.GetId(x) for x in nbrs)
       self.assertEqual(ids,['ZINC00000562',
        'ZINC00000843',
        'ZINC00000969',
        'ZINC00001484',
        'ZINC00001585',
        'ZINC00002094',
        'ZINC00004739',
        'ZINC00005235',
        'ZINC00006300'])


if __name__ == '__main__':
   unittest.main()
