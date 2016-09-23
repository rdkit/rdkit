from rdkit import DataStructs
from rdkit import RDConfig
import unittest, os


def feq(a, b, tol=1e-4):
  return abs(a - b) < tol


class TestCase(unittest.TestCase):

  def setUp(self):
    self.dirname = os.path.join(RDConfig.RDBaseDir, 'Code', 'DataStructs', 'testData')
    self.filename = os.path.join(self.dirname, 'zim.head100.fpb')
    self.fpbr = DataStructs.FPBReader(self.filename)
    self.fpbr.Init()

  def test1Basics(self):
    self.assertEqual(len(self.fpbr), 100)
    self.assertEqual(self.fpbr.GetNumBits(), 2048)
    self.assertEqual(self.fpbr.GetId(0), "ZINC00902219")
    self.assertEqual(self.fpbr.GetId(3), "ZINC04803506")

    fp = self.fpbr.GetFP(0)
    self.assertEqual(fp.GetNumBits(), 2048)
    self.assertEqual(fp.GetNumOnBits(), 17)
    obs = (1, 80, 183, 222, 227, 231, 482, 650, 807, 811, 831, 888, 1335, 1411, 1664, 1820, 1917)
    obl = tuple(fp.GetOnBits())
    self.assertEqual(obs, obl)

    # test operator[]
    fp, nm = self.fpbr[0]
    self.assertEqual(nm, "ZINC00902219")
    self.assertEqual(fp.GetNumOnBits(), 17)

  def test2Tanimoto(self):
    bv = self.fpbr.GetBytes(0)
    self.assertAlmostEqual(self.fpbr.GetTanimoto(0, bv), 1.0, 4)
    self.assertAlmostEqual(self.fpbr.GetTanimoto(1, bv), 0.3704, 4)
    tpl = self.fpbr.GetTanimotoNeighbors(bv)
    self.assertEqual(len(tpl), 1)
    self.assertEqual(tpl[0][1], 0)
    self.assertAlmostEqual(tpl[0][0], 1., 4)
    tpl = self.fpbr.GetTanimotoNeighbors(bv, threshold=0.3)
    self.assertEqual(len(tpl), 5)
    self.assertEqual(tpl[0][1], 0)
    self.assertAlmostEqual(tpl[0][0], 1., 4)
    self.assertEqual(tpl[1][1], 1)
    self.assertAlmostEqual(tpl[1][0], 0.3704, 4)

  def test3Tversky(self):
    bv = self.fpbr.GetBytes(0)
    self.assertAlmostEqual(self.fpbr.GetTversky(0, bv, 1, 1), 1.0, 4)
    self.assertAlmostEqual(self.fpbr.GetTversky(1, bv, 1, 1), 0.3704, 4)
    tpl = self.fpbr.GetTverskyNeighbors(bv, 1, 1)
    self.assertEqual(len(tpl), 1)
    self.assertEqual(tpl[0][1], 0)
    self.assertAlmostEqual(tpl[0][0], 1., 4)
    tpl = self.fpbr.GetTverskyNeighbors(bv, 1, 1, threshold=0.3)
    self.assertEqual(len(tpl), 5)
    self.assertEqual(tpl[0][1], 0)
    self.assertAlmostEqual(tpl[0][0], 1., 4)
    self.assertEqual(tpl[1][1], 1)
    self.assertAlmostEqual(tpl[1][0], 0.3704, 4)

  def test4Contains(self):
    bv = self.fpbr.GetBytes(0)
    nbrs = self.fpbr.GetContainingNeighbors(bv)
    self.assertEqual(len(nbrs), 1)
    self.assertEqual(nbrs[0], 0)

    bv = self.fpbr.GetBytes(1)
    nbrs = self.fpbr.GetContainingNeighbors(bv)
    self.assertEqual(len(nbrs), 4)
    self.assertEqual(nbrs, (1, 2, 3, 4))

  def test5Contains(self):
    " an example based on substructure screening "
    filename = os.path.join(RDConfig.RDBaseDir, 'Code', 'DataStructs', 'testData',
                            'zinc_all_clean.100.patt1k.fpb')
    fpbr = DataStructs.FPBReader(filename)
    fpbr.Init()
    bytes = b'\x00\x00\x00\x00\x00\x00@\x00\x00\x00\x00\x00\x00\x00\x00\x00\x000\x00@\x00 \x00\x00 \x00\x00\x02@\x00\x00\x00\x00\x00\x00\x80\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00`\x07\x00\x04\x00"\x14\x02\x00\x00"\x00\x00\x00\x00\x08\x00\x80\x00\x00@\x00@\x00\x80\x00\x00\x00\x00B\x00\x00\x80\x00\x80\x08\x00\x04\x00@\x00\x00\x00\x00\x10\x00\x00\x00\x00\x00  \x00\x00\x00\x00\x00\x00\x00\x08\x00\x00\x00\x80\x04\x00\x00\x0c\x00\x00\x00@\x88\x10\x10\x00\x00\x88\x00@'
    nbrs = fpbr.GetContainingNeighbors(bytes)
    self.assertEqual(len(nbrs), 9)
    ids = sorted(fpbr.GetId(x) for x in nbrs)
    self.assertEqual(ids, ['ZINC00000562', 'ZINC00000843', 'ZINC00000969', 'ZINC00001484',
                           'ZINC00001585', 'ZINC00002094', 'ZINC00004739', 'ZINC00005235',
                           'ZINC00006300'])

  def test6MultiFPBReaderTani(self):
    basen = os.path.join(RDConfig.RDBaseDir, 'Code', 'DataStructs', 'testData')
    mfpbr = DataStructs.MultiFPBReader()
    self.assertEqual(
      mfpbr.AddReader(DataStructs.FPBReader(os.path.join(basen, "zinc_random200.1.patt.fpb"))), 1)
    self.assertEqual(
      mfpbr.AddReader(DataStructs.FPBReader(os.path.join(basen, "zinc_random200.2.patt.fpb"))), 2)
    self.assertEqual(
      mfpbr.AddReader(DataStructs.FPBReader(os.path.join(basen, "zinc_random200.3.patt.fpb"))), 3)
    self.assertEqual(
      mfpbr.AddReader(DataStructs.FPBReader(os.path.join(basen, "zinc_random200.4.patt.fpb"))), 4)
    mfpbr.Init()
    self.assertEqual(mfpbr.GetNumBits(), 1024)
    self.assertEqual(len(mfpbr), 4)

    fps = "0000000000404000100000001000040000300040222000002004000240000020000000"+\
"8200010200000090000024040860070044003214820000220401054008018000226000"+\
"4800800140000042000080008008020482400000200410800000300430200800400000"+\
"0000080a0000800400010c800200648818100010880040"
    ebv = DataStructs.CreateFromFPSText(fps)
    bytes = DataStructs.BitVectToBinaryText(ebv)
    nbrs = mfpbr.GetTanimotoNeighbors(bytes, threshold=0.6)
    self.assertEqual(len(nbrs), 6)
    self.assertAlmostEqual(nbrs[0][0], 0.66412, 4)
    self.assertEqual(nbrs[0][1], 0)
    self.assertEqual(nbrs[0][2], 3)
    self.assertAlmostEqual(nbrs[1][0], 0.65289, 4)
    self.assertEqual(nbrs[1][1], 1)
    self.assertEqual(nbrs[1][2], 2)
    self.assertAlmostEqual(nbrs[2][0], 0.64341, 4)
    self.assertEqual(nbrs[2][1], 2)
    self.assertEqual(nbrs[2][2], 1)
    self.assertAlmostEqual(nbrs[3][0], 0.61940, 4)
    self.assertEqual(nbrs[3][1], 1)
    self.assertEqual(nbrs[3][2], 0)
    self.assertAlmostEqual(nbrs[4][0], 0.61905, 4)
    self.assertEqual(nbrs[4][1], 0)
    self.assertEqual(nbrs[4][2], 0)
    self.assertAlmostEqual(nbrs[5][0], 0.61344, 4)
    self.assertEqual(nbrs[5][1], 0)
    self.assertEqual(nbrs[5][2], 1)

    # test multi-threaded (won't do anything if the RDKit isn't compiled with threads support)
    nbrs = mfpbr.GetTanimotoNeighbors(bytes, threshold=0.6, numThreads=4)
    self.assertEqual(len(nbrs), 6)
    self.assertAlmostEqual(nbrs[0][0], 0.66412, 4)
    self.assertEqual(nbrs[0][1], 0)
    self.assertEqual(nbrs[0][2], 3)
    self.assertAlmostEqual(nbrs[1][0], 0.65289, 4)
    self.assertEqual(nbrs[1][1], 1)
    self.assertEqual(nbrs[1][2], 2)
    self.assertAlmostEqual(nbrs[2][0], 0.64341, 4)
    self.assertEqual(nbrs[2][1], 2)
    self.assertEqual(nbrs[2][2], 1)
    self.assertAlmostEqual(nbrs[3][0], 0.61940, 4)
    self.assertEqual(nbrs[3][1], 1)
    self.assertEqual(nbrs[3][2], 0)
    self.assertAlmostEqual(nbrs[4][0], 0.61905, 4)
    self.assertEqual(nbrs[4][1], 0)
    self.assertEqual(nbrs[4][2], 0)
    self.assertAlmostEqual(nbrs[5][0], 0.61344, 4)
    self.assertEqual(nbrs[5][1], 0)
    self.assertEqual(nbrs[5][2], 1)

  def test7MultiFPBReaderContains(self):
    basen = os.path.join(RDConfig.RDBaseDir, 'Code', 'DataStructs', 'testData')
    mfpbr = DataStructs.MultiFPBReader()
    self.assertEqual(
      mfpbr.AddReader(DataStructs.FPBReader(os.path.join(basen, "zinc_random200.1.patt.fpb"))), 1)
    self.assertEqual(
      mfpbr.AddReader(DataStructs.FPBReader(os.path.join(basen, "zinc_random200.2.patt.fpb"))), 2)
    self.assertEqual(
      mfpbr.AddReader(DataStructs.FPBReader(os.path.join(basen, "zinc_random200.3.patt.fpb"))), 3)
    self.assertEqual(
      mfpbr.AddReader(DataStructs.FPBReader(os.path.join(basen, "zinc_random200.4.patt.fpb"))), 4)
    mfpbr.Init()
    self.assertEqual(mfpbr.GetNumBits(), 1024)
    self.assertEqual(len(mfpbr), 4)

    fps = "40081010824820021000500010110410003000402b20285000a4040240010030050000"+\
            "080001420040009000003d04086007080c03b31d920004220400074008098010206080"+\
            "00488001080000c64002a00080000200024c2000602410049200340820200002400010"+\
            "02200106090401056801080182006088101000088a0048"
    ebv = DataStructs.CreateFromFPSText(fps)
    bytes = DataStructs.BitVectToBinaryText(ebv)
    nbrs = mfpbr.GetContainingNeighbors(bytes)
    self.assertEqual(len(nbrs), 9)
    self.assertEqual(nbrs[0][0], 160)
    self.assertEqual(nbrs[0][1], 0)
    self.assertEqual(nbrs[1][0], 163)
    self.assertEqual(nbrs[1][1], 0)
    self.assertEqual(nbrs[2][0], 170)
    self.assertEqual(nbrs[2][1], 0)
    self.assertEqual(nbrs[3][0], 180)
    self.assertEqual(nbrs[3][1], 2)
    self.assertEqual(nbrs[4][0], 182)
    self.assertEqual(nbrs[4][1], 3)
    self.assertEqual(nbrs[5][0], 185)
    self.assertEqual(nbrs[5][1], 0)
    self.assertEqual(nbrs[6][0], 189)
    self.assertEqual(nbrs[6][1], 0)
    self.assertEqual(nbrs[7][0], 192)
    self.assertEqual(nbrs[7][1], 3)
    self.assertEqual(nbrs[8][0], 193)
    self.assertEqual(nbrs[8][1], 0)

    nbrs = mfpbr.GetContainingNeighbors(bytes, numThreads=4)
    self.assertEqual(len(nbrs), 9)
    self.assertEqual(nbrs[0][0], 160)
    self.assertEqual(nbrs[0][1], 0)
    self.assertEqual(nbrs[1][0], 163)
    self.assertEqual(nbrs[1][1], 0)
    self.assertEqual(nbrs[2][0], 170)
    self.assertEqual(nbrs[2][1], 0)
    self.assertEqual(nbrs[3][0], 180)
    self.assertEqual(nbrs[3][1], 2)
    self.assertEqual(nbrs[4][0], 182)
    self.assertEqual(nbrs[4][1], 3)
    self.assertEqual(nbrs[5][0], 185)
    self.assertEqual(nbrs[5][1], 0)
    self.assertEqual(nbrs[6][0], 189)
    self.assertEqual(nbrs[6][1], 0)
    self.assertEqual(nbrs[7][0], 192)
    self.assertEqual(nbrs[7][1], 3)
    self.assertEqual(nbrs[8][0], 193)
    self.assertEqual(nbrs[8][1], 0)

  def test8MultiFPBReaderContainsInitOnSearch(self):
    basen = os.path.join(RDConfig.RDBaseDir, 'Code', 'DataStructs', 'testData')
    mfpbr = DataStructs.MultiFPBReader(initOnSearch=True)
    self.assertEqual(
      mfpbr.AddReader(DataStructs.FPBReader(os.path.join(basen, "zinc_random200.1.patt.fpb"))), 1)
    self.assertEqual(
      mfpbr.AddReader(DataStructs.FPBReader(os.path.join(basen, "zinc_random200.2.patt.fpb"))), 2)
    self.assertEqual(
      mfpbr.AddReader(DataStructs.FPBReader(os.path.join(basen, "zinc_random200.3.patt.fpb"))), 3)
    self.assertEqual(
      mfpbr.AddReader(DataStructs.FPBReader(os.path.join(basen, "zinc_random200.4.patt.fpb"))), 4)

    fps = "40081010824820021000500010110410003000402b20285000a4040240010030050000"+\
            "080001420040009000003d04086007080c03b31d920004220400074008098010206080"+\
            "00488001080000c64002a00080000200024c2000602410049200340820200002400010"+\
            "02200106090401056801080182006088101000088a0048"
    ebv = DataStructs.CreateFromFPSText(fps)
    bytes = DataStructs.BitVectToBinaryText(ebv)
    nbrs = mfpbr.GetContainingNeighbors(bytes, numThreads=4)
    self.assertEqual(len(nbrs), 9)
    self.assertEqual(nbrs[0][0], 160)
    self.assertEqual(nbrs[0][1], 0)
    self.assertEqual(nbrs[1][0], 163)
    self.assertEqual(nbrs[1][1], 0)
    self.assertEqual(nbrs[2][0], 170)
    self.assertEqual(nbrs[2][1], 0)
    self.assertEqual(nbrs[3][0], 180)
    self.assertEqual(nbrs[3][1], 2)
    self.assertEqual(nbrs[4][0], 182)
    self.assertEqual(nbrs[4][1], 3)
    self.assertEqual(nbrs[5][0], 185)
    self.assertEqual(nbrs[5][1], 0)
    self.assertEqual(nbrs[6][0], 189)
    self.assertEqual(nbrs[6][1], 0)
    self.assertEqual(nbrs[7][0], 192)
    self.assertEqual(nbrs[7][1], 3)
    self.assertEqual(nbrs[8][0], 193)
    self.assertEqual(nbrs[8][1], 0)

  def test9MultiFPBReaderEdges(self):
    basen = os.path.join(RDConfig.RDBaseDir, 'Code', 'DataStructs', 'testData')
    mfpbr = DataStructs.MultiFPBReader()
    mfpbr.Init()

    fps = "0000000000404000100000001000040000300040222000002004000240000020000000"+\
"8200010200000090000024040860070044003214820000220401054008018000226000"+\
"4800800140000042000080008008020482400000200410800000300430200800400000"+\
"0000080a0000800400010c800200648818100010880040"
    ebv = DataStructs.CreateFromFPSText(fps)
    bytes = DataStructs.BitVectToBinaryText(ebv)
    nbrs = mfpbr.GetTanimotoNeighbors(bytes, threshold=0.6)
    self.assertEqual(len(nbrs), 0)


if __name__ == '__main__':
  unittest.main()
