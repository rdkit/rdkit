from rdkit import RDConfig
import unittest
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def test1(self):
    m = Chem.MolFromSmiles('c1ccc(C)c(C)c1C')
    AllChem.Compute2DCoords(m)
    d = Draw.MolDraw2DSVG(300, 300)
    d.DrawMolecule(m)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("<svg:svg") != -1)
    self.assertTrue(txt.find("</svg:svg>") != -1)

  def test2(self):
    m = Chem.MolFromSmiles('c1ccc(C)c(C)c1C')
    AllChem.Compute2DCoords(m)
    d = Draw.MolDraw2DSVG(300, 300)
    do = d.drawOptions()
    do.atomLabels[3] = 'foolabel'
    d.DrawMolecule(m)
    d.FinishDrawing()
    txt = d.GetDrawingText()
    self.assertTrue(txt.find("foolabel") != -1)

  def testGithubIssue571(self):
    if not hasattr(Draw, 'MolDraw2DCairo'):
      return
    m = Chem.MolFromSmiles('c1ccc(C)c(C)c1C')
    AllChem.Compute2DCoords(m)
    d = Draw.MolDraw2DCairo(300, 300)
    d.DrawMolecule(m)
    d.FinishDrawing()
    txt = d.GetDrawingText()

  def testPrepareForDrawing(self):
    m = Chem.MolFromSmiles('c1ccccc1[C@H](F)Cl')
    nm = rdMolDraw2D.PrepareMolForDrawing(m)
    self.assertEqual(nm.GetNumAtoms(), 9)
    self.assertEqual(nm.GetNumConformers(), 1)
    m = Chem.MolFromSmiles('C1CC[C@H]2NCCCC2C1')
    nm = rdMolDraw2D.PrepareMolForDrawing(m)
    self.assertEqual(nm.GetNumAtoms(), 11)
    self.assertEqual(nm.GetNumConformers(), 1)
    nm = rdMolDraw2D.PrepareMolForDrawing(m, addChiralHs=False)
    self.assertEqual(nm.GetNumAtoms(), 10)
    self.assertEqual(nm.GetNumConformers(), 1)

  def testRepeatedPrepareForDrawingCalls(self):
    m = Chem.MolFromMolBlock("""
          11280715312D 1   1.00000     0.00000     0

 33 36  0     1  0            999 V2000
    7.6125   -5.7917    0.0000 C   0  0  0  0  0  0           0  0  0
    7.0917   -6.0917    0.0000 C   0  0  1  0  0  0           0  0  0
    6.4792   -6.8917    0.0000 C   0  0  2  0  0  0           0  0  0
    8.1292   -6.0792    0.0000 N   0  0  0  0  0  0           0  0  0
    5.5042   -6.8917    0.0000 C   0  0  3  0  0  0           0  0  0
   11.2375   -4.8542    0.0000 N   0  0  0  0  0  0           0  0  0
    9.6792   -5.1667    0.0000 N   0  0  3  0  0  0           0  0  0
    5.9917   -6.5417    0.0000 C   0  0  0  0  0  0           0  0  0
    7.6042   -5.1917    0.0000 O   0  0  0  0  0  0           0  0  0
   10.7167   -5.1625    0.0000 C   0  0  0  0  0  0           0  0  0
    6.2917   -7.4667    0.0000 C   0  0  0  0  0  0           0  0  0
    6.5750   -5.7917    0.0000 C   0  0  0  0  0  0           0  0  0
   10.2000   -4.8667    0.0000 C   0  0  0  0  0  0           0  0  0
    8.6500   -5.7792    0.0000 C   0  0  3  0  0  0           0  0  0
    8.6417   -5.1792    0.0000 C   0  0  0  0  0  0           0  0  0
    9.1667   -6.0750    0.0000 C   0  0  0  0  0  0           0  0  0
    9.6875   -5.7667    0.0000 C   0  0  0  0  0  0           0  0  0
    9.1542   -4.8750    0.0000 C   0  0  0  0  0  0           0  0  0
    5.6917   -7.4667    0.0000 C   0  0  0  0  0  0           0  0  0
    5.2042   -7.4042    0.0000 F   0  0  0  0  0  0           0  0  0
    4.9875   -6.5917    0.0000 F   0  0  0  0  0  0           0  0  0
    7.5167   -6.5167    0.0000 O   0  0  0  0  0  0           0  0  0
   11.7542   -5.1500    0.0000 C   0  0  0  0  0  0           0  0  0
   11.2417   -6.0542    0.0000 C   0  0  0  0  0  0           0  0  0
   10.7250   -5.7625    0.0000 C   0  0  0  0  0  0           0  0  0
    6.5750   -5.1917    0.0000 C   0  0  0  0  0  0           0  0  0
    6.0542   -6.0917    0.0000 C   0  0  0  0  0  0           0  0  0
   11.7667   -5.7542    0.0000 C   0  0  0  0  0  0           0  0  0
   12.2750   -4.8417    0.0000 C   0  0  0  0  0  0           0  0  0
    6.0542   -4.8917    0.0000 C   0  0  0  0  0  0           0  0  0
    5.5375   -5.7917    0.0000 C   0  0  0  0  0  0           0  0  0
    5.5375   -5.1917    0.0000 C   0  0  0  0  0  0           0  0  0
    6.3167   -6.3042    0.0000 H   0  0  0  0  0  0           0  0  0
  2  1  1  0     0  0
  3  2  1  0     0  0
  4  1  1  0     0  0
  5  8  1  0     0  0
  6 10  1  0     0  0
  7 17  1  0     0  0
  8  3  1  0     0  0
  9  1  2  0     0  0
 10 13  1  0     0  0
 11  3  1  0     0  0
 12  2  1  0     0  0
 13  7  1  0     0  0
 14  4  1  0     0  0
 15 14  1  0     0  0
 16 14  1  0     0  0
 17 16  1  0     0  0
 18 15  1  0     0  0
 19 11  1  0     0  0
 20  5  1  0     0  0
 21  5  1  0     0  0
  2 22  1  6     0  0
 23  6  2  0     0  0
 24 25  1  0     0  0
 25 10  2  0     0  0
 26 12  1  0     0  0
 27 12  2  0     0  0
 28 24  2  0     0  0
 29 23  1  0     0  0
 30 26  2  0     0  0
 31 27  1  0     0  0
 32 31  2  0     0  0
  3 33  1  6     0  0
  7 18  1  0     0  0
 19  5  1  0     0  0
 32 30  1  0     0  0
 28 23  1  0     0  0
M  END""")
    nm = Draw.PrepareMolForDrawing(m)
    self.assertEqual(nm.GetBondBetweenAtoms(2, 1).GetBondType(), Chem.BondType.SINGLE)
    self.assertEqual(nm.GetBondBetweenAtoms(2, 1).GetBondDir(), Chem.BondDir.NONE)
    self.assertEqual(nm.GetBondBetweenAtoms(2, 7).GetBondType(), Chem.BondType.SINGLE)
    self.assertEqual(nm.GetBondBetweenAtoms(2, 7).GetBondDir(), Chem.BondDir.BEGINWEDGE)
    nm = Draw.PrepareMolForDrawing(nm)
    self.assertEqual(nm.GetBondBetweenAtoms(2, 1).GetBondType(), Chem.BondType.SINGLE)
    self.assertEqual(nm.GetBondBetweenAtoms(2, 1).GetBondDir(), Chem.BondDir.NONE)
    self.assertEqual(nm.GetBondBetweenAtoms(2, 7).GetBondType(), Chem.BondType.SINGLE)
    self.assertEqual(nm.GetBondBetweenAtoms(2, 7).GetBondDir(), Chem.BondDir.BEGINWEDGE)


if __name__ == "__main__":
  unittest.main()
