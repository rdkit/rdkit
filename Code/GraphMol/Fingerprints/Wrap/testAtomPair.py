from rdkit import Chem
from rdkit.Chem import rdAtomPairGenerator
import unittest

class TestCase(unittest.TestCase):

  def setUp(self):
    pass


  def testAtomPairGenerator(self):
    m = Chem.MolFromSmiles('CCC')
    g = rdAtomPairGenerator.getAtomPairGeneratorWrapped()
    fp = g.getFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 2)

    fp = g.getFoldedFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 2)

    fp = g.getFingerprintAsBitVect(m)
    nzc = fp.GetNumOnBits()
    self.assertEqual(nzc, 3)

    fp = g.getFoldedFingerprintAsBitVect(m)
    nzc = fp.GetNumOnBits()
    self.assertEqual(nzc, 3)

if __name__ == '__main__':
  unittest.main()
