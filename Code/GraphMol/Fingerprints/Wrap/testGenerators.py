from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
import unittest

class TestCase(unittest.TestCase):

  def setUp(self):
    pass


  def testAtomPairGenerator(self):
    m = Chem.MolFromSmiles('CCC')
    g = rdFingerprintGenerator.GetAtomPairGenerator()
    fp = g.GetFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 2)

    fp = g.GetFoldedFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 2)

    fp = g.GetFingerprintAsBitVect(m)
    nzc = fp.GetNumOnBits()
    self.assertEqual(nzc, 3)

    fp = g.GetFoldedFingerprintAsBitVect(m)
    nzc = fp.GetNumOnBits()
    self.assertEqual(nzc, 3)

    g = rdFingerprintGenerator.GetAtomPairGenerator(atomInvariantsGenerator = rdFingerprintGenerator.GetAtomPairAtomInvGen() )
    fp = g.GetFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 2)

    g = rdFingerprintGenerator.GetAtomPairGenerator(minDistance = 2)
    fp = g.GetFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 1)

    g = rdFingerprintGenerator.GetAtomPairGenerator(maxDistance = 1)
    fp = g.GetFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 1)

    g = rdFingerprintGenerator.GetAtomPairGenerator(useCountSimulation = False)
    fp = g.GetFingerprintAsBitVect(m)
    nzc = fp.GetNumOnBits()
    self.assertEqual(nzc, 2)

    invGen = rdFingerprintGenerator.GetAtomPairAtomInvGen(includeChirality = False)
    invGenChirality = rdFingerprintGenerator.GetAtomPairAtomInvGen(includeChirality = True)
    g = rdFingerprintGenerator.GetAtomPairGenerator(includeChirality = False, atomInvariantsGenerator = invGen)
    gChirality = rdFingerprintGenerator.GetAtomPairGenerator(includeChirality = True, atomInvariantsGenerator = invGenChirality)
    fp = g.GetFingerprint(m)
    nz = fp.GetNonzeroElements()
    fpChirality = gChirality.GetFingerprint(m)
    nzChirality = fpChirality.GetNonzeroElements()
    self.assertNotEqual(nz.keys()[0], nzChirality.keys()[0])

  def testMorganGenerator(self):
    m = Chem.MolFromSmiles('CCCCC')
    g = rdFingerprintGenerator.GetMorganGenerator32(3)
    fp = g.GetFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 7)

    m = Chem.MolFromSmiles('CCCCC')
    g = rdFingerprintGenerator.GetMorganGenerator64(3)
    fp = g.GetFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 7)

  def testRDKitFPGenerator(self):
    m = Chem.MolFromSmiles('CCCCC')
    g = rdFingerprintGenerator.GetRDKitFPGenerator32()
    fp = g.GetFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 4)

    m = Chem.MolFromSmiles('CCCCC')
    g = rdFingerprintGenerator.GetRDKitFPGenerator64()
    fp = g.GetFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 4)
  
  def testTopologicalTorsionGenerator(self):

    m = Chem.MolFromSmiles('CCCCC')
    g = rdFingerprintGenerator.GetTopologicalTorsionGenerator64()
    fp = g.GetFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 1)


if __name__ == '__main__':
  unittest.main()
