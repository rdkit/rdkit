from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
import unittest

class TestCase(unittest.TestCase):

  def setUp(self):
    pass


  def testAtomPairGenerator(self):
    m = Chem.MolFromSmiles('CCC')
    g = rdFingerprintGenerator.GetAtomPairGenerator()
    fp = g.GetSparseCountFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 2)

    fp = g.GetCountFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 2)

    fp = g.GetSparseFingerprint(m)
    nzc = fp.GetNumOnBits()
    self.assertEqual(nzc, 3)

    fp = g.GetFingerprint(m)
    nzc = fp.GetNumOnBits()
    self.assertEqual(nzc, 3)

    g = rdFingerprintGenerator.GetAtomPairGenerator(atomInvariantsGenerator = rdFingerprintGenerator.GetAtomPairAtomInvGen() )
    fp = g.GetSparseCountFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 2)

    g = rdFingerprintGenerator.GetAtomPairGenerator(minDistance = 2)
    fp = g.GetSparseCountFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 1)

    g = rdFingerprintGenerator.GetAtomPairGenerator(maxDistance = 1)
    fp = g.GetSparseCountFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 1)

    g = rdFingerprintGenerator.GetAtomPairGenerator(useCountSimulation = False)
    fp = g.GetSparseFingerprint(m)
    nzc = fp.GetNumOnBits()
    self.assertEqual(nzc, 2)

    invGen = rdFingerprintGenerator.GetAtomPairAtomInvGen(includeChirality = False)
    invGenChirality = rdFingerprintGenerator.GetAtomPairAtomInvGen(includeChirality = True)
    g = rdFingerprintGenerator.GetAtomPairGenerator(includeChirality = False, atomInvariantsGenerator = invGen)
    gChirality = rdFingerprintGenerator.GetAtomPairGenerator(includeChirality = True, atomInvariantsGenerator = invGenChirality)
    fp = g.GetSparseCountFingerprint(m)
    nz = fp.GetNonzeroElements()
    fpChirality = gChirality.GetSparseCountFingerprint(m)
    nzChirality = fpChirality.GetNonzeroElements()
    self.assertNotEqual(nz.keys()[0], nzChirality.keys()[0])

  def testMorganGenerator(self):
    m = Chem.MolFromSmiles('CCCCC')
    g = rdFingerprintGenerator.GetMorganGenerator(3)
    fp = g.GetSparseCountFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 7)

  def testRDKitFPGenerator(self):
    m = Chem.MolFromSmiles('CCCCC')
    g = rdFingerprintGenerator.GetRDKitFPGenerator()
    fp = g.GetSparseCountFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 4)
 
  def testTopologicalTorsionGenerator(self):
    m = Chem.MolFromSmiles('CCCCC')
    g = rdFingerprintGenerator.GetTopologicalTorsionGenerator()
    fp = g.GetSparseCountFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 1)


if __name__ == '__main__':
  unittest.main()
