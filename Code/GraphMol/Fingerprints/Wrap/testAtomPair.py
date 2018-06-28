from rdkit import Chem
from rdkit.Chem import rdAtomPairGenerator
import unittest

class TestCase(unittest.TestCase):

  def setUp(self):
    pass


  def testAtomPairGenerator(self):
    m = Chem.MolFromSmiles('CCC')
    g = rdAtomPairGenerator.GetAtomPairGenerator()
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

    g = rdAtomPairGenerator.GetAtomPairGenerator(atomInvariantsGenerator = rdAtomPairGenerator.GetAtomPairAtomInvGen() )
    fp = g.GetFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 2)

    g = rdAtomPairGenerator.GetAtomPairGenerator(minDistance = 2)
    fp = g.GetFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 1)

    g = rdAtomPairGenerator.GetAtomPairGenerator(maxDistance = 1)
    fp = g.GetFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 1)

    g = rdAtomPairGenerator.GetAtomPairGenerator(useCountSimulation = False)
    fp = g.GetFingerprintAsBitVect(m)
    nzc = fp.GetNumOnBits()
    self.assertEqual(nzc, 2)

    invGen = rdAtomPairGenerator.GetAtomPairAtomInvGen(includeChirality = False)
    invGenChirality = rdAtomPairGenerator.GetAtomPairAtomInvGen(includeChirality = True)
    g = rdAtomPairGenerator.GetAtomPairGenerator(includeChirality = False, atomInvariantsGenerator = invGen)
    gChirality = rdAtomPairGenerator.GetAtomPairGenerator(includeChirality = True, atomInvariantsGenerator = invGenChirality)
    fp = g.GetFingerprint(m)
    nz = fp.GetNonzeroElements()
    fpChirality = gChirality.GetFingerprint(m)
    nzChirality = fpChirality.GetNonzeroElements()
    self.assertNotEqual(nz.keys()[0], nzChirality.keys()[0])


if __name__ == '__main__':
  unittest.main()
