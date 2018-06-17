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

    invGen = rdAtomPairGenerator.getAtomPairAtomInvGen()
    g = rdAtomPairGenerator.getAtomPairGeneratorWrapped(atomInvGeneratorWrapper = invGen)
    fp = g.getFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 2)

    g = rdAtomPairGenerator.getAtomPairGeneratorWrapped(minDistance = 2)
    fp = g.getFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 1)

    g = rdAtomPairGenerator.getAtomPairGeneratorWrapped(maxDistance = 1)
    fp = g.getFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 1)

    g = rdAtomPairGenerator.getAtomPairGeneratorWrapped(useCountSimulation = False)
    fp = g.getFingerprintAsBitVect(m)
    nzc = fp.GetNumOnBits()
    self.assertEqual(nzc, 2)

    invGen = rdAtomPairGenerator.getAtomPairAtomInvGen(includeChirality = False)
    invGenChirality = rdAtomPairGenerator.getAtomPairAtomInvGen(includeChirality = True)
    g = rdAtomPairGenerator.getAtomPairGeneratorWrapped(includeChirality = False, atomInvGeneratorWrapper = invGen)
    gChirality = rdAtomPairGenerator.getAtomPairGeneratorWrapped(includeChirality = True, atomInvGeneratorWrapper = invGenChirality)
    fp = g.getFingerprint(m)
    nz = fp.GetNonzeroElements()
    fpChirality = gChirality.getFingerprint(m)
    nzChirality = fpChirality.GetNonzeroElements()
    self.assertNotEqual(nz.keys()[0], nzChirality.keys()[0])


if __name__ == '__main__':
  unittest.main()
