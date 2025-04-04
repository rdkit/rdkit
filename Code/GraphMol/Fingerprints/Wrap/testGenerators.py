import unittest

import numpy as np

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


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

    g = rdFingerprintGenerator.GetAtomPairGenerator(
      atomInvariantsGenerator=rdFingerprintGenerator.GetAtomPairAtomInvGen())
    fp = g.GetSparseCountFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 2)

    g = rdFingerprintGenerator.GetAtomPairGenerator(minDistance=2)
    fp = g.GetSparseCountFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 1)

    g = rdFingerprintGenerator.GetAtomPairGenerator(maxDistance=1)
    fp = g.GetSparseCountFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 1)

    g = rdFingerprintGenerator.GetAtomPairGenerator(countSimulation=False)
    fp = g.GetSparseFingerprint(m)
    nzc = fp.GetNumOnBits()
    self.assertEqual(nzc, 2)

    invGen = rdFingerprintGenerator.GetAtomPairAtomInvGen(includeChirality=False)
    invGenChirality = rdFingerprintGenerator.GetAtomPairAtomInvGen(includeChirality=True)
    g = rdFingerprintGenerator.GetAtomPairGenerator(includeChirality=False,
                                                    atomInvariantsGenerator=invGen)
    gChirality = rdFingerprintGenerator.GetAtomPairGenerator(
      includeChirality=True, atomInvariantsGenerator=invGenChirality)
    fp = g.GetSparseCountFingerprint(m)
    nz = fp.GetNonzeroElements()
    fpChirality = gChirality.GetSparseCountFingerprint(m)
    nzChirality = fpChirality.GetNonzeroElements()
    self.assertNotEqual(nz.keys(), nzChirality.keys())

  def testMorganGenerator(self):
    m = Chem.MolFromSmiles('CCCC(=O)O')
    g = rdFingerprintGenerator.GetMorganGenerator(3)
    fp = g.GetSparseCountFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 14)

    invgen = rdFingerprintGenerator.GetMorganAtomInvGen()
    g = rdFingerprintGenerator.GetMorganGenerator(radius=3, atomInvariantsGenerator=invgen)
    fp = g.GetSparseCountFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 14)

    invgen = rdFingerprintGenerator.GetMorganFeatureAtomInvGen()
    g = rdFingerprintGenerator.GetMorganGenerator(radius=3, atomInvariantsGenerator=invgen)
    fp = g.GetSparseCountFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 13)

    ms = [Chem.MolFromSmiles(x, sanitize=False) for x in ('C1=CC=CN=N1', 'C1C=CC=NN=1')]
    for m in ms:
      m.UpdatePropertyCache()
      Chem.GetSymmSSSR(m)

    g = rdFingerprintGenerator.GetMorganGenerator(radius=2, useBondTypes=True)
    self.assertNotEqual(g.GetSparseCountFingerprint(ms[0]), g.GetSparseCountFingerprint(ms[1]))
    g = rdFingerprintGenerator.GetMorganGenerator(radius=2, useBondTypes=False)
    self.assertEqual(g.GetSparseCountFingerprint(ms[0]), g.GetSparseCountFingerprint(ms[1]))

    binvgen = rdFingerprintGenerator.GetMorganBondInvGen(useBondTypes=False)
    g2 = rdFingerprintGenerator.GetMorganGenerator(radius=2, bondInvariantsGenerator=binvgen)
    self.assertEqual(g.GetSparseCountFingerprint(ms[0]), g2.GetSparseCountFingerprint(ms[0]))
    self.assertEqual(g.GetSparseCountFingerprint(ms[1]), g2.GetSparseCountFingerprint(ms[1]))

  def testRDKitFPGenerator(self):
    m = Chem.MolFromSmiles('CCCCC')
    g = rdFingerprintGenerator.GetRDKitFPGenerator()
    fp = g.GetSparseCountFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 8)

  def testTopologicalTorsionGenerator(self):
    m = Chem.MolFromSmiles('CCCCC')
    g = rdFingerprintGenerator.GetTopologicalTorsionGenerator()
    fp = g.GetSparseCountFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 1)

  def testBulk(self):
    m1 = Chem.MolFromSmiles('CCC')
    m2 = Chem.MolFromSmiles('OCCCCC')
    m3 = Chem.MolFromSmiles('CCCCC')

    g = rdFingerprintGenerator.GetAtomPairGenerator()
    results = rdFingerprintGenerator.GetSparseCountFPs([m1, m2, m3],
                                                       rdFingerprintGenerator.AtomPairFP)
    self.assertEqual(results[0], g.GetSparseCountFingerprint(m1))
    self.assertEqual(results[1], g.GetSparseCountFingerprint(m2))
    self.assertEqual(results[2], g.GetSparseCountFingerprint(m3))
    self.assertEqual(len(results), 3)

    g = rdFingerprintGenerator.GetMorganGenerator(2)
    results = rdFingerprintGenerator.GetSparseCountFPs([m1, m2, m3],
                                                       rdFingerprintGenerator.MorganFP)
    self.assertEqual(results[0], g.GetSparseCountFingerprint(m1))
    self.assertEqual(results[1], g.GetSparseCountFingerprint(m2))
    self.assertEqual(results[2], g.GetSparseCountFingerprint(m3))
    self.assertEqual(len(results), 3)

    g = rdFingerprintGenerator.GetRDKitFPGenerator()
    results = rdFingerprintGenerator.GetSparseCountFPs([m1, m2, m3], rdFingerprintGenerator.RDKitFP)
    self.assertEqual(results[0], g.GetSparseCountFingerprint(m1))
    self.assertEqual(results[1], g.GetSparseCountFingerprint(m2))
    self.assertEqual(results[2], g.GetSparseCountFingerprint(m3))
    self.assertEqual(len(results), 3)

    g = rdFingerprintGenerator.GetTopologicalTorsionGenerator()
    results = rdFingerprintGenerator.GetSparseCountFPs([m1, m2, m3],
                                                       rdFingerprintGenerator.TopologicalTorsionFP)
    self.assertEqual(results[0], g.GetSparseCountFingerprint(m1))
    self.assertEqual(results[1], g.GetSparseCountFingerprint(m2))
    self.assertEqual(results[2], g.GetSparseCountFingerprint(m3))
    self.assertEqual(len(results), 3)

    g = rdFingerprintGenerator.GetAtomPairGenerator()
    results = rdFingerprintGenerator.GetSparseFPs([m1, m2, m3], rdFingerprintGenerator.AtomPairFP)
    self.assertEqual(results[0], g.GetSparseFingerprint(m1))
    self.assertEqual(results[1], g.GetSparseFingerprint(m2))
    self.assertEqual(results[2], g.GetSparseFingerprint(m3))
    self.assertEqual(len(results), 3)

    g = rdFingerprintGenerator.GetMorganGenerator(2)
    results = rdFingerprintGenerator.GetSparseFPs([m1, m2, m3], rdFingerprintGenerator.MorganFP)
    self.assertEqual(results[0], g.GetSparseFingerprint(m1))
    self.assertEqual(results[1], g.GetSparseFingerprint(m2))
    self.assertEqual(results[2], g.GetSparseFingerprint(m3))
    self.assertEqual(len(results), 3)

    g = rdFingerprintGenerator.GetRDKitFPGenerator()
    results = rdFingerprintGenerator.GetSparseFPs([m1, m2, m3], rdFingerprintGenerator.RDKitFP)
    self.assertEqual(results[0], g.GetSparseFingerprint(m1))
    self.assertEqual(results[1], g.GetSparseFingerprint(m2))
    self.assertEqual(results[2], g.GetSparseFingerprint(m3))
    self.assertEqual(len(results), 3)

    g = rdFingerprintGenerator.GetTopologicalTorsionGenerator()
    results = rdFingerprintGenerator.GetSparseFPs([m1, m2, m3],
                                                  rdFingerprintGenerator.TopologicalTorsionFP)
    self.assertEqual(results[0], g.GetSparseFingerprint(m1))
    self.assertEqual(results[1], g.GetSparseFingerprint(m2))
    self.assertEqual(results[2], g.GetSparseFingerprint(m3))
    self.assertEqual(len(results), 3)

    g = rdFingerprintGenerator.GetAtomPairGenerator()
    results = rdFingerprintGenerator.GetCountFPs([m1, m2, m3], rdFingerprintGenerator.AtomPairFP)
    self.assertEqual(results[0], g.GetCountFingerprint(m1))
    self.assertEqual(results[1], g.GetCountFingerprint(m2))
    self.assertEqual(results[2], g.GetCountFingerprint(m3))
    self.assertEqual(len(results), 3)

    g = rdFingerprintGenerator.GetMorganGenerator(2)
    results = rdFingerprintGenerator.GetCountFPs([m1, m2, m3], rdFingerprintGenerator.MorganFP)
    self.assertEqual(results[0], g.GetCountFingerprint(m1))
    self.assertEqual(results[1], g.GetCountFingerprint(m2))
    self.assertEqual(results[2], g.GetCountFingerprint(m3))
    self.assertEqual(len(results), 3)

    g = rdFingerprintGenerator.GetRDKitFPGenerator()
    results = rdFingerprintGenerator.GetCountFPs([m1, m2, m3], rdFingerprintGenerator.RDKitFP)
    self.assertEqual(results[0], g.GetCountFingerprint(m1))
    self.assertEqual(results[1], g.GetCountFingerprint(m2))
    self.assertEqual(results[2], g.GetCountFingerprint(m3))
    self.assertEqual(len(results), 3)

    g = rdFingerprintGenerator.GetTopologicalTorsionGenerator()
    results = rdFingerprintGenerator.GetCountFPs([m1, m2, m3],
                                                 rdFingerprintGenerator.TopologicalTorsionFP)
    self.assertEqual(results[0], g.GetCountFingerprint(m1))
    self.assertEqual(results[1], g.GetCountFingerprint(m2))
    self.assertEqual(results[2], g.GetCountFingerprint(m3))
    self.assertEqual(len(results), 3)

  def testNumBitsPerFeature(self):
    m1 = Chem.MolFromSmiles('CCCO')
    g = rdFingerprintGenerator.GetRDKitFPGenerator(minPath=1, maxPath=2)
    fp = g.GetFingerprint(m1)
    self.assertEqual(fp.GetNumOnBits(), 8)

    g = rdFingerprintGenerator.GetRDKitFPGenerator(minPath=1, maxPath=2, numBitsPerFeature=1)
    fp = g.GetFingerprint(m1)
    self.assertEqual(fp.GetNumOnBits(), 4)

  def testAdditionalOutput(self):
    m1 = Chem.MolFromSmiles('CCO')
    g = rdFingerprintGenerator.GetAtomPairGenerator()
    ao = rdFingerprintGenerator.AdditionalOutput()
    ao.AllocateAtomCounts()
    fp = g.GetFingerprint(m1, additionalOutput=ao)
    self.assertEqual(ao.GetAtomCounts(), (2, 2, 2))
    self.assertIsNone(ao.GetAtomToBits())
    self.assertIsNone(ao.GetBitInfoMap())
    self.assertIsNone(ao.GetBitPaths())

    ao = rdFingerprintGenerator.AdditionalOutput()
    ao.AllocateAtomToBits()
    fp = g.GetFingerprint(m1, additionalOutput=ao)
    self.assertIsNone(ao.GetAtomCounts())
    self.assertEqual(ao.GetAtomToBits(), ((1404, 1916), (1404, 1596), (1596, 1916)))
    self.assertIsNone(ao.GetBitInfoMap())
    self.assertIsNone(ao.GetBitPaths())

    ao = rdFingerprintGenerator.AdditionalOutput()
    ao.AllocateBitInfoMap()
    fp = g.GetFingerprint(m1, additionalOutput=ao)
    self.assertIsNone(ao.GetAtomCounts())
    self.assertIsNone(ao.GetAtomToBits())
    self.assertEqual(ao.GetBitInfoMap(), {1404: ((0, 1), ), 1596: ((1, 2), ), 1916: ((0, 2), )})
    self.assertIsNone(ao.GetBitPaths())

  def testCountBounds(self):
    m = Chem.MolFromSmiles('COc1ccc(CCNC(=O)c2ccccc2C(=O)NCCc2ccc(OC)cc2)cc1')
    fp1 = rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=2048,
                                                     countSimulation=True).GetFingerprint(m)
    fp2 = rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=2048, countSimulation=True,
                                                     countBounds=(1, 8, 16, 32)).GetFingerprint(m)
    self.assertNotEqual(fp1.GetNumOnBits(), fp2.GetNumOnBits())
    fp1 = rdFingerprintGenerator.GetTopologicalTorsionGenerator(
      fpSize=2048, countSimulation=True).GetFingerprint(m)
    fp2 = rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=2048, countSimulation=True,
                                                                countBounds=(1, 8, 16,
                                                                             32)).GetFingerprint(m)
    self.assertNotEqual(fp1.GetNumOnBits(), fp2.GetNumOnBits())
    fp1 = rdFingerprintGenerator.GetMorganGenerator(fpSize=2048,
                                                    countSimulation=True).GetFingerprint(m)
    fp2 = rdFingerprintGenerator.GetMorganGenerator(fpSize=2048, countSimulation=True,
                                                    countBounds=(1, 8, 16, 32)).GetFingerprint(m)
    self.assertNotEqual(fp1.GetNumOnBits(), fp2.GetNumOnBits())
    fp1 = rdFingerprintGenerator.GetAtomPairGenerator(fpSize=2048,
                                                      countSimulation=True).GetFingerprint(m)
    fp2 = rdFingerprintGenerator.GetAtomPairGenerator(fpSize=2048, countSimulation=True,
                                                      countBounds=(1, 8, 16, 32)).GetFingerprint(m)
    self.assertNotEqual(fp1.GetNumOnBits(), fp2.GetNumOnBits())

  def testNumpyFingerprints(self):
    m = Chem.MolFromSmiles('COc1ccc(CCNC(=O)c2ccccc2C(=O)NCCc2ccc(OC)cc2)cc1')
    for fn in (rdFingerprintGenerator.GetRDKitFPGenerator,
               rdFingerprintGenerator.GetMorganGenerator,
               rdFingerprintGenerator.GetAtomPairGenerator,
               rdFingerprintGenerator.GetTopologicalTorsionGenerator):
      gen = fn(fpSize=2048)
      bv = gen.GetFingerprint(m)
      oarr = np.zeros((bv.GetNumBits(), ), 'u1')
      DataStructs.ConvertToNumpyArray(bv, oarr)
      arr = gen.GetFingerprintAsNumPy(m)
      np.testing.assert_array_equal(oarr, arr)

      fp = gen.GetCountFingerprint(m)
      oarr = np.zeros((fp.GetLength(), ), 'u4')
      DataStructs.ConvertToNumpyArray(fp, oarr)
      arr = gen.GetCountFingerprintAsNumPy(m)
      np.testing.assert_array_equal(oarr, arr)

  def testMorganRedundantEnvironments(self):
    m = Chem.MolFromSmiles('CC(=O)O')

    g = rdFingerprintGenerator.GetMorganGenerator(2)
    fp = g.GetSparseCountFingerprint(m)
    self.assertEqual(fp.GetTotalVal(), 8)

    g = rdFingerprintGenerator.GetMorganGenerator(2, includeRedundantEnvironments=True)
    fp = g.GetSparseCountFingerprint(m)
    self.assertEqual(fp.GetTotalVal(), 12)

  def testFingerprintOptions(self):
    m = Chem.MolFromSmiles('CC(=O)O')

    g = rdFingerprintGenerator.GetMorganGenerator(2)
    fp = g.GetSparseCountFingerprint(m)
    self.assertEqual(fp.GetTotalVal(), 8)

    g.GetOptions().numBitsPerFeature = 2
    fp = g.GetSparseCountFingerprint(m)
    self.assertEqual(fp.GetTotalVal(), 16)

    g.GetOptions().includeRedundantEnvironments = True
    fp = g.GetSparseCountFingerprint(m)
    self.assertEqual(fp.GetTotalVal(), 24)

    m = Chem.MolFromSmiles('CC(=O)C')
    g = rdFingerprintGenerator.GetAtomPairGenerator()
    fp = g.GetFingerprint(m)
    self.assertEqual(fp.GetNumOnBits(), 6)
    g.GetOptions().countSimulation = False
    fp = g.GetFingerprint(m)
    self.assertEqual(fp.GetNumOnBits(), 4)
    g.GetOptions().maxDistance = 1
    fp = g.GetFingerprint(m)
    self.assertEqual(fp.GetNumOnBits(), 2)

    m = Chem.MolFromSmiles('OC(C)(C)C')
    g = rdFingerprintGenerator.GetAtomPairGenerator()
    fp = g.GetFingerprint(m)
    self.assertEqual(fp.GetNumOnBits(), 7)
    g.GetOptions().SetCountBounds((1, 2, 3, 4))
    fp = g.GetFingerprint(m)
    self.assertEqual(fp.GetNumOnBits(), 10)

  def testTopologicalTorsionShortestPaths(self):
    m = Chem.MolFromSmiles('CC1CCC1')
    g = rdFingerprintGenerator.GetTopologicalTorsionGenerator()
    fp = g.GetSparseCountFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 3)

    g.GetOptions().onlyShortestPaths = True
    fp = g.GetSparseCountFingerprint(m)
    nz = fp.GetNonzeroElements()
    self.assertEqual(len(nz), 1)

  def testMorganGeneratorMultiMol(self):
    smis = [
      'CC1CCC1',
      'CCC1CCC1',
      'CCCC1CCC1',
      'CC1CC(O)C1',
      'CC1CC(OC)C1',
    ]
    for i in range(4):
      smis = smis + smis
    ms = [Chem.MolFromSmiles(smi) for smi in smis]
    g = rdFingerprintGenerator.GetMorganGenerator()
    ofps = tuple([g.GetFingerprint(m) for m in ms])
    tfps = g.GetFingerprints(ms, numThreads=1)
    for ofp, tfp in zip(ofps, tfps):
      self.assertEqual(ofp, tfp)
    tfps = g.GetFingerprints(ms, numThreads=4)
    for ofp, tfp in zip(ofps, tfps):
      self.assertEqual(ofp, tfp)

    ofps = tuple([g.GetCountFingerprint(m) for m in ms])
    tfps = g.GetCountFingerprints(ms, numThreads=1)
    for ofp, tfp in zip(ofps, tfps):
      self.assertEqual(ofp, tfp)
    tfps = g.GetCountFingerprints(ms, numThreads=4)
    for ofp, tfp in zip(ofps, tfps):
      self.assertEqual(ofp, tfp)

    ofps = tuple([g.GetSparseFingerprint(m) for m in ms])
    tfps = g.GetSparseFingerprints(ms, numThreads=1)
    for ofp, tfp in zip(ofps, tfps):
      self.assertEqual(ofp, tfp)
    tfps = g.GetSparseFingerprints(ms, numThreads=4)
    for ofp, tfp in zip(ofps, tfps):
      self.assertEqual(ofp, tfp)

    ofps = tuple([g.GetSparseCountFingerprint(m) for m in ms])
    tfps = g.GetSparseCountFingerprints(ms, numThreads=1)
    for ofp, tfp in zip(ofps, tfps):
      self.assertEqual(ofp, tfp)
    tfps = g.GetSparseCountFingerprints(ms, numThreads=4)
    for ofp, tfp in zip(ofps, tfps):
      self.assertEqual(ofp, tfp)

  def testFingerprintGeneratorOptionsLifetime(self):
    # this should not result in a seg fault
    import inspect
    inspect.getmembers(rdFingerprintGenerator.GetRDKitFPGenerator().GetOptions())


if __name__ == '__main__':
  unittest.main()
