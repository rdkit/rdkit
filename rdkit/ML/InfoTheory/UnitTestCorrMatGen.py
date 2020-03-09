from rdkit.ML.InfoTheory import rdInfoTheory
from rdkit.ML.Data import DataUtils
from rdkit import DataStructs
import unittest
import random


try:
    from rdkit.ML.InfoTheory import BitClusterer
except ImportError:
    BitClusterer = None


def getValLTM(i, j, mat):
    if i > j:
        id_ = i * (i - 1) // 2 + j
        return mat[id_]
    elif j > i:
        id_ = j * (j - 1) // 2 + i
        return mat[id_]
    else:
        return 0.0


class TestCase(unittest.TestCase):

    def setUp(self):
        # here is what we are going to do to test this out
        # - generate bit vectors of length nbits
        # - turn on a fraction of the first nbits/2 bits at random
        # - for each bit i turned on in the range (0, nbits/2) turn on the bit
        #   nbits/2 + i
        # - basically the first half of a fingerprint is same as the second half of the
        #   fingerprint
        # - if we repeat this process often enough we whould see strong correlation between
        #   the bits i (i < nbits/2) and (nbits/2 + i)
        DataUtils.InitRandomNumbers((100, 23))
        self.nbits = 200
        self.d = 40
        self.nfp = 1000

        self.blist = list(range(self.nbits))

        self.fps = []
        for _ in range(self.nfp):
            fp = DataStructs.ExplicitBitVect(self.nbits)
            obits = list(range(self.nbits // 2))
            random.shuffle(obits, random=random.random)
            for bit in obits[0:self.d]:
                fp.SetBit(bit)
                fp.SetBit(bit + self.nbits // 2)
            self.fps.append(fp)

    def test_getValLTM(self):
        #   - 1 2 4
        #   1 - 3 5
        #   2 3 - 6
        #   4 5 6 -
        mat = list(range(1, 7, 1))
        for i in range(4):
            self.assertEqual(getValLTM(i, i, mat), 0.0)
        self.assertEqual(getValLTM(0, 1, mat), 1)
        self.assertEqual(getValLTM(0, 2, mat), 2)
        self.assertEqual(getValLTM(0, 3, mat), 4)
        self.assertEqual(getValLTM(1, 0, mat), 1)
        self.assertEqual(getValLTM(2, 0, mat), 2)
        self.assertEqual(getValLTM(3, 0, mat), 4)
        self.assertEqual(getValLTM(1, 2, mat), 3)
        self.assertEqual(getValLTM(1, 3, mat), 5)
        self.assertEqual(getValLTM(2, 1, mat), 3)
        self.assertEqual(getValLTM(3, 1, mat), 5)
        self.assertEqual(getValLTM(2, 3, mat), 6)
        self.assertEqual(getValLTM(3, 2, mat), 6)

    def test0CorrMat(self):
        cmg = rdInfoTheory.BitCorrMatGenerator()
        cmg.SetBitList(self.blist)
        for fp in self.fps:
            cmg.CollectVotes(fp)

        corrMat = cmg.GetCorrMatrix()

        avr = 0.0
        navr = 0.0
        for i in range(self.nbits // 2):
            avr += getValLTM(i, i + self.nbits // 2, corrMat)
            navr += getValLTM(i, i + 1, corrMat)

        self.assertEqual(2 * avr / self.nbits, 400.0)
        self.assertEqual(2 * navr / self.nbits, 158.3)

    @unittest.skipIf(BitClusterer is None, 'Cannot import BitClusterer')
    def test1Cluster(self):
        cmg = rdInfoTheory.BitCorrMatGenerator()
        cmg.SetBitList(self.blist)
        for fp in self.fps:
            cmg.CollectVotes(fp)

        corrMat = cmg.GetCorrMatrix()

        bcl = BitClusterer.BitClusterer(self.blist, self.nbits // 2)
        bcl.ClusterBits(corrMat)
        cls = bcl.GetClusters()
        for cl in cls:
            self.assertEqual(len(cl), 2)
            self.assertEqual((cl[0] + self.nbits // 2), cl[1])
        bcl.SetClusters(cls)
        self.assertRaises(AssertionError, bcl.SetClusters, cls[:-1])

        tfp = DataStructs.ExplicitBitVect(self.nbits)
        obits = list(range(0, self.nbits // 4)) + list(range(self.nbits // 2, 3 * self.nbits // 4))
        tfp.SetBitsFromList(obits)
        rvc = bcl.MapToClusterScores(tfp)
        self.assertEqual(len(rvc), self.nbits // 2)
        for i in range(self.nbits // 2):
            if i < self.nbits // 4:
                self.assertEqual(rvc[i], 2)
            else:
                self.assertEqual(rvc[i], 0)

        nfp = bcl.MapToClusterFP(tfp)
        self.assertEqual(len(nfp), self.nbits // 2)
        for i in range(self.nbits // 2):
            if i < self.nbits // 4:
                self.assertTrue(nfp[i])
            else:
                self.assertFalse(nfp[i])


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
