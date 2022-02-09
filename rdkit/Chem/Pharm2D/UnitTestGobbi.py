#
#  Copyright (C) 2002-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""unit testing code for the signatures

"""


import os
import unittest

from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate


class TestCase(unittest.TestCase):

    def setUp(self):
        self.factory = Gobbi_Pharm2D.factory

    def test1Sigs(self):
        probes = [
          ('OCCC=O', {
            'HA': (1, ((0, ), (4, ))),
            'HD': (1, ((0, ), )),
            'LH': (0, None),
            'AR': (0, None),
            'RR': (0, None),
            'X': (0, None),
            'BG': (0, None),
            'AG': (0, None),
          }),
          ('OCCC(=O)O', {
            'HA': (1, ((0, ), (4, ))),
            'HD': (1, ((0, ), (5, ))),
            'LH': (0, None),
            'AR': (0, None),
            'RR': (0, None),
            'X': (0, None),
            'BG': (0, None),
            'AG': (1, ((3, ), )),
          }),
          ('CCCN', {
            'HA': (1, ((3, ), )),
            'HD': (1, ((3, ), )),
            'LH': (0, None),
            'AR': (0, None),
            'RR': (0, None),
            'X': (0, None),
            'BG': (1, ((3, ), )),
            'AG': (0, None),
          }),
          ('CCCCC', {
            'HA': (0, None),
            'HD': (0, None),
            'LH': (1, ((1, ), (3, ))),
            'AR': (0, None),
            'RR': (0, None),
            'X': (0, None),
            'BG': (0, None),
            'AG': (0, None),
          }),
          ('CC1CCC1', {
            'HA': (0, None),
            'HD': (0, None),
            'LH': (1, ((1, ), (3, ))),
            'AR': (0, None),
            'RR': (1, ((1, ), )),
            'X': (0, None),
            'BG': (0, None),
            'AG': (0, None),
          }),
          ('[SiH3]C1CCC1', {
            'HA': (0, None),
            'HD': (0, None),
            'LH': (1, ((1, ), )),
            'AR': (0, None),
            'RR': (1, ((1, ), )),
            'X': (1, ((0, ), )),
            'BG': (0, None),
            'AG': (0, None),
          }),
          ('[SiH3]c1ccccc1', {
            'HA': (0, None),
            'HD': (0, None),
            'LH': (0, None),
            'AR': (1, ((1, ), )),
            'RR': (0, None),
            'X': (1, ((0, ), )),
            'BG': (0, None),
            'AG': (0, None),
          }),
        ]
        for smi, d in probes:
            mol = Chem.MolFromSmiles(smi)
            feats = self.factory.featFactory.GetFeaturesForMol(mol)
            for k in d.keys():
                shouldMatch, mapList = d[k]
                feats = self.factory.featFactory.GetFeaturesForMol(mol, includeOnly=k)
                if shouldMatch:
                    self.assertTrue(feats)
                    self.assertEqual(len(feats), len(mapList))
                    aids = [(x.GetAtomIds()[0], ) for x in feats]
                    aids.sort()
                    self.assertEqual(tuple(aids), mapList)

    def test2Sigs(self):
        probes = [('O=CCC=O', (149, )),
                  ('OCCC=O', (149, 156)),
                  ('OCCC(=O)O', (22, 29, 149, 154, 156, 184, 28822, 30134)), ]
        for smi, tgt in probes:
            sig = Generate.Gen2DFingerprint(Chem.MolFromSmiles(smi), self.factory)
            self.assertEqual(len(sig), 39972)
            bs = tuple(sig.GetOnBits())
            self.assertEqual(len(bs), len(tgt))
            self.assertEqual(bs, tgt)

    def testOrderBug(self):
        sdFile = os.path.join(RDConfig.RDCodeDir, 'Chem', 'Pharm2D', 'test_data', 'orderBug.sdf')
        suppl = Chem.SDMolSupplier(sdFile)
        m1 = next(suppl)
        m2 = next(suppl)
        sig1 = Generate.Gen2DFingerprint(m1, self.factory)
        sig2 = Generate.Gen2DFingerprint(m2, self.factory)
        self.assertEqual(sig1, sig2)

    def testOrderBug2(self):
        from rdkit.Chem import Randomize
        from rdkit import DataStructs
        probes = ['Oc1nc(Oc2ncccc2)ccc1']
        for smi in probes:
            m1 = Chem.MolFromSmiles(smi)
            # m1.Debug()
            sig1 = Generate.Gen2DFingerprint(m1, self.factory)
            csmi = Chem.MolToSmiles(m1)
            m2 = Chem.MolFromSmiles(csmi)
            # m2.Debug()
            sig2 = Generate.Gen2DFingerprint(m2, self.factory)
            self.assertTrue(list(sig1.GetOnBits()) == list(sig2.GetOnBits()), f'{smi} {csmi}')
            self.assertEqual(DataStructs.DiceSimilarity(sig1, sig2), 1.0)
            self.assertEqual(sig1, sig2)
            for _ in range(10):
                m2 = Randomize.RandomizeMol(m1)
                sig2 = Generate.Gen2DFingerprint(m2, self.factory)
                if sig2 != sig1:
                    Generate._verbose = True
                    print('----------------')
                    sig1 = Generate.Gen2DFingerprint(m1, self.factory)
                    print('----------------')
                    sig2 = Generate.Gen2DFingerprint(m2, self.factory)
                    print('----------------')
                    print(Chem.MolToMolBlock(m1))
                    print('----------------')
                    print(Chem.MolToMolBlock(m2))
                    print('----------------')
                    s1 = set(sig1.GetOnBits())
                    s2 = set(sig2.GetOnBits())
                    print(s1.difference(s2))
                self.assertEqual(sig1, sig2)

    def testBitInfo(self):
        m = Chem.MolFromSmiles('OCC=CC(=O)O')
        bi = {}
        sig = Generate.Gen2DFingerprint(m, Gobbi_Pharm2D.factory, bitInfo=bi)
        self.assertEqual(sig.GetNumOnBits(), len(bi))
        self.assertEqual(list(sig.GetOnBits()), sorted(bi.keys()))
        self.assertEqual(sorted(bi.keys()), [23, 30, 150, 154, 157, 185, 28878, 30184])
        self.assertEqual(sorted(bi[28878]), [[(0, ), (5, ), (6, )]])
        self.assertEqual(sorted(bi[157]), [[(0, ), (6, )], [(5, ), (0, )]])


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
