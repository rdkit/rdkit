#
#  Copyright (C) 2002-2008  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import os.path
import unittest
from io import StringIO

from rdkit import Chem, RDConfig
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D, Matcher, SigFactory
from rdkit.TestRunner import redirect_stdout


class TestCase(unittest.TestCase):

    def setUp(self):
        fdefFile = os.path.join(RDConfig.RDCodeDir, 'Chem', 'Pharm2D',
                                'test_data', 'BaseFeatures.fdef')
        featFactory = ChemicalFeatures.BuildFeatureFactory(fdefFile)
        self.factory = SigFactory.SigFactory(featFactory, minPointCount=2, maxPointCount=3)
        self.factory.SetBins([(0, 2), (2, 5), (5, 8)])
        self.factory.Init()

    def test1_simple(self):
        mol = Chem.MolFromSmiles('OCC(=O)CCCN')
        self.factory.skipFeats = ['Donor']
        self.factory.Init()
        self.assertEqual(self.factory.GetSigSize(), 510)
        Generate._verbose = False
        sig = Generate.Gen2DFingerprint(mol, self.factory)
        Generate._verbose = False
        tgt = (1, 2, 11, 52, 117)
        onBits = sig.GetOnBits()
        self.assertEqual(tuple(onBits), tgt)
        self.assertEqual(len(onBits), len(tgt))

        bitMatches = ([((0, ), (3, ))],
                      [((0, ), (7, )), ((3, ), (7, ))],
                      [((0, ), (3, ), (7, ))], )
        for i, bit in enumerate(onBits):
            matches = Matcher.GetAtomsMatchingBit(self.factory, bit, mol)
            # print bit,matches
            # tgt = bitMatches[i]
            # self.assertEqual(matches,tgt)

    def test2Bug28(self):
        smi = r'Cc([s]1)nnc1SCC(\CS2)=C(/C([O-])=O)N3C(=O)[C@H]([C@@H]23)NC(=O)C[n]4cnnn4'
        mol = Chem.MolFromSmiles(smi)
        factory = Gobbi_Pharm2D.factory
        factory.SetBins([(2, 3), (3, 4), (4, 5), (5, 8), (8, 100)])
        sig = Generate.Gen2DFingerprint(mol, factory)
        onBits = sig.GetOnBits()
        for bit in onBits:
            atoms = Matcher.GetAtomsMatchingBit(factory, bit, mol, justOne=1)
            self.assertTrue(len(atoms))

    def test3Roundtrip(self):
        # longer-running Bug 28 test
        nToDo = 20
        with open(os.path.join(RDConfig.RDDataDir, 'NCI', 'first_5K.smi'), 'r') as inF:
            inD = inF.readlines()[:nToDo]
        factory = Gobbi_Pharm2D.factory
        factory.SetBins([(2, 3), (3, 4), (4, 5), (5, 8), (8, 100)])
        for line in inD:
            smi = line.split('\t')[0]
            mol = Chem.MolFromSmiles(smi)
            sig = Generate.Gen2DFingerprint(mol, factory)
            onBits = sig.GetOnBits()
            for bit in onBits:
                atoms = Matcher.GetAtomsMatchingBit(factory, bit, mol, justOne=1)
                assert len(atoms), f'bit {bit} failed to match for smi {smi}'

    def test_exampleCode(self):
        # We make sure that the example code runs
        f = StringIO()
        with redirect_stdout(f):
            Matcher._exampleCode()
        self.assertIn('finished', f.getvalue())


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
