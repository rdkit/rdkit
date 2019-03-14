# $Id$
#
#  Copyright (C) 2005  greg Landrum and Rational Discovery LLC
#   All Rights Reserved
#
import gzip
import os
import unittest

from rdkit import RDConfig
from rdkit.DataStructs import ExplicitBitVect
from rdkit.DataStructs.VectCollection import VectCollection
from rdkit.ML import InfoTheory
from rdkit.ML.DecTree.BuildSigTree import BuildSigTree, _GenerateRandomEnsemble
from rdkit.ML.DecTree.SigTree import SigTreeNode
from rdkit.TestRunner import redirect_stdout
from io import StringIO


class TestCase(unittest.TestCase):

    def setUp(self):
        t1 = SigTreeNode(None, 'root', 0)

        t2 = SigTreeNode(t1, 'nodeL1', 1)
        t1.AddChildNode(t2)
        t3 = SigTreeNode(t2, 'nodeLTerm0', 0, isTerminal=1)
        t4 = SigTreeNode(t2, 'nodeLTerm1', 1, isTerminal=1)
        t2.AddChildNode(t3)
        t2.AddChildNode(t4)

        t2 = SigTreeNode(t1, 'nodeR1', 2)
        t1.AddChildNode(t2)
        t3 = SigTreeNode(t2, 'nodeRTerm0', 1, isTerminal=1)
        t4 = SigTreeNode(t2, 'nodeRTerm1', 0, isTerminal=1)
        t2.AddChildNode(t3)
        t2.AddChildNode(t4)
        self.tree = t1

    def test1(self):
        t1 = self.tree
        bv = ExplicitBitVect(5)

        ex = ['nm', bv]
        self.assertFalse(t1.ClassifyExample(ex))
        bv.SetBit(1)
        self.assertTrue(t1.ClassifyExample(ex))

        bv.SetBit(0)
        self.assertTrue(t1.ClassifyExample(ex))

        bv.SetBit(2)
        self.assertFalse(t1.ClassifyExample(ex))

    def test2(self):
        t1 = self.tree
        vc = VectCollection()

        bv = ExplicitBitVect(5)
        bv.SetBitsFromList([0])
        vc.AddVect(1, bv)

        bv = ExplicitBitVect(5)
        bv.SetBitsFromList([1, 2])
        vc.AddVect(2, bv)

        ex = ['nm', bv, 1]
        self.assertTrue(t1.ClassifyExample(ex))

        bv = ExplicitBitVect(5)
        bv.SetBitsFromList([0, 2])
        vc.AddVect(1, bv)
        ex = ['nm', bv, 1]
        self.assertFalse(t1.ClassifyExample(ex))

    def test3(self):
        examples = []

        bv = ExplicitBitVect(2)
        vc = VectCollection()
        vc.AddVect(1, bv)
        examples.append(['a', vc, 1])

        bv = ExplicitBitVect(2)
        bv.SetBit(1)
        vc = VectCollection()
        vc.AddVect(1, bv)
        examples.append(['c', vc, 0])

        bv = ExplicitBitVect(2)
        bv.SetBit(1)
        vc = VectCollection()
        vc.AddVect(1, bv)
        examples.append(['c2', vc, 0])

        bv = ExplicitBitVect(2)
        bv.SetBit(0)
        vc = VectCollection()
        vc.AddVect(1, bv)
        examples.append(['d', vc, 0])

        bv = ExplicitBitVect(2)
        bv.SetBit(0)
        vc = VectCollection()
        vc.AddVect(1, bv)
        bv = ExplicitBitVect(2)
        bv.SetBit(1)
        vc.AddVect(2, bv)
        examples.append(['d2', vc, 0])

        bv = ExplicitBitVect(2)
        bv.SetBit(0)
        bv.SetBit(1)
        vc = VectCollection()
        vc.AddVect(1, bv)
        examples.append(['d', vc, 1])

        bv = ExplicitBitVect(2)
        bv.SetBit(0)
        bv.SetBit(1)
        vc = VectCollection()
        vc.AddVect(1, bv)
        examples.append(['e', vc, 1])

        f = StringIO()
        with redirect_stdout(f):
            t = BuildSigTree(examples, 2, metric=InfoTheory.InfoType.ENTROPY,
                             maxDepth=2, verbose=True)
        self.assertIn('Build', f.getvalue())

        self.assertEqual(t.GetName(), 'Bit-0')
        self.assertEqual(t.GetLabel(), 0)
        c0 = t.GetChildren()[0]
        self.assertEqual(c0.GetName(), 'Bit-1')
        self.assertEqual(c0.GetLabel(), 1)
        c1 = t.GetChildren()[1]
        self.assertEqual(c1.GetName(), 'Bit-1')
        self.assertEqual(c1.GetLabel(), 1)

        bv = ExplicitBitVect(2)
        bv.SetBit(0)
        vc = VectCollection()
        vc.AddVect(1, bv)
        bv = ExplicitBitVect(2)
        bv.SetBit(1)
        vc.AddVect(2, bv)
        r = t.ClassifyExample(['t', vc, 0])
        self.assertEqual(r, 0)

    def test4(self):
        import pickle
        gz = gzip.open(
          os.path.join(RDConfig.RDCodeDir, 'ML', 'DecTree', 'test_data', 'cdk2-few.pkl.gz'), 'rb')
        examples = pickle.load(gz, encoding='Latin1')
        t = BuildSigTree(examples, 2, maxDepth=3)
        self.assertEqual(t.GetLabel(), 2181)
        self.assertEqual(t.GetChildren()[0].GetLabel(), 2861)
        self.assertEqual(t.GetChildren()[1].GetLabel(), 8182)

    def test_GenerateRandomEnsemble(self):
        ensemble = _GenerateRandomEnsemble(2, 4)
        self.assertEqual(len(ensemble), 2)
        self.assertTrue(all(r < 4 for r in ensemble))

        ensemble = _GenerateRandomEnsemble(4, 4)
        self.assertEqual(len(ensemble), 4)
        self.assertTrue(all(r < 4 for r in ensemble))

        ensemble = _GenerateRandomEnsemble(4, 40)
        self.assertEqual(len(ensemble), 4)
        self.assertTrue(all(r < 40 for r in ensemble))


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
