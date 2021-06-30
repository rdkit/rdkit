#
#  Copyright (C) 2007 Greg Landrum
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import os
import unittest

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Subshape import SubshapeBuilder, SubshapeAligner, SubshapeObjects

debug = False
if not debug:
    SubshapeAligner.logger.setLevel(SubshapeAligner.RDLogger.ERROR)


class TestCase(unittest.TestCase):

    def test1(self):
        # computeCanonicalTransform returns more approximate eigenvalues/eigencvectors
        # when built against the native RDKit PowerEigenSolver, so unit test results
        # differ slightly
        builtAgainstEigen3 = hasattr(AllChem, 'ComputePrincipalAxesAndMomentsFromGyrationMatrix')
        if builtAgainstEigen3:
            expectedSkelPts = 15
            expectedAlgs = [0, 5, 21, 0]
            prunedAlgs = [0, 4, 11, 0]
        else:
            expectedSkelPts = 16
            expectedAlgs = [0, 5, 28, 0]
            prunedAlgs = [0, 4, 12, 0]
        filename = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                'test_data', '5ht3ligs.sdf')
        suppl = Chem.SDMolSupplier(filename)
        builder = SubshapeBuilder.SubshapeBuilder()
        builder.gridDims = (20., 20., 10)
        builder.gridSpacing = 0.5
        builder.winRad = 4.

        ms = []
        shapes = []
        for m in suppl:
            m = Chem.AddHs(m, addCoords=True)
            AllChem.CanonicalizeConformer(m.GetConformer())
            ms.append(m)
            shape = builder(m, terminalPtsOnly=True)
            shapes.append(shape)

        self.assertEqual(len(ms), 4)
        self.assertEqual(len(shapes), 4)
        self.assertEqual([len(x.skelPts) for x in shapes], [5, 5, 5, 5])

        refShape = builder.GenerateSubshapeShape(ms[0])

        self.assertEqual(len(refShape.skelPts), expectedSkelPts)

        aligner = SubshapeAligner.SubshapeAligner()
        aligner.shapeDistTol = .30

        algStore = []
        for i, s1 in enumerate(shapes):
            if not i or not s1:
                algStore.append([])
                continue
            m1 = ms[i]
            alignments = aligner.GetSubshapeAlignments(ms[0], refShape, m1, s1, builder)
            algStore.append(alignments)
        self.assertEqual([len(x) for x in algStore], expectedAlgs)

        algStore = []
        for i, s1 in enumerate(shapes):
            if not i or not s1:
                algStore.append([])
                continue
            m1 = ms[i]
            alignments = list(aligner(ms[0], refShape, m1, s1, builder))
            algStore.append(alignments)
        self.assertEqual([len(x) for x in algStore], expectedAlgs)

        pruned = []
        for i, mi in enumerate(ms):
            alignments = algStore[i]
            pruned.append(SubshapeAligner.ClusterAlignments(
                mi, alignments, builder, neighborTol=0.15))
        self.assertEqual([len(x) for x in pruned], prunedAlgs)


class TestSubshapeObjects(unittest.TestCase):

    def test_SubshapeShape(self):
        shape = SubshapeObjects.SubshapeShape()
        self.assertEqual(shape.shapes, [])


if __name__ == '__main__':
    unittest.main()
