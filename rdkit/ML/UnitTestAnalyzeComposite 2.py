# $Id$
#
#  Copyright (C) 2004-2006  Greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""unit testing code for the AnalyzeComposite functionality

"""
import os
import unittest

from rdkit import RDConfig
from rdkit.ML import AnalyzeComposite
import pickle


class TestCase(unittest.TestCase):

    def setUp(self):
        self.baseDir = os.path.join(RDConfig.RDCodeDir, 'ML', 'test_data')

    def test1_Issue163(self):
        name1 = os.path.join(self.baseDir, 'humanoral.1.pkl')
        try:
            with open(name1, 'rb') as pklF:
                c1 = pickle.load(pklF)
        except Exception:
            c1 = None
        self.assertTrue(c1)
        name2 = os.path.join(self.baseDir, 'humanoral.2.pkl')
        try:
            with open(name2, 'rb') as pklF:
                c2 = pickle.load(pklF)
        except Exception:
            c2 = None
        self.assertTrue(c2)

        try:
            res = sorted(AnalyzeComposite.ProcessIt([c1, c2], verbose=-1))
        except Exception:
            import traceback
            traceback.print_exc()
            ok = 0
        else:
            ok = 1
        self.assertTrue(ok)

        self.assertEqual(res[0][0], 'BALABANJ')
        self.assertEqual(res[1][0], 'BERTZCT')
        self.assertEqual(res[-1][0], 'VSA_ESTATE9')
        for entry in res:
            self.assertEqual(len(entry), 5)


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
