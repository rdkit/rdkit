#
#  Copyright (C) 2003-2017  Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" unit testing code for molecule suppliers

"""
import os
import tempfile
import unittest

from rdkit import Chem, RDLogger
from rdkit import RDConfig


class TestCase(unittest.TestCase):

    def tearDown(self):
        RDLogger.EnableLog('rdApp.error')

    def test1SDSupplier(self):
        fileN = os.path.join(RDConfig.RDCodeDir, 'VLib', 'NodeLib', 'test_data', 'NCI_aids.10.sdf')

        suppl = Chem.SDMolSupplier(fileN)
        ms = [x for x in suppl]
        self.assertEqual(len(ms), 10)

        # test repeating:
        ms = [x for x in suppl]
        self.assertEqual(len(ms), 10)

        # test reset:
        suppl.reset()
        m = next(suppl)
        self.assertEqual(m.GetProp('_Name'), '48')
        self.assertEqual(m.GetProp('NSC'), '48')
        self.assertEqual(m.GetProp('CAS_RN'), '15716-70-8')
        m = next(suppl)
        self.assertEqual(m.GetProp('_Name'), '78')
        self.assertEqual(m.GetProp('NSC'), '78')
        self.assertEqual(m.GetProp('CAS_RN'), '6290-84-2')

        suppl.reset()
        for _ in range(10):
            m = next(suppl)

        with self.assertRaises(StopIteration):
            m = next(suppl)

    def test2SmilesSupplier(self):
        fileN = os.path.join(RDConfig.RDCodeDir, 'VLib', 'NodeLib', 'test_data', 'pgp_20.txt')

        suppl = Chem.SmilesMolSupplier(
            fileN, delimiter='\t', smilesColumn=2, nameColumn=1, titleLine=1)
        ms = [x for x in suppl]
        self.assertEqual(len(ms), 20)

        # test repeating:
        ms = [x for x in suppl]
        self.assertEqual(len(ms), 20)
        # test reset:
        suppl.reset()
        m = next(suppl)
        self.assertEqual(m.GetProp('_Name'), 'ALDOSTERONE')
        self.assertEqual(m.GetProp('ID'), 'RD-PGP-0001')
        m = next(suppl)
        self.assertEqual(m.GetProp('_Name'), 'AMIODARONE')
        self.assertEqual(m.GetProp('ID'), 'RD-PGP-0002')
        suppl.reset()
        for _ in range(20):
            m = next(suppl)
        with self.assertRaises(StopIteration):
            m = next(suppl)

    def test3SmilesSupplier(self):
        txt = """C1CC1,1
CC(=O)O,3
fail,4
CCOC,5
"""
        RDLogger.DisableLog('rdApp.error')

        try:
            with tempfile.NamedTemporaryFile('w+', suffix='.csv', delete=False) as tmp:
                tmp.write(txt)
            suppl = Chem.SmilesMolSupplier(tmp.name, delimiter=',', smilesColumn=0, nameColumn=1,
                                           titleLine=0)
            ms = [x for x in suppl]
            suppl = None
            while ms.count(None):
                ms.remove(None)
            self.assertEqual(len(ms), 3)
        finally:
            os.unlink(tmp.name)


if __name__ == '__main__':
    unittest.main()
