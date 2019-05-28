# $Id$
#
#  Copyright (C) 2007  greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#


import doctest
import unittest

from rdkit import Chem, RDLogger
from rdkit.VLib.NodeLib import SDSupply, SmartsMolFilter, SmartsRemover
from rdkit.VLib.NodeLib import SmilesDupeFilter, SmilesOutput, SmilesSupply
from rdkit.VLib.Supply import SupplyNode
from io import StringIO


def load_tests(loader, tests, ignore):
    """ Add the Doctests from the module """
    tests.addTests(doctest.DocTestSuite(SDSupply, optionflags=doctest.ELLIPSIS))
    tests.addTests(doctest.DocTestSuite(SmartsMolFilter, optionflags=doctest.ELLIPSIS))
    tests.addTests(doctest.DocTestSuite(SmartsRemover, optionflags=doctest.ELLIPSIS))
    tests.addTests(doctest.DocTestSuite(SmilesDupeFilter, optionflags=doctest.ELLIPSIS))
    tests.addTests(doctest.DocTestSuite(SmilesOutput, optionflags=doctest.ELLIPSIS))
    tests.addTests(doctest.DocTestSuite(SmilesSupply, optionflags=doctest.ELLIPSIS))
    #   tests.addTests(doctest.DocTestSuite(DbMolSupply, optionflags=doctest.ELLIPSIS))
    return tests


class Test_NodeLib(unittest.TestCase):

    def tearDown(self):
        RDLogger.EnableLog('rdApp.error')

    def test_SmartsMolFilter(self):
        smis = ['C1CCC1', 'C1CCC1C=O', 'CCCC', 'CCC=O', 'CC(=O)C', 'CCN', 'NCCN', 'NCC=O']
        mols = [Chem.MolFromSmiles(x) for x in smis]
        suppl = SupplyNode(contents=mols)
        self.assertEqual(len(list(suppl)), 8)

        smas = ['C=O', 'CN']
        counts = [1, 2]
        filt = SmartsMolFilter.SmartsFilter(patterns=smas, counts=counts)
        filt.AddParent(suppl)
        self.assertEqual(len(list(filt)), 5)

        suppl.reset()
        filt.SetNegate(True)
        self.assertEqual(len(list(filt)), 3)

        smas = ['C=O', 'CN']
        filt = SmartsMolFilter.SmartsFilter(patterns=smas)
        filt.AddParent(suppl)
        self.assertEqual(len(list(filt)), 6)

        self.assertRaises(ValueError, SmartsMolFilter.SmartsFilter, patterns=smas,
                          counts=['notEnough', ])
        RDLogger.DisableLog('rdApp.error')
        self.assertRaises(ValueError, SmartsMolFilter.SmartsFilter, patterns=['BadSmarts'])
        RDLogger.EnableLog('rdApp.error')

    def test_SmilesOutput(self):
        smis = ['C1CCC1', 'C1CC1', 'C=O', 'CCN']
        mols = [Chem.MolFromSmiles(x) for x in smis]
        for i, mol in enumerate(mols, 100):
            mol.SetProp('ID', str(i))

        suppl1 = SupplyNode(contents=mols)
        suppl2 = SupplyNode(contents='abcd')

        sio = StringIO()
        node = SmilesOutput.OutputNode(idField='ID', dest=sio, delim=', ')
        node.AddParent(suppl1)
        node.AddParent(suppl2)
        list(node)
        self.assertEqual(
            sio.getvalue(), '100, C1CCC1, a\n101, C1CC1, b\n102, C=O, c\n103, CCN, d\n')

    def test_SmartsRemover(self):
        salts = ['[Cl;H1&X1,-]', '[Na+]', '[O;H2,H1&-,X0&-2]', 'BadSmarts']
        RDLogger.DisableLog('rdApp.error')
        self.assertRaises(ValueError, SmartsRemover.SmartsRemover, patterns=salts)
        RDLogger.EnableLog('rdApp.error')

    def test_SmilesDupeFilter(self):
        smis = ['C1CCC1', 'CCCC', 'CCCC', 'C1CCC1']
        mols = [Chem.MolFromSmiles(x) for x in smis]
        suppl = SupplyNode(contents=mols)
        self.assertEqual(len(list(suppl)), 4)

        dupFilter = SmilesDupeFilter.DupeFilter()
        dupFilter.AddParent(suppl)
        self.assertEqual(len(list(dupFilter)), 2)


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
