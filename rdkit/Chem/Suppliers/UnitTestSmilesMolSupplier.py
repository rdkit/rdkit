# $Id$
#
#  Copyright (C) 2003-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""
unit testing code for the Smiles file handling stuff
"""
import unittest

from rdkit import Chem

from rdkit import RDLogger


class TestCase(unittest.TestCase):

    def setUp(self):
        self.smis = ['CC', 'CCC', 'CCCCC', 'CCCCCC', 'CCCCCCC', 'CC', 'CCCCOC']
        self.nMolecules = len(self.smis)

    def tearDown(self):
        RDLogger.EnableLog('rdApp.error')

    def assertMolecule(self, mol, i, msg=''):
        """ Assert that we have a valid molecule """
        self.assertIsNotNone(mol, '{0}read {1} failed'.format(msg, i))
        self.assertGreater(mol.GetNumAtoms(), 0, '{0}no atoms in mol {1}'.format(msg, i))

    def test_SmilesReaderIndex(self):
        # tests lazy reads
        supp = Chem.SmilesMolSupplierFromText('\n'.join(self.smis), ',', 0, -1, 0)
        for i in range(4):
            self.assertMolecule(next(supp), i)

        i = len(supp) - 1
        self.assertMolecule(supp[i], i)

        # Use in a list comprehension
        ms = [Chem.MolToSmiles(mol) for mol in supp]
        self.assertEqual(ms, self.smis)

        self.assertEqual(len(supp), self.nMolecules, 'bad supplier length')

        # Despite iterating through the whole supplier, we can still access by index
        i = self.nMolecules - 3
        self.assertMolecule(supp[i - 1], i, msg='back index: ')

        with self.assertRaises(IndexError):
            _ = supp[self.nMolecules]  # out of bound read must fail

        # and we can access with negative numbers
        mol1 = supp[len(supp) - 1]
        mol2 = supp[-1]
        self.assertEqual(Chem.MolToSmiles(mol1), Chem.MolToSmiles(mol2))

    def test_SmilesReaderIterator(self):
        # tests lazy reads using the iterator interface "
        supp = Chem.SmilesMolSupplierFromText('\n'.join(self.smis), ',', 0, -1, 0)

        nDone = 0
        for mol in supp:
            self.assertMolecule(mol, nDone)
            nDone += 1
        self.assertEqual(nDone, self.nMolecules, 'bad number of molecules')

        self.assertEqual(len(supp), self.nMolecules, 'bad supplier length')

        # Despite iterating through the whole supplier, we can still access by index
        i = self.nMolecules - 3
        self.assertMolecule(supp[i - 1], i, msg='back index: ')

        with self.assertRaises(IndexError):
            _ = supp[self.nMolecules]  # out of bound read must not fail

    def test_SmilesReaderBoundaryConditions(self):
        # Suppress the error message due to the incorrect smiles
        RDLogger.DisableLog('rdApp.error')

        smis = ['CC', 'CCOC', 'fail', 'CCO']
        supp = Chem.SmilesMolSupplierFromText('\n'.join(smis), ',', 0, -1, 0)
        self.assertEqual(len(supp), 4)
        self.assertIsNone(supp[2])
        self.assertIsNotNone(supp[3])

        supp = Chem.SmilesMolSupplierFromText('\n'.join(smis), ',', 0, -1, 0)
        self.assertIsNone(supp[2])
        self.assertIsNotNone(supp[3])
        self.assertEqual(len(supp), 4)
        with self.assertRaises(IndexError):
            supp[4]

        supp = Chem.SmilesMolSupplierFromText('\n'.join(smis), ',', 0, -1, 0)
        self.assertEqual(len(supp), 4)
        self.assertIsNotNone(supp[3])
        with self.assertRaises(IndexError):
            supp[4]

        supp = Chem.SmilesMolSupplierFromText('\n'.join(smis), ',', 0, -1, 0)
        with self.assertRaises(IndexError):
            supp[4]

        self.assertEqual(len(supp), 4)
        self.assertIsNotNone(supp[3])


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
