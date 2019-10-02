"""unit testing code for the DbResultSet object

"""

import random
import unittest

from rdkit import RDConfig, Chem
from rdkit.Chem.Suppliers.DbMolSupplier import ForwardDbMolSupplier, RandomAccessDbMolSupplier
from rdkit.Chem.Suppliers.MolSupplier import MolSupplier
from rdkit.Dbase.DbConnection import DbConnect
from rdkit.Dbase.DbResultSet import DbResultSet, RandomAccessDbResultSet


class TestCase(unittest.TestCase):

    def setUp(self):
        self.dbName = RDConfig.RDTestDatabase
        self.conn = DbConnect(self.dbName)
        self.curs = self.conn.GetCursor()

    def test_MolSupplier(self):
        self.assertRaises(ValueError, MolSupplier)

    def test_general(self):
        # Check for a molecule column
        cmd = 'select * from ten_elements'
        results = DbResultSet(self.curs, self.conn, cmd)
        self.assertRaises(ValueError, ForwardDbMolSupplier, results)

    def test_ForwardDbMolSupplier(self):
        cmd = 'select * from simple_mols order by ID'
        results = DbResultSet(self.curs, self.conn, cmd)
        expected = list(results)

        results = DbResultSet(self.curs, self.conn, cmd)
        supp = ForwardDbMolSupplier(results)
        self.assertEqual(supp.GetColumnNames(), ('ID',))

        for smiles, mol in zip(expected, supp):
            self.assertEqual(Chem.MolToSmiles(Chem.MolFromSmiles(smiles[0])), Chem.MolToSmiles(mol))
            self.assertEqual(smiles[1], mol.GetProp('ID'))
        self.assertRaises(StopIteration, next, supp)

        # We can not use an index for ForwardDbMolSupplier
        with self.assertRaises(TypeError):
            supp[0]

    def test_RandomAccessDbMolSupplier(self):
        cmd = 'select * from simple_mols order by ID'
        results = RandomAccessDbResultSet(self.curs, self.conn, cmd)
        expected = list(results)

        results = RandomAccessDbResultSet(self.curs, self.conn, cmd)
        supp = RandomAccessDbMolSupplier(results)
        self.assertEqual(len(supp), len(expected))
        self.assertEqual(supp.GetColumnNames(), ('ID',))
        for smiles, mol in zip(expected, supp):
            self.assertEqual(Chem.MolToSmiles(Chem.MolFromSmiles(smiles[0])), Chem.MolToSmiles(mol))
            self.assertEqual(smiles[1], mol.GetProp('ID'))

        # Check that we can randomly access the data
        indices = list(range(len(expected)))
        random.shuffle(indices)
        for idx in indices:
            smiles = expected[idx]
            mol = supp[idx]
            self.assertEqual(Chem.MolToSmiles(Chem.MolFromSmiles(smiles[0])), Chem.MolToSmiles(mol))
            self.assertEqual(smiles[1], mol.GetProp('ID'))

        # We get an error if we access outside of the permitted range
        with self.assertRaises(IndexError):
            supp[len(expected)]

        # The DbMolSupplier doesn't support negative indices
        with self.assertRaises(IndexError):
            supp[-1]


if __name__ == '__main__':
    unittest.main()
