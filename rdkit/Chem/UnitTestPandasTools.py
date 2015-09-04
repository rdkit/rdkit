from __future__ import print_function
import os
import unittest

from rdkit.six.moves import cStringIO as StringIO

from rdkit import RDConfig

from  rdkit.Chem import PandasTools
import numpy

methane = """\
Methane
     RDKit          

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
> <prop1>
12.34

> <prop2>
qwe

$$$$
"""

peroxide = """\
Peroxide
     RDKit          

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
> <prop2>
rtz

> <prop3>
yxcv

$$$$
"""


class TestLoadSDF(unittest.TestCase):
    def setUp(self):
        self.gz_filename = os.path.join(RDConfig.RDCodeDir, 'Chem/test_data', 'pandas_load.sdf.gz')
        
    # the doctest tests loading from a ".sdf" file so there's no need for that test here
    
    def test_load_gzip_file(self):
        df = PandasTools.LoadSDF(self.gz_filename)
        self.assertEqual(len(df), 13)
        # The molecule with index 1 is invalid, so it should be missing form the index
        self.assertEqual(list(df.index), [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])

    def test_load_from_sio(self):
        sio = StringIO(methane + peroxide)
        df = PandasTools.LoadSDF(sio)
        self.assertEqual(len(df), 2)
        self.assertEqual(list(df["ID"]), ["Methane", "Peroxide"])
        atom_counts = [mol.GetNumAtoms() for mol in df["ROMol"]]
        self.assertEqual(atom_counts, [1, 2])
    
    def test_load_specify_column_names(self):
        sio = StringIO(methane + peroxide)
        df = PandasTools.LoadSDF(sio, idName="CorpID", molColName="_rdmol")
        self.assertEqual(len(df), 2)
        self.assertEqual(list(df["CorpID"]), ["Methane", "Peroxide"])
        atom_counts = [mol.GetNumAtoms() for mol in df["_rdmol"]]
        self.assertEqual(atom_counts, [1, 2])

    def test_empty_file(self):
        # Should return an empty data frame with no rows or columns
        sio = StringIO()
        df = PandasTools.LoadSDF(sio)
        self.assertEqual(len(df), 0)
        self.assertEqual(len(df.index), 0)

    def test_passed_in_file_is_not_closed(self):
        sio = StringIO(methane)
        df = PandasTools.LoadSDF(sio)
        self.assertEqual(len(df), 1)
        self.assertFalse(sio.closed)

    def test_properties(self):
        sio = StringIO(peroxide + methane)
        df = PandasTools.LoadSDF(sio)
        self.assertEqual(set(df.columns), set("ROMol ID prop1 prop2 prop3".split()))
        prop1 = list(df["prop1"])
        self.assertTrue(numpy.isnan(prop1[0]), prop1[0])
        self.assertEqual(prop1[1], "12.34")

        self.assertEqual(list(df["prop2"]), ["rtz", "qwe"])
        
        prop3 = list(df["prop3"])
        self.assertEqual(prop3[0], "yxcv")
        self.assertTrue(numpy.isnan(prop3[1]), prop3[1])


if __name__ == '__main__':
  unittest.main()
