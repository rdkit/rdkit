from __future__ import print_function
import os
import unittest

from rdkit.six.moves import cStringIO as StringIO
from rdkit import RDConfig

from rdkit.Chem import PandasTools
import numpy
import tempfile, shutil
import gzip

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

class TestWriteSDF(unittest.TestCase):
    def setUp(self):
        sio = StringIO(methane + peroxide)
        self.df = PandasTools.LoadSDF(sio)
        
    def test_default_write_does_not_include_tags(self):
        sio = StringIO()
        PandasTools.WriteSDF(self.df, sio)
        s = sio.getvalue()
        self.assertNotIn(s, "prop2")
        
    def test_identifier_from_a_column(self):
        sio = StringIO()
        PandasTools.WriteSDF(self.df, sio, idName="prop2")
        s = sio.getvalue()
        first_line = s.split("\n", 1)[0]
        self.assertEqual(first_line, "qwe")

    def test_all_numeric_with_no_numeric_columns(self):
        sio = StringIO()
        PandasTools.WriteSDF(self.df, sio, allNumeric=True)
        s = sio.getvalue()
        self.assertFalse(">" in s, s)
        self.assertNotIn("7\n\n", s)  # double-check that the numeric tests don't pass by accident
        self.assertNotIn("8\n\n", s)

    def test_all_numeric_with_numeric_columns(self):
        sio = StringIO()
        df = self.df
        df["len"] = df["ID"].map(len)
        PandasTools.WriteSDF(df, sio, allNumeric=True)
        s = sio.getvalue()
        self.assertEqual(s.count("<len>"), 2)
        self.assertIn("7\n\n", s)
        self.assertIn("8\n\n", s)

    def test_specify_numeric_column(self):
        sio = StringIO()
        df = self.df
        df["len2"] = df["ID"].map(len)
        PandasTools.WriteSDF(df, sio, properties=["len2"])
        s = sio.getvalue()
        self.assertEqual(s.count("<len2>"), 2)
        self.assertIn("7\n\n", s)
        self.assertIn("8\n\n", s)

    def test_specify_numeric_column(self):
        sio = StringIO()
        df = self.df
        df["len2"] = df["ID"].map(len)
        df["len3"] = df["len2"].map(float)
        PandasTools.WriteSDF(df, sio, properties=["len2", "len3"])
        s = sio.getvalue()
        self.assertEqual(s.count("<len2>"), 2)
        self.assertEqual(s.count("<len3>"), 2)
        self.assertIn("7\n\n", s)
        self.assertIn("7.0\n\n", s)
        self.assertIn("8\n\n", s)
        self.assertIn("8.0\n\n", s)
        
    def test_write_to_sdf(self):
        dirname = tempfile.mkdtemp()
        try:
            filename = os.path.join(dirname, "test.sdf")
            PandasTools.WriteSDF(self.df, filename)
            s = open(filename, "U").read()
            self.assertEqual(s.count("\n$$$$\n"), 2)
            self.assertEqual(s.split("\n", 1)[0], "Methane")
        finally:
            shutil.rmtree(dirname)

    def test_write_to_sdf_gz(self):
        dirname = tempfile.mkdtemp()
        try:
            filename = os.path.join(dirname, "test.sdf.gz")
            PandasTools.WriteSDF(self.df, filename)
            s = gzip.open(filename).read()
            self.assertEqual(s.count("\n$$$$\n"), 2)
            self.assertEqual(s.split("\n", 1)[0], "Methane")
        finally:
            shutil.rmtree(dirname)
                    
if __name__ == '__main__':
    from rdkit.six import PY3
    if not PY3: # FIX: The StringIO tests fail on python3
        unittest.main()
