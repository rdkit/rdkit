from __future__ import print_function

import doctest
import gzip
import os
import shutil
import tempfile
import unittest

import numpy

from rdkit.Chem import PandasTools
from rdkit.six import PY3, StringIO, BytesIO
from rdkit import RDConfig, rdBase, Chem

try:
  import IPython
except ImportError:
  IPython = None

# We make sure that we don't mess up the Mol methods for the rest of the tests
PandasTools.UninstallPandasTools()


@unittest.skipIf(PandasTools.pd is None, 'Pandas not installed, skipping')
class TestPandasTools(unittest.TestCase):

  def __init__(self, methodName='runTest'):
    self.df = getTestFrame()
    self.df.index.name = 'IndexName'
    super(TestPandasTools, self).__init__(methodName=methodName)

  def setUp(self):
    PandasTools.InstallPandasTools()
    PandasTools.ChangeMoleculeRendering(renderer='PNG')
    self._molRepresentation = PandasTools.molRepresentation
    self._highlightSubstructures = PandasTools.highlightSubstructures

  def tearDown(self):
    PandasTools.molRepresentation = self._molRepresentation
    PandasTools.highlightSubstructures = self._highlightSubstructures
    PandasTools.UninstallPandasTools()

  def testDoctest(self):
    # We need to do it like this to ensure that default RDkit functionality is restored
    failed, _ = doctest.testmod(PandasTools,
                                optionflags=doctest.ELLIPSIS + doctest.NORMALIZE_WHITESPACE)
    self.assertFalse(failed)

  def test_RestoreMonkeyPatch(self):
    sio = getStreamIO(methane + peroxide)
    df = PandasTools.LoadSDF(sio)
    html = df.to_html()
    self.assertIn('data:image/png;base64', html)
    self.assertIn('table', html)

    PandasTools.UninstallPandasTools()
    html = df.to_html()
    self.assertNotIn('data:image/png;base64', html)
    self.assertIn('rdkit.Chem.rdchem.Mol', html)
    self.assertIn('table', html)

    PandasTools.InstallPandasTools()
    html = df.to_html()
    self.assertIn('data:image/png;base64', html)
    self.assertIn('table', html)

    PandasTools.UninstallPandasTools()
    html = df.to_html()
    self.assertNotIn('data:image/png;base64', html)
    self.assertIn('rdkit.Chem.rdchem.Mol', html)
    self.assertIn('table', html)

  def test_FrameToGridImage(self):
    # This test only makes sure that we get no exception. To see the created images, set
    # interactive to True
    interactive = False
    self.assertTrue(True)
    df = self.df

    result = PandasTools.FrameToGridImage(df)
    if interactive:
      result.show()

    result = PandasTools.FrameToGridImage(df, legendsCol='PUBCHEM_IUPAC_INCHIKEY')
    if interactive:
      result.show()

    result = PandasTools.FrameToGridImage(df, legendsCol=df.index.name)
    if interactive:
      result.show()

  def test_AddMurckoToFrame(self):
    df = self.df.copy()
    self.assertIn('ROMol', df.columns)
    self.assertNotIn('Murcko_SMILES', df.columns)
    PandasTools.AddMurckoToFrame(df)
    self.assertIn('ROMol', df.columns)
    self.assertIn('Murcko_SMILES', df.columns)
    self.assertEqual(df['Murcko_SMILES'][10], 'O=C(CCn1c(-c2ccccc2)n[nH]c1=S)Nc1ccccn1')

    PandasTools.AddMurckoToFrame(df, Generic=True)
    self.assertIn('ROMol', df.columns)
    self.assertIn('Murcko_SMILES', df.columns)
    self.assertEqual(df['Murcko_SMILES'][10], 'CC(CCC1C(C)CCC1C1CCCCC1)CC1CCCCC1')

  def test_SaveSMILESFromFrame(self):
    sio = StringIO()
    PandasTools.SaveSMILESFromFrame(self.df, sio)
    result = sio.getvalue()
    self.assertIn(self.df['SMILES'][10], result)
    self.assertIn(self.df['ID'][10], result)

    sio = StringIO()
    PandasTools.SaveSMILESFromFrame(self.df, sio, NamesCol='PUBCHEM_IUPAC_INCHIKEY')
    result = sio.getvalue()
    self.assertIn(self.df['SMILES'][10], result)
    self.assertIn(self.df['PUBCHEM_IUPAC_INCHIKEY'][10], result)

  @unittest.skipIf(IPython is None, 'Package IPython required for testing')
  def test_svgRendering(self):
    df = PandasTools.LoadSDF(getStreamIO(methane + peroxide))
    self.assertIn('image/png', str(df))
    self.assertNotIn('svg', str(df))

    PandasTools.molRepresentation = 'svg'
    self.assertIn('svg', str(df))
    self.assertNotIn('image/png', str(df))

    # we can use upper case for the molRepresentation
    PandasTools.molRepresentation = 'PNG'
    self.assertNotIn('svg', str(df))
    self.assertIn('image/png', str(df))

  def test_patchHeadFrame(self):
    df = self.df.copy()
    result = str(df.head())
    self.assertIn('35024984', result)
    self.assertNotIn('35024985', result)

  def test_AddMoleculeColumnToFrame(self):
    df = PandasTools.LoadSDF(
      getStreamIO(methane + peroxide), isomericSmiles=True, smilesName='Smiles')
    PandasTools.ChangeMoleculeRendering(frame=df, renderer='String')
    del df['ROMol']
    self.assertNotIn('ROMol', str(df))
    PandasTools.AddMoleculeColumnToFrame(df, includeFingerprints=False)
    self.assertIn('ROMol', str(df))

  def test_molge(self):
    # We want to have the default RDkit functionality for testing
    PandasTools.UninstallPandasTools()
    molge = PandasTools._molge
    mol1 = Chem.MolFromSmiles('CCC')
    mol2 = Chem.MolFromSmiles('CC')
    mol3 = Chem.MolFromSmiles('CN')

    self.assertFalse(molge(mol1, None))
    self.assertFalse(molge(None, mol1))

    self.assertFalse(hasattr(mol1, '_substructfp'))
    self.assertFalse(hasattr(mol2, '_substructfp'))
    self.assertFalse(hasattr(mol1, '__sssAtoms'))
    self.assertFalse(hasattr(mol2, '__sssAtoms'))

    self.assertTrue(molge(mol1, mol2))
    self.assertEqual(mol1.__dict__['__sssAtoms'], [0, 1])
    PandasTools.highlightSubstructures = False
    self.assertTrue(molge(mol1, mol2))
    self.assertEqual(mol1.__dict__['__sssAtoms'], [])

    PandasTools.highlightSubstructures = True
    self.assertFalse(molge(mol2, mol1))
    self.assertEqual(mol2.__dict__['__sssAtoms'], [])

    self.assertFalse(molge(mol1, mol3))
    self.assertEqual(mol1.__dict__['__sssAtoms'], [])


@unittest.skipIf(PandasTools.pd is None, 'Pandas not installed, skipping')
class TestLoadSDF(unittest.TestCase):
  gz_filename = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', 'pandas_load.sdf.gz')

  # the doctest tests loading from a ".sdf" file so there's no need for that test here

  def test_load_gzip_file(self):
    rdBase.DisableLog('rdApp.error')
    df = PandasTools.LoadSDF(self.gz_filename)
    rdBase.EnableLog('rdApp.error')
    self.assertEqual(len(df), 13)
    # The molecule with index 1 is invalid, so it should be missing form the index
    self.assertEqual(list(df.index), [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])

  def test_load_from_sio(self):
    sio = getStreamIO(methane + peroxide)
    df = PandasTools.LoadSDF(sio)
    self.assertEqual(len(df), 2)
    self.assertEqual(list(df["ID"]), ["Methane", "Peroxide"])
    atom_counts = [mol.GetNumAtoms() for mol in df["ROMol"]]
    self.assertEqual(atom_counts, [1, 2])

  def test_load_specify_column_names(self):
    sio = getStreamIO(methane + peroxide)
    df = PandasTools.LoadSDF(sio, idName="CorpID", molColName="_rdmol")
    self.assertEqual(len(df), 2)
    self.assertEqual(list(df["CorpID"]), ["Methane", "Peroxide"])
    atom_counts = [mol.GetNumAtoms() for mol in df["_rdmol"]]
    self.assertEqual(atom_counts, [1, 2])

  def test_empty_file(self):
    # Should return an empty data frame with no rows or columns
    sio = getStreamIO(None)
    df = PandasTools.LoadSDF(sio)
    self.assertEqual(len(df), 0)
    self.assertEqual(len(df.index), 0)

  def test_passed_in_file_is_not_closed(self):
    sio = getStreamIO(methane)
    df = PandasTools.LoadSDF(sio)
    self.assertEqual(len(df), 1)
    self.assertFalse(sio.closed)

  def test_properties(self):
    sio = getStreamIO(peroxide + methane)
    df = PandasTools.LoadSDF(sio)
    self.assertEqual(set(df.columns), set("ROMol ID prop1 prop2 prop3".split()))
    prop1 = list(df["prop1"])
    self.assertTrue(numpy.isnan(prop1[0]), prop1[0])
    self.assertEqual(prop1[1], "12.34")

    self.assertEqual(list(df["prop2"]), ["rtz", "qwe"])

    prop3 = list(df["prop3"])
    self.assertEqual(prop3[0], "yxcv")
    self.assertTrue(numpy.isnan(prop3[1]), prop3[1])

  def test_ignore_mol_column(self):
    sio = getStreamIO(peroxide + methane)
    df = PandasTools.LoadSDF(sio, molColName=None)
    self.assertEqual(set(df.columns), set("ID prop1 prop2 prop3".split()))


@unittest.skipIf(PandasTools.pd is None, 'Pandas not installed, skipping')
class TestWriteSDF(unittest.TestCase):

  def setUp(self):
    self.df = PandasTools.LoadSDF(getStreamIO(methane + peroxide))

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

  def test_specify_numeric_column_2(self):
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
      with open(filename) as f:
        s = f.read()
      self.assertEqual(s.count("\n$$$$\n"), 2)
      self.assertEqual(s.split("\n", 1)[0], "Methane")
    finally:
      shutil.rmtree(dirname)

  def test_write_to_sdf_gz(self):
    dirname = tempfile.mkdtemp()
    try:
      filename = os.path.join(dirname, "test.sdf.gz")
      PandasTools.WriteSDF(self.df, filename)
      with gzip.open(filename) as f:
        s = f.read()
      if PY3:
        s = s.decode('utf-8')
      s = s.replace(os.linesep, '\n')
      self.assertEqual(s.count("\n$$$$\n"), 2)
      self.assertEqual(s.split("\n", 1)[0], "Methane")
    finally:
      shutil.rmtree(dirname)


def getStreamIO(sdfString):
  """ Return a StringIO/BytesIO for the string """
  if PY3:
    sio = BytesIO() if sdfString is None else BytesIO(sdfString.encode('utf-8'))
  else:  # pragma: nocover
    sio = StringIO() if sdfString is None else StringIO(sdfString)
  return sio


def getTestFrame():
  rdBase.DisableLog('rdApp.error')
  sdfFile = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', 'pandas_load.sdf.gz')
  df = PandasTools.LoadSDF(sdfFile, smilesName='SMILES')
  rdBase.EnableLog('rdApp.error')
  return df


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

if __name__ == '__main__':  # pragma: nocover
  unittest.main()
