#
#  Copyright (C) 2003-2023  Greg Landrum and other RDKit contributors
#         All Rights Reserved
#
""" This is a rough coverage test of the python wrapper

it's intended to be shallow, but broad

"""

import doctest
import gc
import gzip
import importlib.util
import logging
import os
import pickle
import sys
import tempfile
import unittest
from contextlib import contextmanager
from datetime import datetime, timedelta
from io import BytesIO, StringIO

import rdkit.Chem.rdDepictor
from rdkit import Chem, DataStructs, RDConfig, __version__, rdBase
from rdkit.Chem import rdqueries
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold

# Boost functions are NOT found by doctest, this "fixes" them
#  by adding the doctests to a fake module
spec = importlib.util.spec_from_loader("TestReplaceCore", loader=None)
TestReplaceCore = importlib.util.module_from_spec(spec)
code = """
from rdkit.Chem import ReplaceCore
def ReplaceCore(*a, **kw):
    '''%s
    '''
    return Chem.ReplaceCore(*a, **kw)
""" % "\n".join([x.lstrip() for x in Chem.ReplaceCore.__doc__.split("\n")])
exec(code, TestReplaceCore.__dict__)


@contextmanager
def log_to_python(level=None):
  """
    Temporarily redirect logging to Python streams, optionally
    setting a specific log level.
  """
  rdBase.LogToPythonLogger()
  pylog = logging.getLogger("rdkit")
  if level is not None:
    original_level = pylog.level
    pylog.setLevel(level)

  yield pylog

  if level is not None:
    pylog.setLevel(original_level)
  rdBase.LogToCppStreams()


@contextmanager
def capture_logging(level=None):
  """
    Temporarily redirect logging to a Python StringIO, optionally
    setting a specific log level.
  """
  log_stream = StringIO()
  stream_handler = logging.StreamHandler(stream=log_stream)

  with log_to_python(level) as pylog:
    pylog.addHandler(stream_handler)

    yield log_stream

    pylog.removeHandler(stream_handler)


def load_tests(loader, tests, ignore):
  tests.addTests(doctest.DocTestSuite(TestReplaceCore))
  return tests


def feq(v1, v2, tol2=1e-4):
  return abs(v1 - v2) <= tol2


def getTotalFormalCharge(mol):
  totalFormalCharge = 0
  for atom in mol.GetAtoms():
    totalFormalCharge += atom.GetFormalCharge()
  return totalFormalCharge


def cmpFormalChargeBondOrder(self, mol1, mol2):
  self.assertEqual(mol1.GetNumAtoms(), mol2.GetNumAtoms())
  self.assertEqual(mol1.GetNumBonds(), mol2.GetNumBonds())
  for i in range(mol1.GetNumAtoms()):
    self.assertEqual(
      mol1.GetAtomWithIdx(i).GetFormalCharge(),
      mol2.GetAtomWithIdx(i).GetFormalCharge())
  for i in range(mol1.GetNumBonds()):
    self.assertEqual(mol1.GetBondWithIdx(i).GetBondType(), mol2.GetBondWithIdx(i).GetBondType())


def setResidueFormalCharge(mol, res, fc):
  for query in res:
    matches = mol.GetSubstructMatches(query)
    for match in matches:
      mol.GetAtomWithIdx(match[-1]).SetFormalCharge(fc)


def getBtList2(resMolSuppl):
  btList2 = []
  while (not resMolSuppl.atEnd()):
    resMol = next(resMolSuppl)
    bt = []
    for bond in resMol.GetBonds():
      bt.append(int(bond.GetBondTypeAsDouble()))
    btList2.append(bt)
  for i in range(len(btList2)):
    same = True
    for j in range(len(btList2[i])):
      if (not i):
        continue
      if (same):
        same = (btList2[i][j] == btList2[i - 1][j])
    if (i and same):
      return None
  return btList2


class TestCase(unittest.TestCase):

  def test0Except(self):

    with self.assertRaises(IndexError):
      Chem.tossit()

  def test1Table(self):

    tbl = Chem.GetPeriodicTable()
    self.assertTrue(tbl)

    self.assertTrue(feq(tbl.GetAtomicWeight(6), 12.011))
    self.assertTrue(feq(tbl.GetAtomicWeight("C"), 12.011))
    self.assertTrue(tbl.GetAtomicNumber('C') == 6)
    self.assertTrue(feq(tbl.GetRvdw(6), 1.7))
    self.assertTrue(feq(tbl.GetRvdw("C"), 1.7))
    self.assertTrue(feq(tbl.GetRcovalent(6), 0.76))
    self.assertTrue(feq(tbl.GetRcovalent("C"), 0.76))
    self.assertTrue(tbl.GetDefaultValence(6) == 4)
    self.assertTrue(tbl.GetDefaultValence("C") == 4)
    self.assertTrue(tuple(tbl.GetValenceList(6)) == (4, ))
    self.assertTrue(tuple(tbl.GetValenceList("C")) == (4, ))
    self.assertTrue(tuple(tbl.GetValenceList(16)) == (2, 4, 6))
    self.assertTrue(tuple(tbl.GetValenceList("S")) == (2, 4, 6))
    self.assertTrue(tbl.GetNOuterElecs(6) == 4)
    self.assertTrue(tbl.GetNOuterElecs("C") == 4)
    self.assertTrue(tbl.GetMostCommonIsotope(6) == 12)
    self.assertTrue(tbl.GetMostCommonIsotope('C') == 12)
    self.assertTrue(tbl.GetMostCommonIsotopeMass(6) == 12.0)
    self.assertTrue(tbl.GetMostCommonIsotopeMass('C') == 12.0)
    self.assertTrue(tbl.GetAbundanceForIsotope(6, 12) == 98.93)
    self.assertTrue(tbl.GetAbundanceForIsotope('C', 12) == 98.93)
    self.assertTrue(feq(tbl.GetRb0(6), 0.77))
    self.assertTrue(feq(tbl.GetRb0("C"), 0.77))
    self.assertTrue(tbl.GetElementSymbol(6) == 'C')
    self.assertTrue(tbl.GetElementName(6) == 'Carbon')

  def test2Atom(self):
    atom = Chem.Atom(6)
    self.assertTrue(atom)
    self.assertTrue(atom.GetAtomicNum() == 6)
    atom.SetAtomicNum(8)
    self.assertTrue(atom.GetAtomicNum() == 8)

    atom = Chem.Atom("C")
    self.assertTrue(atom)
    self.assertTrue(atom.GetAtomicNum() == 6)

  def test3Bond(self):
    # No longer relevant, bonds are not constructible from Python
    pass

  def test4Mol(self):
    mol = Chem.Mol()
    self.assertTrue(mol)

  def test5Smiles(self):
    mol = Chem.MolFromSmiles('n1ccccc1')
    self.assertTrue(mol)
    self.assertTrue(mol.GetNumAtoms() == 6)
    self.assertTrue(mol.GetNumAtoms(1) == 6)
    self.assertTrue(mol.GetNumAtoms(0) == 11)
    at = mol.GetAtomWithIdx(2)
    self.assertTrue(at.GetAtomicNum() == 6)
    at = mol.GetAtomWithIdx(0)
    self.assertTrue(at.GetAtomicNum() == 7)

  def _test6Bookmarks(self):
    mol = Chem.MolFromSmiles('n1ccccc1')
    self.assertTrue(mol)

    self.assertTrue(not mol.HasAtomBookmark(0))
    mol.SetAtomBookmark(mol.GetAtomWithIdx(0), 0)
    mol.SetAtomBookmark(mol.GetAtomWithIdx(1), 1)
    self.assertTrue(mol.HasAtomBookmark(0))
    self.assertTrue(mol.HasAtomBookmark(1))

    if 1:
      self.assertTrue(not mol.HasBondBookmark(0))
      self.assertTrue(not mol.HasBondBookmark(1))
      mol.SetBondBookmark(mol.GetBondWithIdx(0), 0)
      mol.SetBondBookmark(mol.GetBondWithIdx(1), 1)
      self.assertTrue(mol.HasBondBookmark(0))
      self.assertTrue(mol.HasBondBookmark(1))

    at = mol.GetAtomWithBookmark(0)
    self.assertTrue(at)
    self.assertTrue(at.GetAtomicNum() == 7)
    mol.ClearAtomBookmark(0)
    self.assertTrue(not mol.HasAtomBookmark(0))
    self.assertTrue(mol.HasAtomBookmark(1))
    mol.ClearAllAtomBookmarks()
    self.assertTrue(not mol.HasAtomBookmark(0))
    self.assertTrue(not mol.HasAtomBookmark(1))

    mol.SetAtomBookmark(mol.GetAtomWithIdx(1), 1)

    if 1:
      self.assertTrue(mol.HasBondBookmark(0))
      self.assertTrue(mol.HasBondBookmark(1))
      bond = mol.GetBondWithBookmark(0)
      self.assertTrue(bond)
      mol.ClearBondBookmark(0)
      self.assertTrue(not mol.HasBondBookmark(0))
      self.assertTrue(mol.HasBondBookmark(1))
      mol.ClearAllBondBookmarks()
      self.assertTrue(not mol.HasBondBookmark(0))
      self.assertTrue(not mol.HasBondBookmark(1))

      self.assertTrue(mol.HasAtomBookmark(1))

  def test7Atom(self):
    mol = Chem.MolFromSmiles('n1ccccc1C[CH2-]')
    self.assertTrue(mol)
    Chem.SanitizeMol(mol)
    a0 = mol.GetAtomWithIdx(0)
    a1 = mol.GetAtomWithIdx(1)
    a6 = mol.GetAtomWithIdx(6)
    a7 = mol.GetAtomWithIdx(7)

    self.assertTrue(a0.GetAtomicNum() == 7)
    self.assertTrue(a0.GetSymbol() == 'N')
    self.assertTrue(a0.GetIdx() == 0)

    aList = [a0, a1, a6, a7]
    self.assertTrue(a0.GetDegree() == 2)
    self.assertTrue(a1.GetDegree() == 2)
    self.assertTrue(a6.GetDegree() == 2)
    self.assertTrue(a7.GetDegree() == 1)
    self.assertTrue([x.GetDegree() for x in aList] == [2, 2, 2, 1])

    self.assertTrue([x.GetTotalNumHs() for x in aList] == [0, 1, 2, 2])
    self.assertTrue([x.GetNumImplicitHs() for x in aList] == [0, 1, 2, 0])
    self.assertTrue([x.GetExplicitValence() for x in aList] == [3, 3, 2, 3])
    self.assertTrue([x.GetImplicitValence() for x in aList] == [0, 1, 2, 0])
    self.assertTrue([x.GetFormalCharge() for x in aList] == [0, 0, 0, -1])
    self.assertTrue([x.GetNoImplicit() for x in aList] == [0, 0, 0, 1])
    self.assertTrue([x.GetNumExplicitHs() for x in aList] == [0, 0, 0, 2])
    self.assertTrue([x.GetIsAromatic() for x in aList] == [1, 1, 0, 0])
    self.assertTrue([x.GetHybridization() for x in aList] == [
      Chem.HybridizationType.SP2, Chem.HybridizationType.SP2, Chem.HybridizationType.SP3,
      Chem.HybridizationType.SP3
    ], [x.GetHybridization() for x in aList])

  def test8Bond(self):
    mol = Chem.MolFromSmiles('n1ccccc1CC(=O)O')
    self.assertTrue(mol)
    Chem.SanitizeMol(mol)
    # note bond numbering is funny because of ring closure
    b0 = mol.GetBondWithIdx(0)
    b6 = mol.GetBondWithIdx(6)
    b7 = mol.GetBondWithIdx(7)
    b8 = mol.GetBondWithIdx(8)

    bList = [b0, b6, b7, b8]
    self.assertTrue(
      [x.GetBondType() for x in bList] ==
      [Chem.BondType.AROMATIC, Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.SINGLE])
    self.assertTrue([x.GetIsAromatic() for x in bList] == [1, 0, 0, 0])
    self.assertEqual(bList[0].GetBondTypeAsDouble(), 1.5)
    self.assertEqual(bList[1].GetBondTypeAsDouble(), 1.0)
    self.assertEqual(bList[2].GetBondTypeAsDouble(), 2.0)

    self.assertTrue([x.GetIsConjugated() != 0 for x in bList] == [1, 0, 1, 1],
                    [x.GetIsConjugated() != 0 for x in bList])
    self.assertTrue([x.GetBeginAtomIdx() for x in bList] == [0, 6, 7, 7],
                    [x.GetBeginAtomIdx() for x in bList])
    self.assertTrue([x.GetBeginAtom().GetIdx() for x in bList] == [0, 6, 7, 7])
    self.assertTrue([x.GetEndAtomIdx() for x in bList] == [1, 7, 8, 9])
    self.assertTrue([x.GetEndAtom().GetIdx() for x in bList] == [1, 7, 8, 9])

  def test9Smarts(self):
    query1 = Chem.MolFromSmarts('C(=O)O')
    self.assertTrue(query1)
    query2 = Chem.MolFromSmarts('C(=O)[O,N]')
    self.assertTrue(query2)
    query3 = Chem.MolFromSmarts('[$(C(=O)O)]')
    self.assertTrue(query3)

    mol = Chem.MolFromSmiles('CCC(=O)O')
    self.assertTrue(mol)

    self.assertTrue(mol.HasSubstructMatch(query1))
    self.assertTrue(mol.HasSubstructMatch(query2))
    self.assertTrue(mol.HasSubstructMatch(query3))

    mol = Chem.MolFromSmiles('CCC(=O)N')
    self.assertTrue(mol)

    self.assertTrue(not mol.HasSubstructMatch(query1))
    self.assertTrue(mol.HasSubstructMatch(query2))
    self.assertTrue(not mol.HasSubstructMatch(query3))

  def test10Iterators(self):
    mol = Chem.MolFromSmiles('CCOC')
    self.assertTrue(mol)

    for atom in mol.GetAtoms():
      self.assertTrue(atom)
    ats = mol.GetAtoms()
    ats[1]
    with self.assertRaisesRegex(IndexError, ""):
      ats[12]

    for bond in mol.GetBonds():
      self.assertTrue(bond)
    bonds = mol.GetBonds()
    bonds[1]
    with self.assertRaisesRegex(IndexError, ""):
      bonds[12]

  def test11MolOps(self):
    mol = Chem.MolFromSmiles('C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3')
    self.assertTrue(mol)
    smi = Chem.MolToSmiles(mol)
    Chem.SanitizeMol(mol)
    nr = Chem.GetSymmSSSR(mol)
    self.assertTrue((len(nr) == 3))
    nr = Chem.GetSSSR(mol)
    self.assertTrue((len(nr) == 3))

  def test12Smarts(self):
    query1 = Chem.MolFromSmarts('C(=O)O')
    self.assertTrue(query1)
    query2 = Chem.MolFromSmarts('C(=O)[O,N]')
    self.assertTrue(query2)
    query3 = Chem.MolFromSmarts('[$(C(=O)O)]')
    self.assertTrue(query3)

    mol = Chem.MolFromSmiles('CCC(=O)O')
    self.assertTrue(mol)

    self.assertTrue(mol.HasSubstructMatch(query1))
    self.assertTrue(mol.GetSubstructMatch(query1) == (2, 3, 4))
    self.assertTrue(mol.HasSubstructMatch(query2))
    self.assertTrue(mol.GetSubstructMatch(query2) == (2, 3, 4))
    self.assertTrue(mol.HasSubstructMatch(query3))
    self.assertTrue(mol.GetSubstructMatch(query3) == (2, ))

    mol = Chem.MolFromSmiles('CCC(=O)N')
    self.assertTrue(mol)

    self.assertTrue(not mol.HasSubstructMatch(query1))
    self.assertTrue(not mol.GetSubstructMatch(query1))
    self.assertTrue(mol.HasSubstructMatch(query2))
    self.assertTrue(mol.GetSubstructMatch(query2) == (2, 3, 4))
    self.assertTrue(not mol.HasSubstructMatch(query3))

    mol = Chem.MolFromSmiles('OC(=O)CC(=O)O')
    self.assertTrue(mol)
    self.assertTrue(mol.HasSubstructMatch(query1))
    self.assertTrue(mol.GetSubstructMatch(query1) == (1, 2, 0))
    self.assertTrue(mol.GetSubstructMatches(query1) == ((1, 2, 0), (4, 5, 6)))
    self.assertTrue(mol.HasSubstructMatch(query2))
    self.assertTrue(mol.GetSubstructMatch(query2) == (1, 2, 0))
    self.assertTrue(mol.GetSubstructMatches(query2) == ((1, 2, 0), (4, 5, 6)))
    self.assertTrue(mol.HasSubstructMatch(query3))
    self.assertTrue(mol.GetSubstructMatches(query3) == ((1, ), (4, )))

  def test13Smarts(self):
    # previous smarts problems:
    query = Chem.MolFromSmarts('N(=,-C)')
    self.assertTrue(query)
    mol = Chem.MolFromSmiles('N#C')
    self.assertTrue(not mol.HasSubstructMatch(query))
    mol = Chem.MolFromSmiles('N=C')
    self.assertTrue(mol.HasSubstructMatch(query))
    mol = Chem.MolFromSmiles('NC')
    self.assertTrue(mol.HasSubstructMatch(query))

    query = Chem.MolFromSmarts('[Cl,$(O)]')
    mol = Chem.MolFromSmiles('C(=O)O')
    self.assertTrue(len(mol.GetSubstructMatches(query)) == 2)
    mol = Chem.MolFromSmiles('C(=N)N')
    self.assertTrue(len(mol.GetSubstructMatches(query)) == 0)

    query = Chem.MolFromSmarts('[$([O,S]-[!$(*=O)])]')
    mol = Chem.MolFromSmiles('CC(S)C(=O)O')
    self.assertTrue(len(mol.GetSubstructMatches(query)) == 1)
    mol = Chem.MolFromSmiles('C(=O)O')
    self.assertTrue(len(mol.GetSubstructMatches(query)) == 0)

  def test14Hs(self):
    m = Chem.MolFromSmiles('CC(=O)[OH]')
    self.assertEqual(m.GetNumAtoms(), 4)
    m2 = Chem.AddHs(m)
    self.assertEqual(m2.GetNumAtoms(), 8)
    m2 = Chem.RemoveHs(m2)
    self.assertEqual(m2.GetNumAtoms(), 4)

    m = Chem.MolFromSmiles('CC[H]', False)
    self.assertEqual(m.GetNumAtoms(), 3)
    m2 = Chem.MergeQueryHs(m)
    self.assertEqual(m2.GetNumAtoms(), 2)
    self.assertTrue(m2.GetAtomWithIdx(1).HasQuery())

    m = Chem.MolFromSmiles('CC[H]', False)
    self.assertEqual(m.GetNumAtoms(), 3)
    m1 = Chem.RemoveHs(m)
    self.assertEqual(m1.GetNumAtoms(), 2)
    self.assertEqual(m1.GetAtomWithIdx(1).GetNumExplicitHs(), 0)
    m1 = Chem.RemoveHs(m, updateExplicitCount=True)
    self.assertEqual(m1.GetNumAtoms(), 2)
    self.assertEqual(m1.GetAtomWithIdx(1).GetNumExplicitHs(), 1)

    # test merging of mapped hydrogens
    m = Chem.MolFromSmiles('CC[H]', False)
    m.GetAtomWithIdx(2).SetProp("molAtomMapNumber", "1")
    self.assertEqual(m.GetNumAtoms(), 3)
    m2 = Chem.MergeQueryHs(m, mergeUnmappedOnly=True)
    self.assertTrue(m2 is not None)
    self.assertEqual(m2.GetNumAtoms(), 3)
    self.assertFalse(m2.GetAtomWithIdx(1).HasQuery())

    # here the hydrogen is unmapped
    #  should be the same as merging all hydrogens
    m = Chem.MolFromSmiles('CC[H]', False)
    m.GetAtomWithIdx(1).SetProp("molAtomMapNumber", "1")
    self.assertEqual(m.GetNumAtoms(), 3)
    m2 = Chem.MergeQueryHs(m, mergeUnmappedOnly=True)
    self.assertTrue(m2 is not None)
    self.assertEqual(m2.GetNumAtoms(), 2)
    self.assertTrue(m2.GetAtomWithIdx(1).HasQuery())

    # test merging of isotopes, by default deuterium will not be merged
    m = Chem.MolFromSmiles('CC[2H]', False)
    self.assertEqual(m.GetNumAtoms(), 3)
    m2 = Chem.MergeQueryHs(m)
    self.assertTrue(m2 is not None)
    self.assertEqual(m2.GetNumAtoms(), 3)
    self.assertFalse(m2.GetAtomWithIdx(1).HasQuery())

    # here deuterium is merged
    # should be the same as merging all hydrogens
    m = Chem.MolFromSmiles('CC[2H]', False)
    self.assertEqual(m.GetNumAtoms(), 3)
    m2 = Chem.MergeQueryHs(m, mergeIsotopes=True)
    self.assertTrue(m2 is not None)
    self.assertEqual(m2.GetNumAtoms(), 2)
    self.assertTrue(m2.GetAtomWithIdx(1).HasQuery())

    # test github758
    m = Chem.MolFromSmiles('CCC')
    self.assertEqual(m.GetNumAtoms(), 3)
    m = Chem.AddHs(m, onlyOnAtoms=(0, 2))
    self.assertEqual(m.GetNumAtoms(), 9)
    self.assertEqual(m.GetAtomWithIdx(0).GetDegree(), 4)
    self.assertEqual(m.GetAtomWithIdx(2).GetDegree(), 4)
    self.assertEqual(m.GetAtomWithIdx(1).GetDegree(), 2)

  def test15Neighbors(self):
    m = Chem.MolFromSmiles('CC(=O)[OH]')
    self.assertTrue(m.GetNumAtoms() == 4)

    a = m.GetAtomWithIdx(1)
    ns = a.GetNeighbors()
    self.assertTrue(len(ns) == 3)

    bs = a.GetBonds()
    self.assertTrue(len(bs) == 3)

    for b in bs:
      try:
        a2 = b.GetOtherAtom(a)
      except Exception:
        a2 = None
      self.assertTrue(a2)
    self.assertTrue(len(bs) == 3)

  def test16Pickle(self):
    import pickle
    m = Chem.MolFromSmiles('C1=CN=CC=C1')
    pkl = pickle.dumps(m)
    m2 = pickle.loads(pkl)
    self.assertTrue(type(m2) == Chem.Mol)
    smi1 = Chem.MolToSmiles(m)
    smi2 = Chem.MolToSmiles(m2)
    self.assertTrue(smi1 == smi2)

    pkl = pickle.dumps(Chem.RWMol(m))
    m2 = pickle.loads(pkl)
    self.assertTrue(type(m2) == Chem.RWMol)
    smi1 = Chem.MolToSmiles(m)
    smi2 = Chem.MolToSmiles(m2)
    self.assertTrue(smi1 == smi2)

  def test16Props(self):
    m = Chem.MolFromSmiles('C1=CN=CC=C1')
    self.assertTrue(not m.HasProp('prop1'))
    self.assertTrue(not m.HasProp('prop2'))
    self.assertTrue(not m.HasProp('prop2'))
    m.SetProp('prop1', 'foob')
    self.assertTrue(not m.HasProp('prop2'))
    self.assertTrue(m.HasProp('prop1'))
    self.assertTrue(m.GetProp('prop1') == 'foob')
    self.assertTrue(not m.HasProp('propo'))
    try:
      m.GetProp('prop2')
    except KeyError:
      ok = 1
    else:
      ok = 0
    self.assertTrue(ok)

    # test computed properties
    m.SetProp('cprop1', 'foo', 1)
    m.SetProp('cprop2', 'foo2', 1)

    m.ClearComputedProps()
    self.assertTrue(not m.HasProp('cprop1'))
    self.assertTrue(not m.HasProp('cprop2'))

    m.SetDoubleProp("a", 2.0)
    self.assertTrue(m.GetDoubleProp("a") == 2.0)

    try:
      self.assertTrue(m.GetIntProp("a") == 2.0)
      raise Exception("Expected runtime exception")
    except ValueError:
      pass

    try:
      self.assertTrue(m.GetUnsignedProp("a") == 2.0)
      raise Exception("Expected runtime exception")
    except ValueError:
      pass

    m.SetDoubleProp("a", -2)
    self.assertTrue(m.GetDoubleProp("a") == -2.0)
    m.SetIntProp("a", -2)
    self.assertTrue(m.GetIntProp("a") == -2)

    try:
      m.SetUnsignedProp("a", -2)
      raise Exception("Expected failure with negative unsigned number")
    except OverflowError:
      pass

    m.SetBoolProp("a", False)
    self.assertFalse(m.GetBoolProp("a"))

    self.assertEqual(m.GetPropsAsDict(), {'a': False, 'prop1': 'foob'})
    m.SetDoubleProp("b", 1000.0)
    m.SetUnsignedProp("c", 2000)
    m.SetIntProp("d", -2)
    m.SetUnsignedProp("e", 2, True)
    self.assertEqual(m.GetPropsAsDict(False, True), {
      'a': False,
      'c': 2000,
      'b': 1000.0,
      'e': 2,
      'd': -2,
      'prop1': 'foob'
    })
    m = Chem.MolFromSmiles('C1=CN=CC=C1')
    m.SetProp("int", "1000")
    m.SetProp("double", "10000.123")
    self.assertEqual(m.GetPropsAsDict(), {"int": 1000, "double": 10000.123})

    self.assertEqual(type(m.GetPropsAsDict()['int']), int)
    self.assertEqual(type(m.GetPropsAsDict()['double']), float)

  def test17Kekulize(self):
    m = Chem.MolFromSmiles('c1ccccc1')
    smi = Chem.MolToSmiles(m)
    self.assertTrue(smi == 'c1ccccc1')

    Chem.Kekulize(m)
    smi = Chem.MolToSmiles(m)
    self.assertTrue(smi == 'c1ccccc1')

    m = Chem.MolFromSmiles('c1ccccc1')
    smi = Chem.MolToSmiles(m)
    self.assertTrue(smi == 'c1ccccc1')

    Chem.Kekulize(m, 1)
    smi = Chem.MolToSmiles(m)
    self.assertTrue(smi == 'C1=CC=CC=C1', smi)

  def test18Paths(self):

    m = Chem.MolFromSmiles("C1CC2C1CC2")
    #self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==7)
    #print(Chem.FindAllPathsOfLengthN(m,3,useBonds=0))
    self.assertTrue(
      len(Chem.FindAllPathsOfLengthN(m, 2, useBonds=1)) == 10,
      Chem.FindAllPathsOfLengthN(m, 2, useBonds=1))
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m, 3, useBonds=1)) == 14)

    m = Chem.MolFromSmiles('C1CC1C')
    self.assertTrue(m)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m, 1, useBonds=1)) == 4)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m, 2, useBonds=1)) == 5)
    self.assertTrue(
      len(Chem.FindAllPathsOfLengthN(m, 3, useBonds=1)) == 3,
      Chem.FindAllPathsOfLengthN(m, 3, useBonds=1))
    self.assertTrue(
      len(Chem.FindAllPathsOfLengthN(m, 4, useBonds=1)) == 1,
      Chem.FindAllPathsOfLengthN(m, 4, useBonds=1))
    self.assertTrue(
      len(Chem.FindAllPathsOfLengthN(m, 5, useBonds=1)) == 0,
      Chem.FindAllPathsOfLengthN(m, 5, useBonds=1))

    #
    #  Hexane example from Hall-Kier Rev.Comp.Chem. paper
    #  Rev. Comp. Chem. vol 2, 367-422, (1991)
    #
    m = Chem.MolFromSmiles("CCCCCC")
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m, 1, useBonds=1)) == 5)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m, 2, useBonds=1)) == 4)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m, 3, useBonds=1)) == 3)

    m = Chem.MolFromSmiles("CCC(C)CC")
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m, 1, useBonds=1)) == 5)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m, 2, useBonds=1)) == 5)
    self.assertTrue(
      len(Chem.FindAllPathsOfLengthN(m, 3, useBonds=1)) == 4,
      Chem.FindAllPathsOfLengthN(m, 3, useBonds=1))

    m = Chem.MolFromSmiles("CCCC(C)C")
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m, 1, useBonds=1)) == 5)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m, 2, useBonds=1)) == 5)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m, 3, useBonds=1)) == 3)

    m = Chem.MolFromSmiles("CC(C)C(C)C")
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m, 1, useBonds=1)) == 5)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m, 2, useBonds=1)) == 6)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m, 3, useBonds=1)) == 4)

    m = Chem.MolFromSmiles("CC(C)(C)CC")
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m, 1, useBonds=1)) == 5)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m, 2, useBonds=1)) == 7)
    self.assertTrue(
      len(Chem.FindAllPathsOfLengthN(m, 3, useBonds=1)) == 3,
      Chem.FindAllPathsOfLengthN(m, 3, useBonds=1))

    m = Chem.MolFromSmiles("C1CCCCC1")
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m, 1, useBonds=1)) == 6)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m, 2, useBonds=1)) == 6)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m, 3, useBonds=1)) == 6)

    m = Chem.MolFromSmiles("C1CC2C1CC2")
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m, 1, useBonds=1)) == 7)
    self.assertTrue(
      len(Chem.FindAllPathsOfLengthN(m, 2, useBonds=1)) == 10,
      Chem.FindAllPathsOfLengthN(m, 2, useBonds=1))
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m, 3, useBonds=1)) == 14)

    m = Chem.MolFromSmiles("CC2C1CCC12")
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m, 1, useBonds=1)) == 7)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m, 2, useBonds=1)) == 11)
    # FIX: this result disagrees with the paper (which says 13),
    #   but it seems right
    self.assertTrue(
      len(Chem.FindAllPathsOfLengthN(m, 3, useBonds=1)) == 15,
      Chem.FindAllPathsOfLengthN(m, 3, useBonds=1))

  def test19Subgraphs(self):
    m = Chem.MolFromSmiles('C1CC1C')
    self.assertTrue(m)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 1, 0)) == 4)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 2)) == 5)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 3)) == 4)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 4)) == 1)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 5)) == 0)

    #
    #  Hexane example from Hall-Kier Rev.Comp.Chem. paper
    #  Rev. Comp. Chem. vol 2, 367-422, (1991)
    #
    m = Chem.MolFromSmiles("CCCCCC")
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 1)) == 5)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 2)) == 4)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 3)) == 3)

    l = Chem.FindAllSubgraphsOfLengthMToN(m, 1, 3)
    self.assertEqual(len(l), 3)
    self.assertEqual(len(l[0]), 5)
    self.assertEqual(len(l[1]), 4)
    self.assertEqual(len(l[2]), 3)
    self.assertRaises(ValueError, lambda: Chem.FindAllSubgraphsOfLengthMToN(m, 4, 3))

    m = Chem.MolFromSmiles("CCC(C)CC")
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 1)) == 5)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 2)) == 5)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 3)) == 5)

    m = Chem.MolFromSmiles("CCCC(C)C")
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 1)) == 5)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 2)) == 5)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 3)) == 4)

    m = Chem.MolFromSmiles("CC(C)C(C)C")
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 1)) == 5)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 2)) == 6)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 3)) == 6)

    m = Chem.MolFromSmiles("CC(C)(C)CC")
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 1)) == 5)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 2)) == 7)
    self.assertTrue(
      len(Chem.FindAllSubgraphsOfLengthN(m, 3)) == 7, Chem.FindAllSubgraphsOfLengthN(m, 3))

    m = Chem.MolFromSmiles("C1CCCCC1")
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 1)) == 6)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 2)) == 6)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 3)) == 6)
    #self.assertTrue(len(Chem.FindUniqueSubgraphsOfLengthN(m,1))==1)
    self.assertTrue(len(Chem.FindUniqueSubgraphsOfLengthN(m, 2)) == 1)
    self.assertTrue(len(Chem.FindUniqueSubgraphsOfLengthN(m, 3)) == 1)

    m = Chem.MolFromSmiles("C1CC2C1CC2")
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 1)) == 7)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 2)) == 10)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 3)) == 16)

    m = Chem.MolFromSmiles("CC2C1CCC12")
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 1)) == 7)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m, 2)) == 11)
    self.assertTrue(
      len(Chem.FindAllSubgraphsOfLengthN(m, 3)) == 18, len(Chem.FindAllSubgraphsOfLengthN(m, 3)))

  def test20IsInRing(self):
    m = Chem.MolFromSmiles('C1CCC1C')
    self.assertTrue(m)
    self.assertTrue(m.GetAtomWithIdx(0).IsInRingSize(4))
    self.assertTrue(m.GetAtomWithIdx(1).IsInRingSize(4))
    self.assertTrue(m.GetAtomWithIdx(2).IsInRingSize(4))
    self.assertTrue(m.GetAtomWithIdx(3).IsInRingSize(4))
    self.assertTrue(not m.GetAtomWithIdx(4).IsInRingSize(4))

    self.assertTrue(not m.GetAtomWithIdx(0).IsInRingSize(3))
    self.assertTrue(not m.GetAtomWithIdx(1).IsInRingSize(3))
    self.assertTrue(not m.GetAtomWithIdx(2).IsInRingSize(3))
    self.assertTrue(not m.GetAtomWithIdx(3).IsInRingSize(3))
    self.assertTrue(not m.GetAtomWithIdx(4).IsInRingSize(3))

    self.assertTrue(m.GetBondWithIdx(0).IsInRingSize(4))
    self.assertTrue(not m.GetBondWithIdx(3).IsInRingSize(4))
    self.assertTrue(not m.GetBondWithIdx(0).IsInRingSize(3))
    self.assertTrue(not m.GetBondWithIdx(3).IsInRingSize(3))

  def test21Robustification(self):
    ok = False
    # FIX: at the moment I can't figure out how to catch the
    # actual exception that BPL is throwing when it gets
    # invalid arguments (Boost.Python.ArgumentError)
    try:
      Chem.MolFromSmiles('C=O').HasSubstructMatch(Chem.MolFromSmarts('fiib'))
    #except ValueError:
    #  ok=True
    except Exception:
      ok = True
    self.assertTrue(ok)

  def test22DeleteSubstruct(self):
    query = Chem.MolFromSmarts('C(=O)O')
    mol = Chem.MolFromSmiles('CCC(=O)O')
    nmol = Chem.DeleteSubstructs(mol, query)

    self.assertTrue(Chem.MolToSmiles(nmol) == 'CC')

    mol = Chem.MolFromSmiles('CCC(=O)O.O=CO')
    # now delete only fragments
    nmol = Chem.DeleteSubstructs(mol, query, 1)
    self.assertTrue(Chem.MolToSmiles(nmol) == 'CCC(=O)O', Chem.MolToSmiles(nmol))

    mol = Chem.MolFromSmiles('CCC(=O)O.O=CO')
    nmol = Chem.DeleteSubstructs(mol, query, 0)
    self.assertTrue(Chem.MolToSmiles(nmol) == 'CC')

    mol = Chem.MolFromSmiles('CCCO')
    nmol = Chem.DeleteSubstructs(mol, query, 0)
    self.assertTrue(Chem.MolToSmiles(nmol) == 'CCCO')

    # Issue 96 prevented this from working:
    mol = Chem.MolFromSmiles('CCC(=O)O.O=CO')
    nmol = Chem.DeleteSubstructs(mol, query, 1)
    self.assertTrue(Chem.MolToSmiles(nmol) == 'CCC(=O)O')
    nmol = Chem.DeleteSubstructs(nmol, query, 1)
    self.assertTrue(Chem.MolToSmiles(nmol) == 'CCC(=O)O')
    nmol = Chem.DeleteSubstructs(nmol, query, 0)
    self.assertTrue(Chem.MolToSmiles(nmol) == 'CC')

  def test23MolFileParsing(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'triazine.mol')
    with open(fileN, 'r') as inF:
      inD = inF.read()
    m1 = Chem.MolFromMolBlock(inD)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 9)

    m1 = Chem.MolFromMolFile(fileN)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 9)

    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'triazine.mof')
    self.assertRaises(IOError, lambda: Chem.MolFromMolFile(fileN))

    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'list-query.mol')
    query = Chem.MolFromMolFile(fileN)
    smi = Chem.MolToSmiles(query)
    self.assertEqual(smi, '*1ccccc1')
    smi = Chem.MolToSmarts(query)
    self.assertEqual(smi, '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6,#7,#15]:1')
    smi = Chem.MolToSmarts(query, rootedAtAtom=5)
    self.assertEqual(smi, '[#6,#7,#15]1:[#6]:[#6]:[#6]:[#6]:[#6]:1')

    query = Chem.MolFromMolFile(fileN, sanitize=False)
    smi = Chem.MolToSmiles(query)
    self.assertEqual(smi, '*1=CC=CC=C1')
    query.UpdatePropertyCache()
    smi = Chem.MolToSmarts(query)
    self.assertEqual(smi, '[#6]1=[#6]-[#6]=[#6]-[#6]=[#6,#7,#15]-1')
    smi = Chem.MolToSmarts(query, rootedAtAtom=3)
    self.assertEqual(smi, '[#6]1=[#6]-[#6]=[#6]-[#6,#7,#15]=[#6]-1')
    smi = "C1=CC=CC=C1"
    mol = Chem.MolFromSmiles(smi, 0)
    self.assertTrue(mol.HasSubstructMatch(query))
    Chem.SanitizeMol(mol)
    self.assertTrue(not mol.HasSubstructMatch(query))

    mol = Chem.MolFromSmiles('N1=CC=CC=C1', 0)
    self.assertTrue(mol.HasSubstructMatch(query))
    mol = Chem.MolFromSmiles('S1=CC=CC=C1', 0)
    self.assertTrue(not mol.HasSubstructMatch(query))
    mol = Chem.MolFromSmiles('P1=CC=CC=C1', 0)
    self.assertTrue(mol.HasSubstructMatch(query))

    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'issue123.mol')
    mol = Chem.MolFromMolFile(fileN)
    self.assertTrue(mol)
    self.assertEqual(mol.GetNumAtoms(), 23)
    mol = Chem.MolFromMolFile(fileN, removeHs=False)
    self.assertTrue(mol)
    self.assertEqual(mol.GetNumAtoms(), 39)

  # test23 was for Chem.DaylightFingerprint, which is deprecated

  def test24RDKFingerprint(self):
    from rdkit import DataStructs
    m1 = Chem.MolFromSmiles('C1=CC=CC=C1')
    fp1 = Chem.RDKFingerprint(m1)
    self.assertTrue(len(fp1) == 2048)
    m2 = Chem.MolFromSmiles('C1=CC=CC=C1')
    fp2 = Chem.RDKFingerprint(m2)

    tmp = DataStructs.TanimotoSimilarity(fp1, fp2)
    self.assertTrue(tmp == 1.0, tmp)

    m2 = Chem.MolFromSmiles('C1=CC=CC=N1')
    fp2 = Chem.RDKFingerprint(m2)
    self.assertTrue(len(fp2) == 2048)
    tmp = DataStructs.TanimotoSimilarity(fp1, fp2)
    self.assertTrue(tmp < 1.0, tmp)
    self.assertTrue(tmp > 0.0, tmp)

    fp3 = Chem.RDKFingerprint(m1, tgtDensity=0.3)
    self.assertTrue(len(fp3) < 2048)

    m1 = Chem.MolFromSmiles('C1=CC=CC=C1')
    fp1 = Chem.RDKFingerprint(m1)
    m2 = Chem.MolFromSmiles('C1=CC=CC=N1')
    fp2 = Chem.RDKFingerprint(m2)
    self.assertNotEqual(fp1, fp2)

    atomInvariants = [1] * 6
    fp1 = Chem.RDKFingerprint(m1, atomInvariants=atomInvariants)
    fp2 = Chem.RDKFingerprint(m2, atomInvariants=atomInvariants)
    self.assertEqual(fp1, fp2)

    m2 = Chem.MolFromSmiles('C1CCCCN1')
    fp1 = Chem.RDKFingerprint(m1, atomInvariants=atomInvariants, useBondOrder=False)
    fp2 = Chem.RDKFingerprint(m2, atomInvariants=atomInvariants, useBondOrder=False)
    self.assertEqual(fp1, fp2)

    # rooted at atom
    m1 = Chem.MolFromSmiles('CCCCCO')
    fp1 = Chem.RDKFingerprint(m1, 1, 4, nBitsPerHash=1, fromAtoms=[0])
    self.assertEqual(fp1.GetNumOnBits(), 4)
    m1 = Chem.MolFromSmiles('CCCCCO')
    fp1 = Chem.RDKFingerprint(m1, 1, 4, nBitsPerHash=1, fromAtoms=[0, 5])
    self.assertEqual(fp1.GetNumOnBits(), 8)

    # test sf.net issue 270:
    fp1 = Chem.RDKFingerprint(m1, atomInvariants=[x.GetAtomicNum() + 10 for x in m1.GetAtoms()])

    # atomBits
    m1 = Chem.MolFromSmiles('CCCO')
    l = []
    fp1 = Chem.RDKFingerprint(m1, minPath=1, maxPath=2, nBitsPerHash=1, atomBits=l)
    self.assertEqual(fp1.GetNumOnBits(), 4)
    self.assertEqual(len(l), m1.GetNumAtoms())
    self.assertEqual(len(l[0]), 2)
    self.assertEqual(len(l[1]), 3)
    self.assertEqual(len(l[2]), 4)
    self.assertEqual(len(l[3]), 2)

  def test25SDMolSupplier(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'NCI_aids_few.sdf')
    #fileN = "../FileParsers/test_data/NCI_aids_few.sdf"
    sdSup = Chem.SDMolSupplier(fileN)
    molNames = [
      "48", "78", "128", "163", "164", "170", "180", "186", "192", "203", "210", "211", "213",
      "220", "229", "256"
    ]

    chgs192 = {8: 1, 11: 1, 15: -1, 18: -1, 20: 1, 21: 1, 23: -1, 25: -1}
    i = 0
    for mol in sdSup:
      self.assertTrue(mol)
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      i += 1
      if (mol.GetProp("_Name") == "192"):
        # test parsed charges on one of the molecules
        for id in chgs192.keys():
          self.assertTrue(mol.GetAtomWithIdx(id).GetFormalCharge() == chgs192[id])
    self.assertRaises(StopIteration, lambda: next(sdSup))
    sdSup.reset()

    ns = [mol.GetProp("_Name") for mol in sdSup]
    self.assertTrue(ns == molNames)

    sdSup = Chem.SDMolSupplier(fileN, 0)
    for mol in sdSup:
      self.assertTrue(not mol.HasProp("numArom"))

    sdSup = Chem.SDMolSupplier(fileN)
    self.assertTrue(len(sdSup) == 16)
    mol = sdSup[5]
    self.assertTrue(mol.GetProp("_Name") == "170")

    # test handling of H removal:
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'withHs.sdf')
    sdSup = Chem.SDMolSupplier(fileN)
    m = next(sdSup)
    self.assertTrue(m)
    self.assertTrue(m.GetNumAtoms() == 23)
    m = next(sdSup)
    self.assertTrue(m)
    self.assertTrue(m.GetNumAtoms() == 28)

    sdSup = Chem.SDMolSupplier(fileN, removeHs=False)
    m = next(sdSup)
    self.assertTrue(m)
    self.assertTrue(m.GetNumAtoms() == 39)
    m = next(sdSup)
    self.assertTrue(m)
    self.assertTrue(m.GetNumAtoms() == 30)

    with open(fileN, 'rb') as dFile:
      d = dFile.read()
    sdSup.SetData(d)
    m = next(sdSup)
    self.assertTrue(m)
    self.assertTrue(m.GetNumAtoms() == 23)
    m = next(sdSup)
    self.assertTrue(m)
    self.assertTrue(m.GetNumAtoms() == 28)

    sdSup.SetData(d, removeHs=False)
    m = next(sdSup)
    self.assertTrue(m)
    self.assertTrue(m.GetNumAtoms() == 39)
    m = next(sdSup)
    self.assertTrue(m)
    self.assertTrue(m.GetNumAtoms() == 30)

    # test strictParsing1:
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'strictLax1.sdf')
    #strict from file
    sdSup = Chem.SDMolSupplier(fileN, strictParsing=True)

    i = 0
    for mol in sdSup:
      self.assertTrue(mol.HasProp("_Name"))
      if (i == 0):
        self.assertTrue(not mol.HasProp("ID"))
      self.assertTrue(not mol.HasProp("ANOTHER_PROPERTY"))
      i += 1
    self.assertTrue(i == 2)

    #lax from file
    sdSup = Chem.SDMolSupplier(fileN, strictParsing=False)

    i = 0
    for mol in sdSup:
      self.assertTrue(mol.HasProp("_Name"))
      self.assertTrue(mol.HasProp("ID"))
      self.assertTrue(mol.HasProp("ANOTHER_PROPERTY"))
      i += 1
    self.assertTrue(i == 2)

    #strict from text
    with open(fileN, 'rb') as dFile:
      d = dFile.read()
    sdSup = Chem.SDMolSupplier()
    sdSup.SetData(d, strictParsing=True)

    i = 0
    for mol in sdSup:
      self.assertTrue(mol.HasProp("_Name"))
      if (i == 0):
        self.assertTrue(not mol.HasProp("ID"))
      self.assertTrue(not mol.HasProp("ANOTHER_PROPERTY"))
      i += 1
    self.assertTrue(i == 2)

    #lax from text
    sdSup = Chem.SDMolSupplier()
    sdSup.SetData(d, strictParsing=False)

    i = 0
    for mol in sdSup:
      self.assertTrue(mol.HasProp("_Name"))
      self.assertTrue(mol.HasProp("ID"))
      self.assertTrue(mol.HasProp("ANOTHER_PROPERTY"))
      i += 1
    self.assertTrue(i == 2)

    # test strictParsing2:
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'strictLax2.sdf')
    #strict from file
    sdSup = Chem.SDMolSupplier(fileN, strictParsing=True)

    i = 0
    for mol in sdSup:
      self.assertTrue(mol.HasProp("_Name"))
      self.assertTrue(mol.HasProp("ID"))
      self.assertTrue(mol.GetProp("ID") == "Lig1")
      self.assertTrue(mol.HasProp("ANOTHER_PROPERTY"))
      self.assertTrue(
        mol.GetProp("ANOTHER_PROPERTY") == "No blank line before dollars\n"
        "$$$$\n"
        "Structure1\n"
        "csChFnd70/05230312262D")
      i += 1
    self.assertTrue(i == 1)

    #lax from file
    sdSup = Chem.SDMolSupplier(fileN, strictParsing=False)

    i = 0
    for mol in sdSup:
      self.assertTrue(mol.HasProp("_Name"))
      self.assertTrue(mol.HasProp("ID"))
      self.assertTrue(mol.GetProp("ID") == "Lig2")
      self.assertTrue(mol.HasProp("ANOTHER_PROPERTY"))
      self.assertTrue(mol.GetProp("ANOTHER_PROPERTY") == "Value2")
      i += 1
    self.assertTrue(i == 1)

    #strict from text
    with open(fileN, 'rb') as dFile:
      d = dFile.read()
    sdSup = Chem.SDMolSupplier()
    sdSup.SetData(d, strictParsing=True)

    i = 0
    for mol in sdSup:
      self.assertTrue(mol.HasProp("_Name"))
      self.assertTrue(mol.HasProp("ID"))
      self.assertTrue(mol.GetProp("ID") == "Lig1")
      self.assertTrue(mol.HasProp("ANOTHER_PROPERTY"))
      self.assertTrue(
        mol.GetProp("ANOTHER_PROPERTY") == "No blank line before dollars\n"
        "$$$$\n"
        "Structure1\n"
        "csChFnd70/05230312262D")
      i += 1
    self.assertTrue(i == 1)

    #lax from text
    sdSup = Chem.SDMolSupplier()
    sdSup.SetData(d, strictParsing=False)

    i = 0
    for mol in sdSup:
      self.assertTrue(mol.HasProp("_Name"))
      self.assertTrue(mol.HasProp("ID"))
      self.assertTrue(mol.GetProp("ID") == "Lig2")
      self.assertTrue(mol.HasProp("ANOTHER_PROPERTY"))
      self.assertTrue(mol.GetProp("ANOTHER_PROPERTY") == "Value2")
      i += 1
    self.assertTrue(i == 1)

  def test26SmiMolSupplier(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'first_200.tpsa.csv')
    #fileN = "../FileParsers/test_data/first_200.tpsa.csv"
    smiSup = Chem.SmilesMolSupplier(fileN, ",", 0, -1)
    mol = smiSup[16]
    self.assertTrue(mol.GetProp("TPSA") == "46.25")

    mol = smiSup[8]
    self.assertTrue(mol.GetProp("TPSA") == "65.18")

    self.assertTrue(len(smiSup) == 200)

    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'fewSmi.csv')
    #fileN = "../FileParsers/test_data/fewSmi.csv"
    smiSup = Chem.SmilesMolSupplier(fileN, delimiter=",", smilesColumn=1, nameColumn=0, titleLine=0)
    names = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
    i = 0
    for mol in smiSup:
      self.assertTrue(mol.GetProp("_Name") == names[i])
      i += 1

    mol = smiSup[3]

    self.assertTrue(mol.GetProp("_Name") == "4")
    self.assertTrue(mol.GetProp("Column_2") == "82.78")

    # and test doing a supplier from text:
    with open(fileN, 'r') as inF:
      inD = inF.read()
    smiSup.SetData(inD, delimiter=",", smilesColumn=1, nameColumn=0, titleLine=0)
    names = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
    i = 0
    # iteration interface:
    for mol in smiSup:
      self.assertTrue(mol.GetProp("_Name") == names[i])
      i += 1
    self.assertTrue(i == 10)
    # random access:
    mol = smiSup[3]
    self.assertTrue(len(smiSup) == 10)
    self.assertTrue(mol.GetProp("_Name") == "4")
    self.assertTrue(mol.GetProp("Column_2") == "82.78")

    # issue 113:
    smiSup.SetData(inD, delimiter=",", smilesColumn=1, nameColumn=0, titleLine=0)
    self.assertTrue(len(smiSup) == 10)

    # and test failure handling:
    inD = """mol-1,CCC
mol-2,CCCC
mol-3,fail
mol-4,CCOC
    """
    smiSup.SetData(inD, delimiter=",", smilesColumn=1, nameColumn=0, titleLine=0)
    # there are 4 entries in the supplier:
    self.assertTrue(len(smiSup) == 4)
    # but the 3rd is a None:
    self.assertTrue(smiSup[2] is None)


    text="Id SMILES Column_2\n"+\
    "mol-1 C 1.0\n"+\
    "mol-2 CC 4.0\n"+\
    "mol-4 CCCC 16.0"
    smiSup.SetData(text, delimiter=" ", smilesColumn=1, nameColumn=0, titleLine=1)
    self.assertTrue(len(smiSup) == 3)
    self.assertTrue(smiSup[0])
    self.assertTrue(smiSup[1])
    self.assertTrue(smiSup[2])
    m = [x for x in smiSup]
    self.assertTrue(smiSup[2])
    self.assertTrue(len(m) == 3)
    self.assertTrue(m[0].GetProp("Column_2") == "1.0")

    # test simple parsing and Issue 114:
    smis = ['CC', 'CCC', 'CCOC', 'CCCOCC', 'CCCOCCC']
    inD = '\n'.join(smis)
    smiSup.SetData(inD, delimiter=",", smilesColumn=0, nameColumn=-1, titleLine=0)
    self.assertTrue(len(smiSup) == 5)
    m = [x for x in smiSup]
    self.assertTrue(smiSup[4])
    self.assertTrue(len(m) == 5)

    # order dependence:
    smiSup.SetData(inD, delimiter=",", smilesColumn=0, nameColumn=-1, titleLine=0)
    self.assertTrue(smiSup[4])
    self.assertTrue(len(smiSup) == 5)

    # this was a nasty BC:
    # asking for a particular entry with a higher index than what we've
    # already seen resulted in a duplicate:
    smis = ['CC', 'CCC', 'CCOC', 'CCCCOC']
    inD = '\n'.join(smis)
    smiSup.SetData(inD, delimiter=",", smilesColumn=0, nameColumn=-1, titleLine=0)
    m = next(smiSup)
    m = smiSup[3]
    self.assertTrue(len(smiSup) == 4)

    with self.assertRaisesRegex(Exception, ""):
      smiSup[4]

    smiSup.SetData(inD, delimiter=",", smilesColumn=0, nameColumn=-1, titleLine=0)
    with self.assertRaisesRegex(Exception, ""):
      smiSup[4]

    sys.stderr.write(
      '>>> This may result in an infinite loop.  It should finish almost instantly\n')
    self.assertEqual(len(smiSup), 4)
    sys.stderr.write('<<< OK, it finished.\n')

  def test27SmilesWriter(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'fewSmi.csv')
    #fileN = "../FileParsers/test_data/fewSmi.csv"

    smiSup = Chem.SmilesMolSupplier(fileN, delimiter=",", smilesColumn=1, nameColumn=0, titleLine=0)
    propNames = []
    propNames.append("Column_2")
    ofile = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Wrap', 'test_data',
                         'outSmiles.txt')
    writer = Chem.SmilesWriter(ofile)
    writer.SetProps(propNames)
    for mol in smiSup:
      writer.write(mol)
    writer.flush()

  def test28SmilesReverse(self):
    names = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
    props = [
      "34.14", "25.78", "106.51", "82.78", "60.16", "87.74", "37.38", "77.28", "65.18", "0.00"
    ]
    ofile = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Wrap', 'test_data',
                         'outSmiles.txt')
    #ofile = "test_data/outSmiles.csv"
    smiSup = Chem.SmilesMolSupplier(ofile)
    i = 0
    for mol in smiSup:
      #print([repr(x) for x in mol.GetPropNames()])
      self.assertTrue(mol.GetProp("_Name") == names[i])
      self.assertTrue(mol.GetProp("Column_2") == props[i])
      i += 1

  def writerSDFile(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'NCI_aids_few.sdf')
    #fileN = "../FileParsers/test_data/NCI_aids_few.sdf"
    ofile = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Wrap', 'test_data',
                         'outNCI_few.sdf')
    writer = Chem.SDWriter(ofile)
    sdSup = Chem.SDMolSupplier(fileN)
    for mol in sdSup:
      writer.write(mol)
    writer.flush()

  def test29SDWriterLoop(self):
    self.writerSDFile()
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Wrap', 'test_data',
                         'outNCI_few.sdf')
    sdSup = Chem.SDMolSupplier(fileN)
    molNames = [
      "48", "78", "128", "163", "164", "170", "180", "186", "192", "203", "210", "211", "213",
      "220", "229", "256"
    ]
    chgs192 = {8: 1, 11: 1, 15: -1, 18: -1, 20: 1, 21: 1, 23: -1, 25: -1}
    i = 0

    for mol in sdSup:
      #print('mol:',mol)
      #print('\t',molNames[i])
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      i += 1
      if (mol.GetProp("_Name") == "192"):
        # test parsed charges on one of the molecules
        for id in chgs192.keys():
          self.assertTrue(mol.GetAtomWithIdx(id).GetFormalCharge() == chgs192[id])

  def test30Issues109and110(self):
    """ issues 110 and 109 were both related to handling of explicit Hs in
       SMILES input.

    """
    m1 = Chem.MolFromSmiles('N12[CH](SC(C)(C)[CH]1C(O)=O)[CH](C2=O)NC(=O)[CH](N)c3ccccc3')
    self.assertTrue(m1.GetNumAtoms() == 24)
    m2 = Chem.MolFromSmiles(
      'C1C=C([CH](N)C(=O)N[C]2([H])[C]3([H])SC(C)(C)[CH](C(=O)O)N3C(=O)2)C=CC=1')
    self.assertTrue(m2.GetNumAtoms() == 24)

    smi1 = Chem.MolToSmiles(m1)
    smi2 = Chem.MolToSmiles(m2)
    self.assertTrue(smi1 == smi2)

    m1 = Chem.MolFromSmiles('[H]CCl')
    self.assertTrue(m1.GetNumAtoms() == 2)
    self.assertTrue(m1.GetAtomWithIdx(0).GetNumExplicitHs() == 1)
    m1 = Chem.MolFromSmiles('[H][CH2]Cl')
    self.assertTrue(m1.GetNumAtoms() == 2)
    self.assertTrue(m1.GetAtomWithIdx(0).GetNumExplicitHs() == 3)
    m2 = Chem.AddHs(m1)
    self.assertTrue(m2.GetNumAtoms() == 5)
    m2 = Chem.RemoveHs(m2)
    self.assertTrue(m2.GetNumAtoms() == 2)

  def test31ChiralitySmiles(self):
    m1 = Chem.MolFromSmiles('F[C@](Br)(I)Cl')
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 5)
    self.assertTrue(Chem.MolToSmiles(m1, 1) == 'F[C@](Cl)(Br)I', Chem.MolToSmiles(m1, 1))

    m1 = Chem.MolFromSmiles('CC1C[C@@]1(Cl)F')
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 6)
    self.assertTrue(Chem.MolToSmiles(m1, 1) == 'CC1C[C@]1(F)Cl', Chem.MolToSmiles(m1, 1))

    m1 = Chem.MolFromSmiles('CC1C[C@]1(Cl)F')
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 6)
    self.assertTrue(Chem.MolToSmiles(m1, 1) == 'CC1C[C@@]1(F)Cl', Chem.MolToSmiles(m1, 1))

  def test31aChiralitySubstructs(self):
    m1 = Chem.MolFromSmiles('CC1C[C@@]1(Cl)F')
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 6)
    self.assertTrue(Chem.MolToSmiles(m1, 1) == 'CC1C[C@]1(F)Cl', Chem.MolToSmiles(m1, 1))

    m2 = Chem.MolFromSmiles('CC1C[C@]1(Cl)F')
    self.assertTrue(m2 is not None)
    self.assertTrue(m2.GetNumAtoms() == 6)
    self.assertTrue(Chem.MolToSmiles(m2, 1) == 'CC1C[C@@]1(F)Cl', Chem.MolToSmiles(m2, 1))

    self.assertTrue(m1.HasSubstructMatch(m1))
    self.assertTrue(m1.HasSubstructMatch(m2))
    self.assertTrue(m1.HasSubstructMatch(m1, useChirality=True))
    self.assertTrue(not m1.HasSubstructMatch(m2, useChirality=True))

  def _test32MolFilesWithChirality(self):
    inD = """chiral1.mol
  ChemDraw10160313232D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0553    0.6188    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.0553   -0.2062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7697   -0.6188    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6592   -0.6188    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -0.7697   -0.2062    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  2  4  1  1
  2  5  1  0
M  END
"""
    m1 = Chem.MolFromMolBlock(inD)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 5)
    smi = Chem.MolToSmiles(m1)
    self.assertTrue(smi == 'F[C@](Cl)(Br)I', smi)

    inD = """chiral2.cdxml
  ChemDraw10160314052D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0553    0.6188    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.0553   -0.2062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7697   -0.6188    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6592   -0.6188    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -0.7697   -0.2062    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  2  4  1  6
  2  5  1  0
M  END
"""
    m1 = Chem.MolFromMolBlock(inD)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 5)
    self.assertTrue(Chem.MolToSmiles(m1, 1) == 'F[C@@](Cl)(Br)I')

    inD = """chiral1.mol
  ChemDraw10160313232D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0553    0.6188    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.0553   -0.2062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7697   -0.2062    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -0.6592   -0.6188    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
    0.7697   -0.6188    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  2  4  1  1
  2  5  1  0
M  END
"""
    m1 = Chem.MolFromMolBlock(inD)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 5)
    self.assertTrue(Chem.MolToSmiles(m1, 1) == 'F[C@](Cl)(Br)I')

    inD = """chiral1.mol
  ChemDraw10160313232D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0553   -0.2062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7697   -0.2062    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -0.6592   -0.6188    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
    0.7697   -0.6188    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0
    0.0553    0.6188    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  1
  1  4  1  0
  1  5  1  0
M  END
"""
    m1 = Chem.MolFromMolBlock(inD)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 5)
    self.assertTrue(Chem.MolToSmiles(m1, 1) == 'F[C@](Cl)(Br)I')

    inD = """chiral3.mol
  ChemDraw10160314362D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.4125    0.6188    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.4125   -0.2062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3020   -0.6188    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -0.4125   -0.2062    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  1
  2  4  1  0
M  END

"""
    m1 = Chem.MolFromMolBlock(inD)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 4)
    self.assertTrue(Chem.MolToSmiles(m1, 1) == 'F[C@H](Cl)Br')

    inD = """chiral4.mol
  ChemDraw10160314362D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.4125    0.6188    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.4125   -0.2062    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3020   -0.6188    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -0.4125   -0.2062    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  1
  2  4  1  0
M  END

"""
    m1 = Chem.MolFromMolBlock(inD)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 4)
    self.assertTrue(Chem.MolToSmiles(m1, 1) == 'FN(Cl)Br')

    inD = """chiral5.mol
  ChemDraw10160314362D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.4125    0.6188    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.4125   -0.2062    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3020   -0.6188    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -0.4125   -0.2062    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  1
  2  4  1  0
M  CHG  1   2   1
M  END

"""
    m1 = Chem.MolFromMolBlock(inD)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 4)
    self.assertTrue(Chem.MolToSmiles(m1, 1) == 'F[N@H+](Cl)Br')

    inD = """Case 10-14-3
  ChemDraw10140308512D

  4  3  0  0  0  0  0  0  0  0999 V2000
   -0.8250   -0.4125    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8250   -0.4125    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.4125    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  2  4  1  1
M  END
"""
    m1 = Chem.MolFromMolBlock(inD)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 4)
    self.assertTrue(Chem.MolToSmiles(m1, 1) == 'F[C@H](Cl)Br')

    inD = """Case 10-14-4
  ChemDraw10140308512D

  4  3  0  0  0  0  0  0  0  0999 V2000
   -0.8250   -0.4125    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8250   -0.4125    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.4125    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  1
  2  4  1  0
M  END

"""
    m1 = Chem.MolFromMolBlock(inD)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 4)
    self.assertTrue(Chem.MolToSmiles(m1, 1) == 'F[C@H](Cl)Br')

    inD = """chiral4.mol
  ChemDraw10160315172D

  6  6  0  0  0  0  0  0  0  0999 V2000
   -0.4422    0.1402    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4422   -0.6848    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2723   -0.2723    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8547    0.8547    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6848    0.4422    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.8547   -0.8547    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  1  1  0
  1  4  1  0
  3  5  1  1
  3  6  1  0
M  END
"""
    m1 = Chem.MolFromMolBlock(inD)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 6)
    self.assertTrue(Chem.MolToSmiles(m1, 1) == 'CC1C[C@@]1(F)Cl', Chem.MolToSmiles(m1, 1))

    inD = """chiral4.mol
  ChemDraw10160315172D

  6  6  0  0  0  0  0  0  0  0999 V2000
   -0.4422    0.1402    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4422   -0.6848    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2723   -0.2723    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8547    0.8547    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6848    0.4422    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.8547   -0.8547    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  1  1  0
  1  4  1  0
  3  5  1  6
  3  6  1  0
M  END
"""
    m1 = Chem.MolFromMolBlock(inD)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 6)
    self.assertTrue(Chem.MolToSmiles(m1, 1) == 'CC1C[C@]1(F)Cl', Chem.MolToSmiles(m1, 1))

  def test33Issue65(self):
    """ issue 65 relates to handling of [H] in SMARTS

    """
    m1 = Chem.MolFromSmiles('OC(O)(O)O')
    m2 = Chem.MolFromSmiles('OC(O)O')
    m3 = Chem.MolFromSmiles('OCO')
    q1 = Chem.MolFromSmarts('OC[H]', 1)
    q2 = Chem.MolFromSmarts('O[C;H1]', 1)
    q3 = Chem.MolFromSmarts('O[C;H1][H]', 1)

    self.assertTrue(not m1.HasSubstructMatch(q1))
    self.assertTrue(not m1.HasSubstructMatch(q2))
    self.assertTrue(not m1.HasSubstructMatch(q3))

    self.assertTrue(m2.HasSubstructMatch(q1))
    self.assertTrue(m2.HasSubstructMatch(q2))
    self.assertTrue(m2.HasSubstructMatch(q3))

    self.assertTrue(m3.HasSubstructMatch(q1))
    self.assertTrue(not m3.HasSubstructMatch(q2))
    self.assertTrue(not m3.HasSubstructMatch(q3))

    m1H = Chem.AddHs(m1)
    m2H = Chem.AddHs(m2)
    m3H = Chem.AddHs(m3)
    q1 = Chem.MolFromSmarts('OC[H]')
    q2 = Chem.MolFromSmarts('O[C;H1]')
    q3 = Chem.MolFromSmarts('O[C;H1][H]')

    self.assertTrue(not m1H.HasSubstructMatch(q1))
    self.assertTrue(not m1H.HasSubstructMatch(q2))
    self.assertTrue(not m1H.HasSubstructMatch(q3))

    #m2H.Debug()
    self.assertTrue(m2H.HasSubstructMatch(q1))
    self.assertTrue(m2H.HasSubstructMatch(q2))
    self.assertTrue(m2H.HasSubstructMatch(q3))

    self.assertTrue(m3H.HasSubstructMatch(q1))
    self.assertTrue(not m3H.HasSubstructMatch(q2))
    self.assertTrue(not m3H.HasSubstructMatch(q3))

  def test34Issue124(self):
    """ issue 124 relates to calculation of the distance matrix

    """
    m = Chem.MolFromSmiles('CC=C')
    d = Chem.GetDistanceMatrix(m, 0)
    self.assertTrue(feq(d[0, 1], 1.0))
    self.assertTrue(feq(d[0, 2], 2.0))
    # force an update:
    d = Chem.GetDistanceMatrix(m, 1, 0, 1)
    self.assertTrue(feq(d[0, 1], 1.0))
    self.assertTrue(feq(d[0, 2], 1.5))

  def test35ChiralityPerception(self):
    """ Test perception of chirality and CIP encoding
    """
    m = Chem.MolFromSmiles('F[C@]([C@])(Cl)Br')
    Chem.AssignStereochemistry(m, 1)
    self.assertTrue(m.GetAtomWithIdx(1).HasProp('_CIPCode'))
    self.assertFalse(m.GetAtomWithIdx(2).HasProp('_CIPCode'))
    Chem.RemoveStereochemistry(m)
    self.assertFalse(m.GetAtomWithIdx(1).HasProp('_CIPCode'))

    m = Chem.MolFromSmiles('F[C@H](C)C')
    Chem.AssignStereochemistry(m, 1)
    self.assertTrue(m.GetAtomWithIdx(1).GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED)
    self.assertFalse(m.GetAtomWithIdx(1).HasProp('_CIPCode'))

    m = Chem.MolFromSmiles('F\\C=C/Cl')
    self.assertTrue(m.GetBondWithIdx(0).GetStereo() == Chem.BondStereo.STEREONONE)
    self.assertTrue(m.GetBondWithIdx(1).GetStereo() == Chem.BondStereo.STEREOZ)
    atoms = m.GetBondWithIdx(1).GetStereoAtoms()
    self.assertTrue(0 in atoms)
    self.assertTrue(3 in atoms)
    self.assertTrue(m.GetBondWithIdx(2).GetStereo() == Chem.BondStereo.STEREONONE)
    Chem.RemoveStereochemistry(m)
    self.assertTrue(m.GetBondWithIdx(1).GetStereo() == Chem.BondStereo.STEREONONE)

    m = Chem.MolFromSmiles('F\\C=CCl')
    self.assertTrue(m.GetBondWithIdx(1).GetStereo() == Chem.BondStereo.STEREONONE)

  def checkDefaultBondProperties(self, m):
    for bond in m.GetBonds():
      self.assertIn(bond.GetBondType(), [Chem.BondType.SINGLE, Chem.BondType.DOUBLE])
      self.assertEqual(bond.GetBondDir(), Chem.BondDir.NONE)
      self.assertEqual(list(bond.GetStereoAtoms()), [])
      self.assertEqual(bond.GetStereo(), Chem.BondStereo.STEREONONE)

  def assertHasDoubleBondStereo(self, smi):
    m = Chem.MolFromSmiles(smi)

    self.checkDefaultBondProperties(m)

    Chem.FindPotentialStereoBonds(m)

    for bond in m.GetBonds():
      self.assertIn(bond.GetBondType(), [Chem.BondType.SINGLE, Chem.BondType.DOUBLE])
      self.assertEqual(bond.GetBondDir(), Chem.BondDir.NONE)

      if bond.GetBondType() == Chem.BondType.DOUBLE:
        self.assertEqual(bond.GetStereo(), Chem.BondStereo.STEREOANY)
        self.assertEqual(len(list(bond.GetStereoAtoms())), 2)
      else:
        self.assertEqual(list(bond.GetStereoAtoms()), [])
        self.assertEqual(bond.GetStereo(), Chem.BondStereo.STEREONONE)

  def testFindPotentialStereoBonds(self):
    self.assertHasDoubleBondStereo("FC=CF")
    self.assertHasDoubleBondStereo("FC(Cl)=C(Br)I")
    self.assertHasDoubleBondStereo("FC=CC=CC=CCl")
    self.assertHasDoubleBondStereo("C1CCCCC1C=CC1CCCCC1")

  def assertDoesNotHaveDoubleBondStereo(self, smi):
    m = Chem.MolFromSmiles(smi)
    self.checkDefaultBondProperties(m)
    Chem.FindPotentialStereoBonds(m)
    self.checkDefaultBondProperties(m)

  def testFindPotentialStereoBondsShouldNotFindThisDoubleBondAsStereo(self):
    self.assertDoesNotHaveDoubleBondStereo("FC(F)=CF")
    self.assertDoesNotHaveDoubleBondStereo("C=C")
    self.assertDoesNotHaveDoubleBondStereo("C1CCCCC1C(C1CCCCC1)=CC1CCCCC1")

  def assertDoubleBondStereo(self, smi, stereo):
    mol = Chem.MolFromSmiles(smi)

    bond = mol.GetBondWithIdx(1)
    self.assertEqual(bond.GetBondType(), Chem.BondType.DOUBLE)
    self.assertEqual(bond.GetStereo(), stereo)
    self.assertEqual(list(bond.GetStereoAtoms()), [0, 3])

  def allStereoBonds(self, bonds):
    for bond in bonds:
      self.assertEqual(len(list(bond.GetStereoAtoms())), 2)

  def testBondSetStereo(self):
    for testAssignStereo in [False, True]:
      mol = Chem.MolFromSmiles("FC=CF")
      Chem.FindPotentialStereoBonds(mol)

      for bond in mol.GetBonds():
        if (bond.GetBondType() == Chem.BondType.DOUBLE
            and bond.GetStereo() == Chem.BondStereo.STEREOANY):
          break
      self.assertEqual(bond.GetBondType(), Chem.BondType.DOUBLE)
      self.assertEqual(bond.GetStereo(), Chem.BondStereo.STEREOANY)
      self.assertEqual(list(bond.GetStereoAtoms()), [0, 3])

      bond.SetStereo(Chem.BondStereo.STEREOTRANS)
      self.assertEqual(bond.GetStereo(), Chem.BondStereo.STEREOTRANS)
      if testAssignStereo:  # should be invariant of Chem.AssignStereochemistry being called
        Chem.AssignStereochemistry(mol, force=True)
      smi = Chem.MolToSmiles(mol, isomericSmiles=True)
      self.allStereoBonds([bond])
      self.assertEqual(smi, "F/C=C/F")
      self.assertDoubleBondStereo(smi, Chem.BondStereo.STEREOE)

      bond.SetStereo(Chem.BondStereo.STEREOCIS)
      self.assertEqual(bond.GetStereo(), Chem.BondStereo.STEREOCIS)
      if testAssignStereo:
        Chem.AssignStereochemistry(mol, force=True)
      smi = Chem.MolToSmiles(mol, isomericSmiles=True)
      self.allStereoBonds([bond])
      self.assertEqual(smi, r"F/C=C\F")
      self.assertDoubleBondStereo(smi, Chem.BondStereo.STEREOZ)

  def recursive_enumerate_stereo_bonds(self, mol, done_bonds, bonds):
    if not bonds:
      yield done_bonds, Chem.Mol(mol)
      return

    bond = bonds[0]
    child_bonds = bonds[1:]
    self.assertEqual(len(list(bond.GetStereoAtoms())), 2)
    bond.SetStereo(Chem.BondStereo.STEREOTRANS)
    for isomer in self.recursive_enumerate_stereo_bonds(mol, done_bonds + [Chem.BondStereo.STEREOE],
                                                        child_bonds):
      yield isomer

    self.assertEqual(len(list(bond.GetStereoAtoms())), 2)
    bond.SetStereo(Chem.BondStereo.STEREOCIS)
    for isomer in self.recursive_enumerate_stereo_bonds(mol, done_bonds + [Chem.BondStereo.STEREOZ],
                                                        child_bonds):
      yield isomer

  def testBondSetStereoDifficultCase(self):
    unspec_smiles = "CCC=CC(CO)=C(C)CC"
    mol = Chem.MolFromSmiles(unspec_smiles)
    Chem.FindPotentialStereoBonds(mol)

    stereo_bonds = []
    for bond in mol.GetBonds():
      if bond.GetStereo() == Chem.BondStereo.STEREOANY:
        stereo_bonds.append(bond)

    isomers = set()
    for bond_stereo, isomer in self.recursive_enumerate_stereo_bonds(mol, [], stereo_bonds):
      self.allStereoBonds(stereo_bonds)
      isosmi = Chem.MolToSmiles(isomer, isomericSmiles=True)
      self.allStereoBonds(stereo_bonds)

      self.assertNotIn(isosmi, isomers)
      isomers.add(isosmi)

      isomol = Chem.MolFromSmiles(isosmi)
      round_trip_stereo = [
        b.GetStereo() for b in isomol.GetBonds() if b.GetStereo() != Chem.BondStereo.STEREONONE
      ]

      self.assertEqual(bond_stereo, round_trip_stereo)

    self.assertEqual(len(isomers), 4)

  def getNumUnspecifiedBondStereo(self, smi):
    mol = Chem.MolFromSmiles(smi)
    Chem.FindPotentialStereoBonds(mol)

    count = 0
    for bond in mol.GetBonds():
      if bond.GetStereo() == Chem.BondStereo.STEREOANY:
        count += 1

    return count

  def testBondSetStereoReallyDifficultCase(self):
    # this one is much trickier because a double bond can gain and
    # lose it's stereochemistry based upon whether 2 other double
    # bonds have the same or different stereo chemistry.

    unspec_smiles = "CCC=CC(C=CCC)=C(CO)CC"
    mol = Chem.MolFromSmiles(unspec_smiles)
    Chem.FindPotentialStereoBonds(mol)

    stereo_bonds = []
    for bond in mol.GetBonds():
      if bond.GetStereo() == Chem.BondStereo.STEREOANY:
        stereo_bonds.append(bond)

    self.assertEqual(len(stereo_bonds), 2)

    isomers = set()
    for bond_stereo, isomer in self.recursive_enumerate_stereo_bonds(mol, [], stereo_bonds):
      isosmi = Chem.MolToSmiles(isomer, isomericSmiles=True)
      isomers.add(isosmi)

    self.assertEqual(len(isomers), 3)

    # one of these then gains a new stereo bond due to the
    # introduction of a new symmetry
    counts = {}
    for isosmi in isomers:
      num_unspecified = self.getNumUnspecifiedBondStereo(isosmi)
      counts[num_unspecified] = counts.get(num_unspecified, 0) + 1

    # 2 of the isomers don't have any unspecified bond stereo centers
    # left, 1 does
    self.assertEqual(counts, {0: 2, 1: 1})

  def assertBondSetStereoIsAlwaysEquivalent(self, all_smiles, desired_stereo, bond_idx):
    refSmiles = None
    for smi in all_smiles:
      mol = Chem.MolFromSmiles(smi)

      doubleBond = None
      for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
          doubleBond = bond

      self.assertTrue(doubleBond is not None)

      Chem.FindPotentialStereoBonds(mol)
      doubleBond.SetStereo(desired_stereo)

      isosmi = Chem.MolToSmiles(mol, isomericSmiles=True)

      if refSmiles is None:
        refSmiles = isosmi

      self.assertEqual(refSmiles, isosmi)

  def testBondSetStereoAllHalogens(self):
    # can't get much more brutal than this test
    from itertools import combinations, permutations
    halogens = ['F', 'Cl', 'Br', 'I']

    # binary double bond stereo
    for unique_set in combinations(halogens, 2):
      all_smiles = []
      for fmt in ['%sC=C%s', 'C(%s)=C%s']:
        for ordering in permutations(unique_set):
          all_smiles.append(fmt % ordering)

      #print(fmt, all_smiles)
      for desired_stereo in [Chem.BondStereo.STEREOTRANS, Chem.BondStereo.STEREOCIS]:
        self.assertBondSetStereoIsAlwaysEquivalent(all_smiles, desired_stereo, 1)

    # tertiary double bond stereo
    for unique_set in combinations(halogens, 3):
      for mono_side in unique_set:
        halogens_left = list(unique_set)
        halogens_left.remove(mono_side)
        for binary_side in combinations(halogens_left, 2):
          all_smiles = []

          for binary_side_permutation in permutations(binary_side):
            all_smiles.append('%sC=C(%s)%s' % ((mono_side, ) + binary_side_permutation))
            all_smiles.append('C(%s)=C(%s)%s' % ((mono_side, ) + binary_side_permutation))

            all_smiles.append('%sC(%s)=C%s' % (binary_side_permutation + (mono_side, )))
            all_smiles.append('C(%s)(%s)=C%s' % (binary_side_permutation + (mono_side, )))

          #print(all_smiles)
          for desired_stereo in [Chem.BondStereo.STEREOTRANS, Chem.BondStereo.STEREOCIS]:
            self.assertBondSetStereoIsAlwaysEquivalent(all_smiles, desired_stereo, 1)

    # quaternary double bond stereo
    for unique_ordering in permutations(halogens):
      left_side = unique_ordering[:2]
      rght_side = unique_ordering[2:]

      all_smiles = []
      for left_side_permutation in permutations(left_side):
        for rght_side_permutation in permutations(rght_side):
          for smifmt in ['%sC(%s)=C(%s)%s', 'C(%s)(%s)=C(%s)%s']:
            all_smiles.append(smifmt % (left_side_permutation + rght_side_permutation))

      #print(all_smiles)
      for desired_stereo in [Chem.BondStereo.STEREOTRANS, Chem.BondStereo.STEREOCIS]:
        self.assertBondSetStereoIsAlwaysEquivalent(all_smiles, desired_stereo, 1)

  def testBondSetStereoAtoms(self):
    # use this difficult molecule that only generates 4 isomers, but
    # assume all double bonds are stereo!
    unspec_smiles = "CCC=CC(C=CCC)=C(CO)CC"
    mol = Chem.MolFromSmiles(unspec_smiles)

    def getNbr(atom, exclude):
      for nbr in atom.GetNeighbors():
        if nbr.GetIdx() not in exclude:
          return nbr
      raise ValueError("No neighbor found!")

    double_bonds = []
    for bond in mol.GetBonds():
      if bond.GetBondType() == 2:
        double_bonds.append(bond)

        exclude = {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}
        bgnNbr = getNbr(bond.GetBeginAtom(), exclude)
        endNbr = getNbr(bond.GetEndAtom(), exclude)

        bond.SetStereoAtoms(bgnNbr.GetIdx(), endNbr.GetIdx())

    self.assertEqual(len(double_bonds), 3)

    import itertools
    stereos = [Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ]
    isomers = set()
    for stereo_config in itertools.product(stereos, repeat=len(double_bonds)):
      for bond, stereo in zip(double_bonds, stereo_config):
        bond.SetStereo(stereo)
      smi = Chem.MolToSmiles(mol, True)
      isomers.add(smi)

    # the dependent double bond stereo isn't picked up by this, should it?
    self.assertEqual(len(isomers), 6)

    # round tripping them through one more time does pick up the dependency, so meh?
    round_trip_isomers = set()
    for smi in isomers:
      isosmi = Chem.MolToSmiles(Chem.MolFromSmiles(smi), True)
      round_trip_isomers.add(isosmi)

    self.assertEqual(len(round_trip_isomers), 4)

  def test36SubstructMatchStr(self):
    """ test the _SubstructMatchStr function """
    query = Chem.MolFromSmarts('[n,p]1ccccc1')
    self.assertTrue(query)
    mol = Chem.MolFromSmiles('N1=CC=CC=C1')
    self.assertTrue(mol.HasSubstructMatch(query))
    self.assertTrue(Chem._HasSubstructMatchStr(mol.ToBinary(), query))
    mol = Chem.MolFromSmiles('S1=CC=CC=C1')
    self.assertTrue(not Chem._HasSubstructMatchStr(mol.ToBinary(), query))
    self.assertTrue(not mol.HasSubstructMatch(query))
    mol = Chem.MolFromSmiles('P1=CC=CC=C1')
    self.assertTrue(mol.HasSubstructMatch(query))
    self.assertTrue(Chem._HasSubstructMatchStr(mol.ToBinary(), query))

  def test37SanitException(self):
    mol = Chem.MolFromSmiles('CC(C)(C)(C)C', 0)
    self.assertTrue(mol)
    self.assertRaises(ValueError, lambda: Chem.SanitizeMol(mol))

  def test38TDTSuppliers(self):
    data = """$SMI<Cc1nnc(N)nc1C>
CAS<17584-12-2>
|
$SMI<Cc1n[nH]c(=O)nc1N>
CAS<~>
|
$SMI<Cc1n[nH]c(=O)[nH]c1=O>
CAS<932-53-6>
|
$SMI<Cc1nnc(NN)nc1O>
CAS<~>
|"""
    suppl = Chem.TDTMolSupplier()
    suppl.SetData(data, "CAS")
    i = 0
    for mol in suppl:
      self.assertTrue(mol)
      self.assertTrue(mol.GetNumAtoms())
      self.assertTrue(mol.HasProp("CAS"))
      self.assertTrue(mol.HasProp("_Name"))
      self.assertTrue(mol.GetProp("CAS") == mol.GetProp("_Name"))
      self.assertTrue(mol.GetNumConformers() == 0)
      i += 1
    self.assertTrue(i == 4)
    self.assertTrue(len(suppl) == 4)

  def test38Issue266(self):
    """ test issue 266: generation of kekulized smiles"""
    mol = Chem.MolFromSmiles('c1ccccc1')
    Chem.Kekulize(mol)
    smi = Chem.MolToSmiles(mol)
    self.assertTrue(smi == 'c1ccccc1')
    smi = Chem.MolToSmiles(mol, kekuleSmiles=True)
    self.assertTrue(smi == 'C1=CC=CC=C1')

  def test39Issue273(self):
    """ test issue 273: MolFileComments and MolFileInfo props ending up in SD files

    """
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Wrap', 'test_data',
                         'outNCI_few.sdf')
    suppl = Chem.SDMolSupplier(fileN)
    ms = [x for x in suppl]
    for m in ms:
      self.assertTrue(m.HasProp('_MolFileInfo'))
      self.assertTrue(m.HasProp('_MolFileComments'))
    fName = tempfile.NamedTemporaryFile(suffix='.sdf', delete=False).name
    w = Chem.SDWriter(fName)
    w.SetProps(ms[0].GetPropNames())
    for m in ms:
      w.write(m)
    w = None

    with open(fName, 'r') as txtFile:
      txt = txtFile.read()
    os.unlink(fName)
    self.assertTrue(txt.find('MolFileInfo') == -1)
    self.assertTrue(txt.find('MolFileComments') == -1)

  def test40SmilesRootedAtAtom(self):
    """ test the rootAtAtom functionality

    """
    smi = 'CN(C)C'
    m = Chem.MolFromSmiles(smi)

    self.assertTrue(Chem.MolToSmiles(m) == 'CN(C)C')
    self.assertTrue(Chem.MolToSmiles(m, rootedAtAtom=1) == 'N(C)(C)C')

  def test41SetStreamIndices(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'NCI_aids_few.sdf')
    allIndices = []
    ifs = open(fileN, 'rb')
    addIndex = True
    line = True
    pos = 0
    while (line):
      if (addIndex):
        pos = ifs.tell()
      line = ifs.readline().decode('utf-8')
      if (line):
        if (addIndex):
          allIndices.append(pos)
        addIndex = (line[:4] == '$$$$')
    ifs.close()
    indices = allIndices
    sdSup = Chem.SDMolSupplier(fileN)
    molNames = [
      "48", "78", "128", "163", "164", "170", "180", "186", "192", "203", "210", "211", "213",
      "220", "229", "256"
    ]

    sdSup._SetStreamIndices(indices)
    self.assertTrue(len(sdSup) == 16)
    mol = sdSup[5]
    self.assertTrue(mol.GetProp("_Name") == "170")

    i = 0
    for mol in sdSup:
      self.assertTrue(mol)
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      i += 1

    ns = [mol.GetProp("_Name") for mol in sdSup]
    self.assertTrue(ns == molNames)

    # this can also be used to skip molecules in the file:
    indices = [allIndices[0], allIndices[2], allIndices[5]]
    sdSup._SetStreamIndices(indices)
    self.assertTrue(len(sdSup) == 3)
    mol = sdSup[2]
    self.assertTrue(mol.GetProp("_Name") == "170")

    # or to reorder them:
    indices = [allIndices[0], allIndices[5], allIndices[2]]
    sdSup._SetStreamIndices(indices)
    self.assertTrue(len(sdSup) == 3)
    mol = sdSup[1]
    self.assertTrue(mol.GetProp("_Name") == "170")

  def test42LifeTheUniverseAndEverything(self):
    self.assertTrue(True)

  def test43TplFileParsing(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'cmpd2.tpl')
    m1 = Chem.MolFromTPLFile(fileN)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 12)
    self.assertTrue(m1.GetNumConformers() == 2)

    m1 = Chem.MolFromTPLFile(fileN, skipFirstConf=True)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 12)
    self.assertTrue(m1.GetNumConformers() == 1)

    with open(fileN, 'r') as blockFile:
      block = blockFile.read()
    m1 = Chem.MolFromTPLBlock(block)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 12)
    self.assertTrue(m1.GetNumConformers() == 2)

  def test44TplFileWriting(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'cmpd2.tpl')
    m1 = Chem.MolFromTPLFile(fileN)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 12)
    self.assertTrue(m1.GetNumConformers() == 2)

    block = Chem.MolToTPLBlock(m1)
    m1 = Chem.MolFromTPLBlock(block)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 12)
    self.assertTrue(m1.GetNumConformers() == 2)

  def test45RingInfo(self):
    """ test the RingInfo class

    """
    smi = 'CNC'
    m = Chem.MolFromSmiles(smi)
    ri = m.GetRingInfo()
    self.assertTrue(ri)
    self.assertTrue(ri.NumRings() == 0)
    self.assertFalse(ri.IsAtomInRingOfSize(0, 3))
    self.assertFalse(ri.IsAtomInRingOfSize(1, 3))
    self.assertFalse(ri.IsAtomInRingOfSize(2, 3))
    self.assertFalse(ri.IsBondInRingOfSize(1, 3))
    self.assertFalse(ri.IsBondInRingOfSize(2, 3))
    if hasattr(Chem, 'FindRingFamilies'):
      self.assertEqual(ri.AtomRingFamilies(), ())
    if hasattr(Chem, 'FindRingFamilies'):
      self.assertEqual(ri.BondRingFamilies(), ())

    smi = 'C1CC2C1C2'
    m = Chem.MolFromSmiles(smi)
    ri = m.GetRingInfo()
    self.assertTrue(ri)
    self.assertEqual(ri.NumRings(), 2)
    self.assertFalse(ri.IsAtomInRingOfSize(0, 3))
    self.assertTrue(ri.IsAtomInRingOfSize(0, 4))
    self.assertFalse(ri.IsBondInRingOfSize(0, 3))
    self.assertTrue(ri.IsBondInRingOfSize(0, 4))
    self.assertTrue(ri.IsAtomInRingOfSize(2, 4))
    self.assertTrue(ri.IsAtomInRingOfSize(2, 3))
    self.assertTrue(ri.IsBondInRingOfSize(2, 3))
    self.assertTrue(ri.IsBondInRingOfSize(2, 4))
    self.assertEqual(ri.AtomRings(), ((0, 3, 2, 1), (4, 3, 2)))
    self.assertEqual(ri.BondRings(), ((4, 2, 1, 0), (3, 2, 5)))
    self.assertEqual(len(ri.AtomMembers(2)), 2)
    self.assertEqual(ri.AtomRingSizes(2), (4, 3))
    self.assertEqual(ri.AtomRingSizes(99), ())
    self.assertTrue(ri.AreAtomsInSameRing(2, 3))
    self.assertFalse(ri.AreAtomsInSameRing(1, 4))
    self.assertTrue(ri.AreAtomsInSameRingOfSize(2, 3, 3))
    self.assertTrue(ri.AreAtomsInSameRingOfSize(2, 3, 4))
    self.assertFalse(ri.AreAtomsInSameRingOfSize(2, 3, 5))
    self.assertEqual(len(ri.BondMembers(2)), 2)
    self.assertEqual(len(ri.BondMembers(0)), 1)
    self.assertTrue(ri.IsRingFused(0))
    self.assertTrue(ri.IsRingFused(1))
    self.assertTrue(ri.AreRingsFused(0, 1))
    self.assertTrue(ri.NumFusedBonds(0) == 1)
    self.assertTrue(ri.NumFusedBonds(1) == 1)
    self.assertEqual(ri.BondRingSizes(2), (4, 3))
    self.assertEqual(ri.BondRingSizes(0), (4, ))
    self.assertEqual(ri.BondRingSizes(99), ())
    self.assertTrue(ri.AreBondsInSameRing(1, 2))
    self.assertTrue(ri.AreBondsInSameRing(2, 5))
    self.assertFalse(ri.AreBondsInSameRing(1, 3))
    self.assertTrue(ri.AreBondsInSameRingOfSize(1, 2, 4))
    self.assertTrue(ri.AreBondsInSameRingOfSize(2, 5, 3))
    self.assertFalse(ri.AreBondsInSameRingOfSize(1, 2, 3))
    self.assertFalse(ri.AreBondsInSameRingOfSize(1, 3, 4))

    if hasattr(Chem, 'FindRingFamilies'):
      ri = m.GetRingInfo()
      self.assertFalse(ri.AreRingFamiliesInitialized())
      Chem.FindRingFamilies(m)
      ri = m.GetRingInfo()
      self.assertTrue(ri.AreRingFamiliesInitialized())
      self.assertEqual(ri.NumRingFamilies(), 2)
      self.assertEqual(sorted(ri.AtomRingFamilies()), [(0, 1, 2, 3), (2, 3, 4)])
      self.assertEqual(sorted(ri.BondRingFamilies()), [(0, 1, 2, 4), (2, 3, 5)])

  def test46ReplaceCore(self):
    """ test the ReplaceCore functionality

    """

    core = Chem.MolFromSmiles('C=O')

    smi = 'CCC=O'
    m = Chem.MolFromSmiles(smi)
    r = Chem.ReplaceCore(m, core)
    self.assertTrue(r)
    self.assertEqual(Chem.MolToSmiles(r, True), '[1*]CC')

    smi = 'C1CC(=O)CC1'
    m = Chem.MolFromSmiles(smi)
    r = Chem.ReplaceCore(m, core)
    self.assertTrue(r)
    self.assertEqual(Chem.MolToSmiles(r, True), '[1*]CCCC[2*]')

    smi = 'C1CC(=N)CC1'
    m = Chem.MolFromSmiles(smi)
    r = Chem.ReplaceCore(m, core)
    self.assertFalse(r)

    # smiles, smarts, replaceDummies, labelByIndex, useChirality
    expected = {
      ('C1O[C@@]1(OC)NC', 'C1O[C@]1(*)*', False, False, False): '[1*]OC.[2*]NC',
      ('C1O[C@@]1(OC)NC', 'C1O[C@]1(*)*', False, False, True): '[1*]NC.[2*]OC',
      ('C1O[C@@]1(OC)NC', 'C1O[C@]1(*)*', False, True, False): '[3*]OC.[4*]NC',
      ('C1O[C@@]1(OC)NC', 'C1O[C@]1(*)*', False, True, True): '[3*]NC.[4*]OC',
      ('C1O[C@@]1(OC)NC', 'C1O[C@]1(*)*', True, False, False): '[1*]C.[2*]C',
      ('C1O[C@@]1(OC)NC', 'C1O[C@]1(*)*', True, False, True): '[1*]C.[2*]C',
      ('C1O[C@@]1(OC)NC', 'C1O[C@]1(*)*', True, True, False): '[3*]C.[4*]C',
      ('C1O[C@@]1(OC)NC', 'C1O[C@]1(*)*', True, True, True): '[3*]C.[4*]C',
      ('C1O[C@]1(OC)NC', 'C1O[C@]1(*)*', False, False, False): '[1*]OC.[2*]NC',
      ('C1O[C@]1(OC)NC', 'C1O[C@]1(*)*', False, False, True): '[1*]OC.[2*]NC',
      ('C1O[C@]1(OC)NC', 'C1O[C@]1(*)*', False, True, False): '[3*]OC.[4*]NC',
      ('C1O[C@]1(OC)NC', 'C1O[C@]1(*)*', False, True, True): '[3*]OC.[4*]NC',
      ('C1O[C@]1(OC)NC', 'C1O[C@]1(*)*', True, False, False): '[1*]C.[2*]C',
      ('C1O[C@]1(OC)NC', 'C1O[C@]1(*)*', True, False, True): '[1*]C.[2*]C',
      ('C1O[C@]1(OC)NC', 'C1O[C@]1(*)*', True, True, False): '[3*]C.[4*]C',
      ('C1O[C@]1(OC)NC', 'C1O[C@]1(*)*', True, True, True): '[3*]C.[4*]C',
    }

    for (smiles, smarts, replaceDummies, labelByIndex,
         useChirality), expected_smiles in expected.items():
      mol = Chem.MolFromSmiles(smiles)
      core = Chem.MolFromSmarts(smarts)
      nm = Chem.ReplaceCore(mol, core, replaceDummies=replaceDummies, labelByIndex=labelByIndex,
                            useChirality=useChirality)

      if Chem.MolToSmiles(nm, True) != expected_smiles:
        print(
          "ReplaceCore(%r, %r, replaceDummies=%r, labelByIndex=%r, useChirality=%r" %
          (smiles, smarts, replaceDummies, labelByIndex, useChirality), file=sys.stderr)
        print("expected: %s\ngot: %s" % (expected_smiles, Chem.MolToSmiles(nm, True)),
              file=sys.stderr)
        self.assertEqual(expected_smiles, Chem.MolToSmiles(nm, True))

      matchVect = mol.GetSubstructMatch(core, useChirality=useChirality)
      nm = Chem.ReplaceCore(mol, core, matchVect, replaceDummies=replaceDummies,
                            labelByIndex=labelByIndex)
      if Chem.MolToSmiles(nm, True) != expected_smiles:
        print(
          "ReplaceCore(%r, %r, %r, replaceDummies=%r, labelByIndex=%rr" %
          (smiles, smarts, matchVect, replaceDummies, labelByIndex), file=sys.stderr)
        print("expected: %s\ngot: %s" % (expected_smiles, Chem.MolToSmiles(nm, True)),
              file=sys.stderr)
        self.assertEqual(expected_smiles, Chem.MolToSmiles(nm, True))

    mol = Chem.MolFromSmiles("C")
    smarts = Chem.MolFromSmarts("C")
    try:
      Chem.ReplaceCore(mol, smarts, (3, ))
      self.asssertFalse(True)
    except Exception:
      pass

    mol = Chem.MolFromSmiles("C")
    smarts = Chem.MolFromSmarts("C")
    try:
      Chem.ReplaceCore(mol, smarts, (0, 0))
      self.asssertFalse(True)
    except Exception:
      pass

  def test47RWMols(self):
    """ test the RWMol class

    """
    mol = Chem.MolFromSmiles('C1CCC1')
    self.assertTrue(type(mol) == Chem.Mol)

    for rwmol in [Chem.EditableMol(mol), Chem.RWMol(mol)]:
      self.assertTrue(type(rwmol) in [Chem.EditableMol, Chem.RWMol])
      newAt = Chem.Atom(8)
      rwmol.ReplaceAtom(0, newAt)
      self.assertTrue(Chem.MolToSmiles(rwmol.GetMol()) == 'C1COC1')

      rwmol.RemoveBond(0, 1)
      self.assertTrue(Chem.MolToSmiles(rwmol.GetMol()) == 'CCCO')
      a = Chem.Atom(7)
      idx = rwmol.AddAtom(a)
      self.assertEqual(rwmol.GetMol().GetNumAtoms(), 5)
      self.assertEqual(idx, 4)

      idx = rwmol.AddBond(0, 4, order=Chem.BondType.SINGLE)
      self.assertEqual(idx, 4)

      self.assertTrue(Chem.MolToSmiles(rwmol.GetMol()) == 'CCCON')
      rwmol.AddBond(4, 1, order=Chem.BondType.SINGLE)
      self.assertTrue(Chem.MolToSmiles(rwmol.GetMol()) == 'C1CNOC1')

      rwmol.RemoveAtom(3)
      self.assertTrue(Chem.MolToSmiles(rwmol.GetMol()) == 'CCNO')

      # practice shooting ourselves in the foot:
      m = Chem.MolFromSmiles('c1ccccc1')
      em = Chem.EditableMol(m)
      em.RemoveAtom(0)
      m2 = em.GetMol()
      self.assertRaises(ValueError, lambda: Chem.SanitizeMol(m2))
      m = Chem.MolFromSmiles('c1ccccc1')
      em = Chem.EditableMol(m)
      em.RemoveBond(0, 1)
      m2 = em.GetMol()
      self.assertRaises(ValueError, lambda: Chem.SanitizeMol(m2))

      # boundary cases:

      # removing non-existent bonds:
      m = Chem.MolFromSmiles('c1ccccc1')
      em = Chem.EditableMol(m)
      em.RemoveBond(0, 2)
      m2 = em.GetMol()
      Chem.SanitizeMol(m2)
      self.assertTrue(Chem.MolToSmiles(m2) == 'c1ccccc1')

      # removing non-existent atoms:
      m = Chem.MolFromSmiles('c1ccccc1')
      em = Chem.EditableMol(m)
      self.assertRaises(RuntimeError, lambda: em.RemoveAtom(12))

      # confirm that an RWMol can be constructed without arguments
      m = Chem.RWMol()

    # test replaceAtom/Bond preserving properties
    mol = Chem.MolFromSmiles('C1CCC1')
    mol2 = Chem.MolFromSmiles('C1CCC1')
    mol.GetAtomWithIdx(0).SetProp("foo", "bar")
    mol.GetBondWithIdx(0).SetProp("foo", "bar")
    newBond = mol2.GetBondWithIdx(0)
    self.assertTrue(type(mol) == Chem.Mol)

    for rwmol in [Chem.EditableMol(mol), Chem.RWMol(mol)]:
      newAt = Chem.Atom(8)
      rwmol.ReplaceAtom(0, newAt)
      self.assertTrue(Chem.MolToSmiles(rwmol.GetMol()) == 'C1COC1')
      self.assertFalse(rwmol.GetMol().GetAtomWithIdx(0).HasProp("foo"))

    for rwmol in [Chem.EditableMol(mol), Chem.RWMol(mol)]:
      newAt = Chem.Atom(8)
      rwmol.ReplaceAtom(0, newAt, preserveProps=True)
      self.assertTrue(Chem.MolToSmiles(rwmol.GetMol()) == 'C1COC1')
      self.assertTrue(rwmol.GetMol().GetAtomWithIdx(0).HasProp("foo"))
      self.assertEqual(rwmol.GetMol().GetAtomWithIdx(0).GetProp("foo"), "bar")

    for rwmol in [Chem.EditableMol(mol), Chem.RWMol(mol)]:
      rwmol.ReplaceBond(0, newBond)
      self.assertTrue(Chem.MolToSmiles(rwmol.GetMol()) == 'C1CCC1')
      self.assertFalse(rwmol.GetMol().GetBondWithIdx(0).HasProp("foo"))

    for rwmol in [Chem.EditableMol(mol), Chem.RWMol(mol)]:
      rwmol.ReplaceBond(0, newBond, preserveProps=True)
      self.assertTrue(Chem.MolToSmiles(rwmol.GetMol()) == 'C1CCC1')
      self.assertTrue(rwmol.GetMol().GetBondWithIdx(0).HasProp("foo"))
      self.assertEqual(rwmol.GetMol().GetBondWithIdx(0).GetProp("foo"), "bar")

  def test47SmartsPieces(self):
    """ test the GetAtomSmarts and GetBondSmarts functions

    """
    m = Chem.MolFromSmarts("[C,N]C")
    self.assertTrue(m.GetAtomWithIdx(0).GetSmarts() == '[C,N]')
    self.assertTrue(m.GetAtomWithIdx(1).GetSmarts() == 'C')
    self.assertEqual(m.GetBondBetweenAtoms(0, 1).GetSmarts(), '')

    m = Chem.MolFromSmarts("[$(C=O)]-O")
    self.assertTrue(m.GetAtomWithIdx(0).GetSmarts() == '[$(C=O)]')
    self.assertTrue(m.GetAtomWithIdx(1).GetSmarts() == 'O')
    self.assertTrue(m.GetBondBetweenAtoms(0, 1).GetSmarts() == '-')

    m = Chem.MolFromSmiles("CO")
    self.assertTrue(m.GetAtomWithIdx(0).GetSmarts() == 'C')
    self.assertTrue(m.GetAtomWithIdx(1).GetSmarts() == 'O')
    self.assertTrue(m.GetBondBetweenAtoms(0, 1).GetSmarts() == '')
    self.assertTrue(m.GetBondBetweenAtoms(0, 1).GetSmarts(allBondsExplicit=True) == '-')

    m = Chem.MolFromSmiles("C=O")
    self.assertTrue(m.GetAtomWithIdx(0).GetSmarts() == 'C')
    self.assertTrue(m.GetAtomWithIdx(1).GetSmarts() == 'O')
    self.assertTrue(m.GetBondBetweenAtoms(0, 1).GetSmarts() == '=')

    m = Chem.MolFromSmiles('C[C@H](F)[15NH3+]')
    self.assertEqual(m.GetAtomWithIdx(0).GetSmarts(), 'C')
    self.assertEqual(m.GetAtomWithIdx(0).GetSmarts(allHsExplicit=True), '[CH3]')
    self.assertEqual(m.GetAtomWithIdx(3).GetSmarts(), '[15NH3+]')
    self.assertEqual(m.GetAtomWithIdx(3).GetSmarts(allHsExplicit=True), '[15NH3+]')

  def test48Issue1928819(self):
    """ test a crash involving looping directly over mol suppliers
    """
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'NCI_aids_few.sdf')
    ms = [x for x in Chem.SDMolSupplier(fileN)]
    self.assertEqual(len(ms), 16)
    count = 0
    for m in Chem.SDMolSupplier(fileN):
      count += 1
    self.assertEqual(count, 16)

    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'fewSmi.csv')
    count = 0
    for m in Chem.SmilesMolSupplier(fileN, titleLine=False, smilesColumn=1, delimiter=','):
      count += 1
    self.assertEqual(count, 10)

    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'acd_few.tdt')
    count = 0
    for m in Chem.TDTMolSupplier(fileN):
      count += 1
    self.assertEqual(count, 10)

  def test49Issue1932365(self):
    """ test aromatic Se and Te from smiles/smarts
    """
    m = Chem.MolFromSmiles('c1ccc[se]1')
    self.assertTrue(m)
    self.assertTrue(m.GetAtomWithIdx(0).GetIsAromatic())
    self.assertTrue(m.GetAtomWithIdx(4).GetIsAromatic())
    m = Chem.MolFromSmiles('c1ccc[te]1')
    self.assertTrue(m)
    self.assertTrue(m.GetAtomWithIdx(0).GetIsAromatic())
    self.assertTrue(m.GetAtomWithIdx(4).GetIsAromatic())
    m = Chem.MolFromSmiles('C1=C[Se]C=C1')
    self.assertTrue(m)
    self.assertTrue(m.GetAtomWithIdx(0).GetIsAromatic())
    self.assertTrue(m.GetAtomWithIdx(2).GetIsAromatic())
    m = Chem.MolFromSmiles('C1=C[Te]C=C1')
    self.assertTrue(m)
    self.assertTrue(m.GetAtomWithIdx(0).GetIsAromatic())
    self.assertTrue(m.GetAtomWithIdx(2).GetIsAromatic())

    p = Chem.MolFromSmarts('[se]')
    self.assertTrue(Chem.MolFromSmiles('c1ccc[se]1').HasSubstructMatch(p))
    self.assertFalse(Chem.MolFromSmiles('C1=CCC[Se]1').HasSubstructMatch(p))

    p = Chem.MolFromSmarts('[te]')
    self.assertTrue(Chem.MolFromSmiles('c1ccc[te]1').HasSubstructMatch(p))
    self.assertFalse(Chem.MolFromSmiles('C1=CCC[Te]1').HasSubstructMatch(p))

  def test50Issue1968608(self):
    """ test sf.net issue 1968608
    """
    smarts = Chem.MolFromSmarts("[r5]")
    mol = Chem.MolFromSmiles("N12CCC36C1CC(C(C2)=CCOC4CC5=O)C4C3N5c7ccccc76")
    count = len(mol.GetSubstructMatches(smarts, uniquify=0))
    self.assertTrue(count == 9)

  def test51RadicalHandling(self):
    """ test handling of atoms with radicals
    """
    mol = Chem.MolFromSmiles("[C]C")
    self.assertTrue(mol)
    atom = mol.GetAtomWithIdx(0)
    self.assertTrue(atom.GetNumRadicalElectrons() == 3)
    self.assertTrue(atom.GetNoImplicit())
    atom.SetNoImplicit(False)
    atom.SetNumRadicalElectrons(1)
    mol.UpdatePropertyCache()
    self.assertTrue(atom.GetNumRadicalElectrons() == 1)
    self.assertTrue(atom.GetNumImplicitHs() == 2)

    mol = Chem.MolFromSmiles("[c]1ccccc1")
    self.assertTrue(mol)
    atom = mol.GetAtomWithIdx(0)
    self.assertTrue(atom.GetNumRadicalElectrons() == 1)
    self.assertTrue(atom.GetNoImplicit())

    mol = Chem.MolFromSmiles("[n]1ccccc1")
    self.assertTrue(mol)
    atom = mol.GetAtomWithIdx(0)
    self.assertTrue(atom.GetNumRadicalElectrons() == 0)
    self.assertTrue(atom.GetNoImplicit())

  def test52MolFrags(self):
    """ test GetMolFrags functionality
    """
    mol = Chem.MolFromSmiles("C.CC")
    self.assertTrue(mol)
    fs = Chem.GetMolFrags(mol)
    self.assertTrue(len(fs) == 2)
    self.assertTrue(len(fs[0]) == 1)
    self.assertTrue(tuple(fs[0]) == (0, ))
    self.assertTrue(len(fs[1]) == 2)
    self.assertTrue(tuple(fs[1]) == (1, 2))

    fs = Chem.GetMolFrags(mol, True)
    self.assertTrue(len(fs) == 2)
    self.assertTrue(fs[0].GetNumAtoms() == 1)
    self.assertTrue(fs[1].GetNumAtoms() == 2)

    mol = Chem.MolFromSmiles("CCC")
    self.assertTrue(mol)
    fs = Chem.GetMolFrags(mol)
    self.assertTrue(len(fs) == 1)
    self.assertTrue(len(fs[0]) == 3)
    self.assertTrue(tuple(fs[0]) == (0, 1, 2))
    fs = Chem.GetMolFrags(mol, True)
    self.assertTrue(len(fs) == 1)
    self.assertTrue(fs[0].GetNumAtoms() == 3)

    mol = Chem.MolFromSmiles("CO")
    em = Chem.EditableMol(mol)
    em.RemoveBond(0, 1)
    nm = em.GetMol()
    fs = Chem.GetMolFrags(nm, asMols=True)
    self.assertEqual([x.GetNumAtoms(onlyExplicit=False) for x in fs], [5, 3])
    fs = Chem.GetMolFrags(nm, asMols=True, sanitizeFrags=False)
    self.assertEqual([x.GetNumAtoms(onlyExplicit=False) for x in fs], [4, 2])

    mol = Chem.MolFromSmiles("CC.CCC")
    fs = Chem.GetMolFrags(mol, asMols=True)
    self.assertEqual([x.GetNumAtoms() for x in fs], [2, 3])
    frags = []
    fragsMolAtomMapping = []
    fs = Chem.GetMolFrags(mol, asMols=True, frags=frags, fragsMolAtomMapping=fragsMolAtomMapping)
    self.assertEqual(mol.GetNumAtoms(onlyExplicit=True), len(frags))
    fragsCheck = []
    for i, f in enumerate(fs):
      fragsCheck.extend([i] * f.GetNumAtoms(onlyExplicit=True))
    self.assertEqual(frags, fragsCheck)
    fragsMolAtomMappingCheck = []
    i = 0
    for f in fs:
      n = f.GetNumAtoms(onlyExplicit=True)
      fragsMolAtomMappingCheck.append(tuple(range(i, i + n)))
      i += n
    self.assertEqual(fragsMolAtomMapping, fragsMolAtomMappingCheck)

  def test53Matrices(self):
    """ test adjacency and distance matrices

    """
    m = Chem.MolFromSmiles('CC=C')
    d = Chem.GetDistanceMatrix(m, 0)
    self.assertTrue(feq(d[0, 1], 1.0))
    self.assertTrue(feq(d[0, 2], 2.0))
    self.assertTrue(feq(d[1, 0], 1.0))
    self.assertTrue(feq(d[2, 0], 2.0))
    a = Chem.GetAdjacencyMatrix(m, 0)
    self.assertTrue(a[0, 1] == 1)
    self.assertTrue(a[0, 2] == 0)
    self.assertTrue(a[1, 2] == 1)
    self.assertTrue(a[1, 0] == 1)
    self.assertTrue(a[2, 0] == 0)

    m = Chem.MolFromSmiles('C1CC1')
    d = Chem.GetDistanceMatrix(m, 0)
    self.assertTrue(feq(d[0, 1], 1.0))
    self.assertTrue(feq(d[0, 2], 1.0))
    a = Chem.GetAdjacencyMatrix(m, 0)
    self.assertTrue(a[0, 1] == 1)
    self.assertTrue(a[0, 2] == 1)
    self.assertTrue(a[1, 2] == 1)

    m = Chem.MolFromSmiles('CC.C')
    d = Chem.GetDistanceMatrix(m, 0)
    self.assertTrue(feq(d[0, 1], 1.0))
    self.assertTrue(d[0, 2] > 1000)
    self.assertTrue(d[1, 2] > 1000)
    a = Chem.GetAdjacencyMatrix(m, 0)
    self.assertTrue(a[0, 1] == 1)
    self.assertTrue(a[0, 2] == 0)
    self.assertTrue(a[1, 2] == 0)

  def test54Mol2Parser(self):
    """ test the mol2 parser
    """
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'pyrazole_pyridine.mol2')
    m = Chem.MolFromMol2File(fileN)
    self.assertTrue(m.GetNumAtoms() == 5)
    self.assertTrue(Chem.MolToSmiles(m) == 'c1cn[nH]c1', Chem.MolToSmiles(m))

    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         '3505.mol2')
    m = Chem.MolFromMol2File(fileN)
    self.assertTrue(m.GetBondBetweenAtoms(3, 12) is not None)
    self.assertEqual(m.GetBondBetweenAtoms(3, 12).GetBondType(), Chem.BondType.SINGLE)
    self.assertEqual(m.GetAtomWithIdx(12).GetFormalCharge(), 0)

    m = Chem.MolFromMol2File(fileN, cleanupSubstructures=False)
    self.assertTrue(m.GetBondBetweenAtoms(3, 12) is not None)
    self.assertEqual(m.GetBondBetweenAtoms(3, 12).GetBondType(), Chem.BondType.DOUBLE)
    self.assertEqual(m.GetAtomWithIdx(12).GetFormalCharge(), 1)

  def test55LayeredFingerprint(self):
    m1 = Chem.MolFromSmiles('CC(C)C')
    fp1 = Chem.LayeredFingerprint(m1)
    self.assertEqual(len(fp1), 2048)
    atomCounts = [0] * m1.GetNumAtoms()
    fp2 = Chem.LayeredFingerprint(m1, atomCounts=atomCounts)
    self.assertEqual(fp1, fp2)
    self.assertEqual(atomCounts, [4, 7, 4, 4])

    fp2 = Chem.LayeredFingerprint(m1, atomCounts=atomCounts)
    self.assertEqual(fp1, fp2)
    self.assertEqual(atomCounts, [8, 14, 8, 8])

    pbv = DataStructs.ExplicitBitVect(2048)
    fp3 = Chem.LayeredFingerprint(m1, setOnlyBits=pbv)
    self.assertEqual(fp3.GetNumOnBits(), 0)

    fp3 = Chem.LayeredFingerprint(m1, setOnlyBits=fp2)
    self.assertEqual(fp3, fp2)

    m2 = Chem.MolFromSmiles('CC')
    fp4 = Chem.LayeredFingerprint(m2)
    atomCounts = [0] * m1.GetNumAtoms()
    fp3 = Chem.LayeredFingerprint(m1, setOnlyBits=fp4, atomCounts=atomCounts)
    self.assertEqual(atomCounts, [1, 3, 1, 1])

    m2 = Chem.MolFromSmiles('CCC')
    fp4 = Chem.LayeredFingerprint(m2)
    atomCounts = [0] * m1.GetNumAtoms()
    fp3 = Chem.LayeredFingerprint(m1, setOnlyBits=fp4, atomCounts=atomCounts)
    self.assertEqual(atomCounts, [3, 6, 3, 3])

  def test56LazySDMolSupplier(self):
    if not hasattr(Chem, 'CompressedSDMolSupplier'):
      return

    self.assertRaises(ValueError, lambda: Chem.CompressedSDMolSupplier('nosuchfile.sdf.gz'))

    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'NCI_aids_few.sdf.gz')
    sdSup = Chem.CompressedSDMolSupplier(fileN)
    molNames = [
      "48", "78", "128", "163", "164", "170", "180", "186", "192", "203", "210", "211", "213",
      "220", "229", "256"
    ]

    chgs192 = {8: 1, 11: 1, 15: -1, 18: -1, 20: 1, 21: 1, 23: -1, 25: -1}
    i = 0
    for mol in sdSup:
      self.assertTrue(mol)
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      i += 1
      if (mol.GetProp("_Name") == "192"):
        # test parsed charges on one of the molecules
        for id in chgs192.keys():
          self.assertTrue(mol.GetAtomWithIdx(id).GetFormalCharge() == chgs192[id])
    self.assertEqual(i, 16)

    sdSup = Chem.CompressedSDMolSupplier(fileN)
    ns = [mol.GetProp("_Name") for mol in sdSup]
    self.assertTrue(ns == molNames)

    sdSup = Chem.CompressedSDMolSupplier(fileN, 0)
    for mol in sdSup:
      self.assertTrue(not mol.HasProp("numArom"))

  def test57AddRecursiveQuery(self):
    q1 = Chem.MolFromSmiles('CC')
    q2 = Chem.MolFromSmiles('CO')
    Chem.AddRecursiveQuery(q1, q2, 1)

    m1 = Chem.MolFromSmiles('OCC')
    self.assertTrue(m1.HasSubstructMatch(q2))
    self.assertTrue(m1.HasSubstructMatch(q1))
    self.assertTrue(m1.HasSubstructMatch(q1))
    self.assertTrue(m1.GetSubstructMatch(q1) == (2, 1))

    q3 = Chem.MolFromSmiles('CS')
    Chem.AddRecursiveQuery(q1, q3, 1)

    self.assertFalse(m1.HasSubstructMatch(q3))
    self.assertFalse(m1.HasSubstructMatch(q1))

    m2 = Chem.MolFromSmiles('OC(S)C')
    self.assertTrue(m2.HasSubstructMatch(q1))
    self.assertTrue(m2.GetSubstructMatch(q1) == (3, 1))

    m3 = Chem.MolFromSmiles('SCC')
    self.assertTrue(m3.HasSubstructMatch(q3))
    self.assertFalse(m3.HasSubstructMatch(q1))

    q1 = Chem.MolFromSmiles('CC')
    Chem.AddRecursiveQuery(q1, q2, 1)
    Chem.AddRecursiveQuery(q1, q3, 1, False)
    self.assertTrue(m3.HasSubstructMatch(q1))
    self.assertTrue(m3.GetSubstructMatch(q1) == (2, 1))

  def test58Issue2983794(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Wrap', 'test_data',
                         'issue2983794.sdf')
    m1 = Chem.MolFromMolFile(fileN)
    self.assertTrue(m1)
    em = Chem.EditableMol(m1)
    em.RemoveAtom(0)
    m2 = em.GetMol()
    Chem.Kekulize(m2)

  def test59Issue3007178(self):
    m = Chem.MolFromSmiles('CCC')
    a = m.GetAtomWithIdx(0)
    m = None
    self.assertEqual(Chem.MolToSmiles(a.GetOwningMol()), 'CCC')
    a = None
    m = Chem.MolFromSmiles('CCC')
    b = m.GetBondWithIdx(0)
    m = None
    self.assertEqual(Chem.MolToSmiles(b.GetOwningMol()), 'CCC')

  def test60SmilesWriterClose(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'fewSmi.csv')
    smiSup = Chem.SmilesMolSupplier(fileN, delimiter=",", smilesColumn=1, nameColumn=0, titleLine=0)
    ms = [x for x in smiSup]

    ofile = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Wrap', 'test_data',
                         'outSmiles.txt')
    writer = Chem.SmilesWriter(ofile)
    for mol in ms:
      writer.write(mol)
    writer.close()

    newsup = Chem.SmilesMolSupplier(ofile)
    newms = [x for x in newsup]
    self.assertEqual(len(ms), len(newms))

  def test61PathToSubmol(self):
    m = Chem.MolFromSmiles('CCCCCC1C(O)CC(O)N1C=CCO')
    env = Chem.FindAtomEnvironmentOfRadiusN(m, 2, 11)
    self.assertEqual(len(env), 8)
    amap = {}
    submol = Chem.PathToSubmol(m, env, atomMap=amap)
    self.assertEqual(submol.GetNumAtoms(), len(amap.keys()))
    self.assertEqual(submol.GetNumAtoms(), 9)
    smi = Chem.MolToSmiles(submol, rootedAtAtom=amap[11])
    self.assertEqual(smi[0], 'N')
    refsmi = Chem.MolToSmiles(Chem.MolFromSmiles('N(C=C)(C(C)C)C(O)C'))
    csmi = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
    self.assertEqual(refsmi, csmi)

  def test62SmilesAndSmartsReplacements(self):
    mol = Chem.MolFromSmiles('C{branch}C', replacements={'{branch}': 'C1(CC1)'})
    self.assertEqual(mol.GetNumAtoms(), 5)
    mol = Chem.MolFromSmarts('C{branch}C', replacements={'{branch}': 'C1(CC1)'})
    self.assertEqual(mol.GetNumAtoms(), 5)
    mol = Chem.MolFromSmiles('C{branch}C{acid}', replacements={
      '{branch}': 'C1(CC1)',
      '{acid}': "C(=O)O"
    })
    self.assertEqual(mol.GetNumAtoms(), 8)

  def test63Issue3313539(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'rgroups1.mol')
    m = Chem.MolFromMolFile(fileN)
    self.assertTrue(m is not None)
    at = m.GetAtomWithIdx(3)
    self.assertTrue(at is not None)
    self.assertTrue(at.HasProp('_MolFileRLabel'))
    p = at.GetProp('_MolFileRLabel')
    self.assertEqual(p, '2')
    self.assertEqual(Chem.GetAtomRLabel(at), 2)

    at = m.GetAtomWithIdx(4)
    self.assertTrue(at is not None)
    self.assertTrue(at.HasProp('_MolFileRLabel'))
    p = at.GetProp('_MolFileRLabel')
    self.assertEqual(p, '1')
    self.assertEqual(Chem.GetAtomRLabel(at), 1)

  def test64MoleculeCleanup(self):
    m = Chem.MolFromSmiles('CN(=O)=O', False)
    self.assertTrue(m)
    self.assertTrue(
      m.GetAtomWithIdx(1).GetFormalCharge() == 0 and m.GetAtomWithIdx(2).GetFormalCharge() == 0
      and m.GetAtomWithIdx(3).GetFormalCharge() == 0)
    self.assertTrue(
      m.GetBondBetweenAtoms(1, 3).GetBondType() == Chem.BondType.DOUBLE
      and m.GetBondBetweenAtoms(1, 2).GetBondType() == Chem.BondType.DOUBLE)
    Chem.Cleanup(m)
    m.UpdatePropertyCache()
    self.assertTrue(
      m.GetAtomWithIdx(1).GetFormalCharge() == 1 and
      (m.GetAtomWithIdx(2).GetFormalCharge() == -1 or m.GetAtomWithIdx(3).GetFormalCharge() == -1))
    self.assertTrue(
      m.GetBondBetweenAtoms(1, 3).GetBondType() == Chem.BondType.SINGLE
      or m.GetBondBetweenAtoms(1, 2).GetBondType() == Chem.BondType.SINGLE)

  def test65StreamSupplier(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'NCI_aids_few.sdf.gz')
    molNames = [
      "48", "78", "128", "163", "164", "170", "180", "186", "192", "203", "210", "211", "213",
      "220", "229", "256"
    ]
    inf = gzip.open(fileN)
    if 0:
      sb = Chem.streambuf(inf)
      suppl = Chem.ForwardSDMolSupplier(sb)
    else:
      suppl = Chem.ForwardSDMolSupplier(inf)

    i = 0
    while not suppl.atEnd():
      mol = next(suppl)
      self.assertTrue(mol)
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      i += 1
    self.assertEqual(i, 16)

    # make sure we have object ownership preserved
    inf = gzip.open(fileN)
    suppl = Chem.ForwardSDMolSupplier(inf)
    inf = None
    i = 0
    while not suppl.atEnd():
      mol = next(suppl)
      self.assertTrue(mol)
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      i += 1
    self.assertEqual(i, 16)

  @unittest.skipIf(not hasattr(Chem, 'MaeMolSupplier'), "not built with MAEParser support")
  def testMaeStreamSupplier(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'NCI_aids_few.maegz')
    molNames = [
      "48", "78", "128", "163", "164", "170", "180", "186", "192", "203", "210", "211", "213",
      "220", "229", "256"
    ]
    inf = gzip.open(fileN)
    suppl = Chem.MaeMolSupplier(inf)

    i = 0
    while not suppl.atEnd():
      mol = next(suppl)
      self.assertTrue(mol)
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      i += 1
    self.assertEqual(i, 16)

    # make sure we have object ownership preserved
    inf = gzip.open(fileN)
    suppl = Chem.MaeMolSupplier(inf)
    inf = None
    i = 0
    while not suppl.atEnd():
      mol = next(suppl)
      self.assertTrue(mol)
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      i += 1
    self.assertEqual(i, 16)

  @unittest.skipIf(not hasattr(Chem, 'MaeMolSupplier'), "not built with MAEParser support")
  def testMaeFileSupplier(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'NCI_aids_few.mae')
    molNames = [
      "48", "78", "128", "163", "164", "170", "180", "186", "192", "203", "210", "211", "213",
      "220", "229", "256"
    ]
    suppl = Chem.MaeMolSupplier(fileN)

    i = 0
    while not suppl.atEnd():
      mol = next(suppl)
      self.assertTrue(mol)
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      i += 1
    self.assertEqual(i, 16)

  @unittest.skipIf(not hasattr(Chem, 'MaeMolSupplier'), "not built with MAEParser support")
  def testMaeFileSupplierException(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'bad_ppty.mae')
    err_msg_substr = "Bad format for property"

    ok = False
    suppl = Chem.MaeMolSupplier(fileN)
    for i in range(5):
      try:
        mol = next(suppl)
      except RuntimeError as e:
        self.assertEqual(i, 1)
        self.assertTrue(err_msg_substr in str(e))
        ok = True
        break
      else:
        self.assertTrue(mol)
        self.assertTrue(mol.HasProp("_Name"))
        self.assertTrue(mol.GetNumAtoms() == 1)

    self.assertFalse(suppl.atEnd())
    self.assertTrue(ok)

  @unittest.skipIf(not hasattr(Chem, 'MaeMolSupplier'), "not built with MAEParser support")
  def testMaeFileSupplierSetData(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'NCI_aids_few.mae')
    molNames = [
      "48", "78", "128", "163", "164", "170", "180", "186", "192", "203", "210", "211", "213",
      "220", "229", "256"
    ]
    suppl = Chem.MaeMolSupplier()

    with open(fileN) as f:
      data = f.read()

    suppl.SetData(data)

    self.assertEqual(len(suppl), 16)

    for i, mol in enumerate(suppl):
      self.assertTrue(mol)
      self.assertEqual(mol.GetProp("_Name"), molNames[i])

    self.assertEqual(i, 15)

    # Do it again, to check the reset() method
    suppl.reset()

    self.assertEqual(len(suppl), 16)

    for i, mol in enumerate(suppl):
      self.assertTrue(mol)
      self.assertEqual(mol.GetProp("_Name"), molNames[i])

    self.assertEqual(i, 15)

  @unittest.skipIf(not hasattr(Chem, 'MaeMolSupplier'), "not built with MAEParser support")
  def testMaeFileSupplierGetItem(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'NCI_aids_few.mae')
    molNames = [
      "48", "78", "128", "163", "164", "170", "180", "186", "192", "203", "210", "211", "213",
      "220", "229", "256"
    ]
    with Chem.MaeMolSupplier(fileN) as suppl:
      num_mols = len(suppl)
      self.assertEqual(num_mols, 16)

      for i in range(num_mols):
        self.assertEqual(suppl[i].GetProp("_Name"), molNames[i])

        j = -(i + 1)
        self.assertEqual(suppl[j].GetProp("_Name"), molNames[j])

      for i in (num_mols, 2 * num_mols, -(num_mols + 1), -2 * (num_mols + 1)):
        self.assertRaises(IndexError, lambda: suppl[i])

  @unittest.skipIf(not hasattr(Chem, 'MaeMolSupplier'), "not built with MAEParser support")
  def testMaeFileSupplierExceptionMsgs(self):
    maeBlock = "f_m_ct {\n  s_m_title\n  :::\n  " "\n  }\n}"

    with Chem.MaeMolSupplier() as suppl:
      suppl.SetData(maeBlock)
      self.assertRaisesRegex(RuntimeError, r'File parsing error: Indexed block not found: m_atom',
                             lambda: next(suppl))

  def test66StreamSupplierIter(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'NCI_aids_few.sdf.gz')
    inf = gzip.open(fileN)
    if 0:
      sb = Chem.streambuf(inf)
      suppl = Chem.ForwardSDMolSupplier(sb)
    else:
      suppl = Chem.ForwardSDMolSupplier(inf)

    molNames = [
      "48", "78", "128", "163", "164", "170", "180", "186", "192", "203", "210", "211", "213",
      "220", "229", "256"
    ]
    i = 0
    for mol in suppl:
      self.assertTrue(mol)
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      i += 1
    self.assertEqual(i, 16)

  def test67StreamSupplierStringIO(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'NCI_aids_few.sdf.gz')
    sio = BytesIO(gzip.open(fileN).read())
    suppl = Chem.ForwardSDMolSupplier(sio)
    molNames = [
      "48", "78", "128", "163", "164", "170", "180", "186", "192", "203", "210", "211", "213",
      "220", "229", "256"
    ]
    i = 0
    for mol in suppl:
      self.assertTrue(mol)
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      i += 1
    self.assertEqual(i, 16)

  def test68ForwardSupplierUsingFilename(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'NCI_aids_few.sdf')
    suppl = Chem.ForwardSDMolSupplier(fileN)
    molNames = [
      "48", "78", "128", "163", "164", "170", "180", "186", "192", "203", "210", "211", "213",
      "220", "229", "256"
    ]
    i = 0
    for mol in suppl:
      self.assertTrue(mol)
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      i += 1
    self.assertEqual(i, 16)

    self.assertRaises(IOError, lambda: Chem.ForwardSDMolSupplier('nosuchfile.sdf'))

  def test69StreamSupplierStreambuf(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'NCI_aids_few.sdf.gz')
    sb = rdBase.streambuf(gzip.open(fileN))
    suppl = Chem.ForwardSDMolSupplier(sb)

    molNames = [
      "48", "78", "128", "163", "164", "170", "180", "186", "192", "203", "210", "211", "213",
      "220", "229", "256"
    ]
    i = 0
    for mol in suppl:
      self.assertTrue(mol)
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      i += 1
    self.assertEqual(i, 16)

  def test70StreamSDWriter(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'NCI_aids_few.sdf.gz')
    inf = gzip.open(fileN)
    suppl = Chem.ForwardSDMolSupplier(inf)
    osio = StringIO()
    w = Chem.SDWriter(osio)
    molNames = [
      "48", "78", "128", "163", "164", "170", "180", "186", "192", "203", "210", "211", "213",
      "220", "229", "256"
    ]
    i = 0
    for mol in suppl:
      self.assertTrue(mol)
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      w.write(mol)
      i += 1
    self.assertEqual(i, 16)
    w.flush()
    w = None
    txt = osio.getvalue().encode()
    isio = BytesIO(txt)
    suppl = Chem.ForwardSDMolSupplier(isio)
    i = 0
    for mol in suppl:
      self.assertTrue(mol)
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      i += 1
    self.assertEqual(i, 16)

  def test71StreamSmilesWriter(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'esters.sdf')
    suppl = Chem.ForwardSDMolSupplier(fileN)
    osio = StringIO()
    w = Chem.SmilesWriter(osio)
    ms = [x for x in suppl]
    w.SetProps(ms[0].GetPropNames())
    i = 0
    for mol in ms:
      self.assertTrue(mol)
      w.write(mol)
      i += 1
    self.assertEqual(i, 6)
    w.flush()
    w = None
    txt = osio.getvalue()
    self.assertEqual(txt.count('ID'), 1)
    self.assertEqual(txt.count('\n'), 7)

  def test72StreamTDTWriter(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'esters.sdf')
    suppl = Chem.ForwardSDMolSupplier(fileN)
    osio = StringIO()
    w = Chem.TDTWriter(osio)
    ms = [x for x in suppl]
    w.SetProps(ms[0].GetPropNames())
    i = 0
    for mol in ms:
      self.assertTrue(mol)
      w.write(mol)
      i += 1
    self.assertEqual(i, 6)
    w.flush()
    w = None
    txt = osio.getvalue()
    self.assertEqual(txt.count('ID'), 6)
    self.assertEqual(txt.count('NAME'), 6)

  def test73SanitizationOptions(self):
    m = Chem.MolFromSmiles('c1ccccc1', sanitize=False)
    res = Chem.SanitizeMol(m, catchErrors=True)
    self.assertEqual(res, 0)

    m = Chem.MolFromSmiles('c1cccc1', sanitize=False)
    res = Chem.SanitizeMol(m, catchErrors=True)
    self.assertEqual(res, Chem.SanitizeFlags.SANITIZE_KEKULIZE)

    m = Chem.MolFromSmiles('CC(C)(C)(C)C', sanitize=False)
    res = Chem.SanitizeMol(m, catchErrors=True)
    self.assertEqual(res, Chem.SanitizeFlags.SANITIZE_PROPERTIES)

    m = Chem.MolFromSmiles('c1cccc1', sanitize=False)
    res = Chem.SanitizeMol(
      m, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
      catchErrors=True)
    self.assertEqual(res, Chem.SanitizeFlags.SANITIZE_NONE)

  def test74Issue3510149(self):
    mol = Chem.MolFromSmiles("CCC1CNCC1CC")
    atoms = mol.GetAtoms()
    mol = None
    for atom in atoms:
      idx = atom.GetIdx()
      p = atom.GetOwningMol().GetNumAtoms()

    mol = Chem.MolFromSmiles("CCC1CNCC1CC")
    bonds = mol.GetBonds()
    mol = None
    for bond in bonds:
      idx = bond.GetIdx()
      p = atom.GetOwningMol().GetNumAtoms()

    mol = Chem.MolFromSmiles("CCC1CNCC1CC")
    bond = mol.GetBondBetweenAtoms(0, 1)
    mol = None
    idx = bond.GetBeginAtomIdx()
    p = bond.GetOwningMol().GetNumAtoms()

    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'NCI_aids_few.sdf')
    sdSup = Chem.SDMolSupplier(fileN)
    mol = next(sdSup)
    nats = mol.GetNumAtoms()
    conf = mol.GetConformer()
    mol = None
    self.assertEqual(nats, conf.GetNumAtoms())
    conf.GetOwningMol().GetProp("_Name")

  def test75AllBondsExplicit(self):
    m = Chem.MolFromSmiles("CCC")
    smi = Chem.MolToSmiles(m)
    self.assertEqual(smi, "CCC")
    smi = Chem.MolToSmiles(m, allBondsExplicit=True)
    self.assertEqual(smi, "C-C-C")

    m = Chem.MolFromSmiles("c1ccccc1")
    smi = Chem.MolToSmiles(m)
    self.assertEqual(smi, "c1ccccc1")
    smi = Chem.MolToSmiles(m, allBondsExplicit=True)
    self.assertEqual(smi, "c1:c:c:c:c:c:1")

  def test76VeryLargeMolecule(self):
    # this is sf.net issue 3524984
    smi = '[C@H](F)(Cl)' + 'c1cc[nH]c1' * 500 + '[C@H](F)(Cl)'
    m = Chem.MolFromSmiles(smi)
    self.assertTrue(m)
    self.assertEqual(m.GetNumAtoms(), 2506)
    scs = Chem.FindMolChiralCenters(m)
    self.assertEqual(len(scs), 2)

  def test77MolFragmentToSmiles(self):
    smi = "OC1CC1CC"
    m = Chem.MolFromSmiles(smi)
    fsmi = Chem.MolFragmentToSmiles(m, [1, 2, 3])
    self.assertEqual(fsmi, "C1CC1")
    fsmi = Chem.MolFragmentToSmiles(m, [1, 2, 3], bondsToUse=[1, 2, 5])
    self.assertEqual(fsmi, "C1CC1")
    fsmi = Chem.MolFragmentToSmiles(m, [1, 2, 3], bondsToUse=[1, 2])
    self.assertEqual(fsmi, "CCC")
    fsmi = Chem.MolFragmentToSmiles(m, [1, 2, 3], atomSymbols=["", "[A]", "[C]", "[B]", "", ""])
    self.assertEqual(fsmi, "[A]1[B][C]1")
    fsmi = Chem.MolFragmentToSmiles(m, [1, 2, 3], bondSymbols=["", "%", "%", "", "", "%"])
    self.assertEqual(fsmi, "C1%C%C%1")

    smi = "c1ccccc1C"
    m = Chem.MolFromSmiles(smi)
    fsmi = Chem.MolFragmentToSmiles(m, range(6))
    self.assertEqual(fsmi, "c1ccccc1")
    Chem.Kekulize(m)
    fsmi = Chem.MolFragmentToSmiles(m, range(6), kekuleSmiles=True)
    self.assertEqual(fsmi, "C1=CC=CC=C1")
    fsmi = Chem.MolFragmentToSmiles(m, range(6), atomSymbols=["[C]"] * 7, kekuleSmiles=True)
    self.assertEqual(fsmi, "[C]1=[C][C]=[C][C]=[C]1")

    self.assertRaises(ValueError, lambda: Chem.MolFragmentToSmiles(m, []))

  def test78AtomAndBondProps(self):
    m = Chem.MolFromSmiles('c1ccccc1')
    at = m.GetAtomWithIdx(0)
    self.assertFalse(at.HasProp('foo'))
    at.SetProp('foo', 'bar')
    self.assertTrue(at.HasProp('foo'))
    self.assertEqual(at.GetProp('foo'), 'bar')
    bond = m.GetBondWithIdx(0)
    self.assertFalse(bond.HasProp('foo'))
    bond.SetProp('foo', 'bar')
    self.assertTrue(bond.HasProp('foo'))
    self.assertEqual(bond.GetProp('foo'), 'bar')

  def test79AddRecursiveStructureQueries(self):
    qs = {'carbonyl': Chem.MolFromSmiles('CO'), 'amine': Chem.MolFromSmiles('CN')}
    q = Chem.MolFromSmiles('CCC')
    q.GetAtomWithIdx(0).SetProp('query', 'carbonyl,amine')
    Chem.MolAddRecursiveQueries(q, qs, 'query')
    m = Chem.MolFromSmiles('CCCO')
    self.assertTrue(m.HasSubstructMatch(q))
    m = Chem.MolFromSmiles('CCCN')
    self.assertTrue(m.HasSubstructMatch(q))
    m = Chem.MolFromSmiles('CCCC')
    self.assertFalse(m.HasSubstructMatch(q))

  def test80ParseMolQueryDefFile(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'ChemTransforms', 'testData',
                         'query_file1.txt')
    d = Chem.ParseMolQueryDefFile(fileN, standardize=False)
    self.assertTrue('CarboxylicAcid' in d)
    m = Chem.MolFromSmiles('CC(=O)O')
    self.assertTrue(m.HasSubstructMatch(d['CarboxylicAcid']))
    self.assertFalse(m.HasSubstructMatch(d['CarboxylicAcid.Aromatic']))

    d = Chem.ParseMolQueryDefFile(fileN)
    self.assertTrue('carboxylicacid' in d)
    self.assertFalse('CarboxylicAcid' in d)

  def test81Issue275(self):
    smi = Chem.MolToSmiles(Chem.MurckoDecompose(
      Chem.MolFromSmiles('CCCCC[C@H]1CC[C@H](C(=O)O)CC1')))
    self.assertEqual(smi, 'C1CCCCC1')

  def test82Issue288(self):
    m = Chem.MolFromSmiles('CC*')
    m.GetAtomWithIdx(2).SetProp('molAtomMapNumber', '30')

    smi = Chem.MolToSmiles(m)
    self.assertEqual(smi, 'CC[*:30]')
    # try newer api
    m = Chem.MolFromSmiles('CC*')
    m.GetAtomWithIdx(2).SetAtomMapNum(30)
    smi = Chem.MolToSmiles(m)
    self.assertEqual(smi, 'CC[*:30]')

  def test83GitHubIssue19(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'empty2.sdf')
    with self.assertRaises(OSError):
      sdSup = Chem.SDMolSupplier(fileN)

    sdSup = Chem.SDMolSupplier()
    sdSup.SetData('')
    self.assertTrue(sdSup.atEnd())
    self.assertRaises(IndexError, lambda: sdSup[0])

    sdSup.SetData('')
    self.assertRaises(IndexError, lambda: sdSup[0])
    self.assertEqual(len(sdSup), 0)

  def test84PDBBasics(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         '1CRN.pdb')
    m = Chem.MolFromPDBFile(fileN, proximityBonding=False)
    self.assertEqual(m.GetNumAtoms(), 327)
    self.assertEqual(m.GetNumBonds(), 3)
    m = Chem.MolFromPDBFile(fileN)
    self.assertTrue(m is not None)
    self.assertEqual(m.GetNumAtoms(), 327)
    self.assertEqual(m.GetNumBonds(), 337)
    self.assertTrue(m.GetAtomWithIdx(0).GetPDBResidueInfo())
    self.assertEqual(m.GetAtomWithIdx(0).GetPDBResidueInfo().GetName(), " N  ")
    self.assertEqual(m.GetAtomWithIdx(0).GetPDBResidueInfo().GetResidueName(), "THR")
    self.assertAlmostEqual(m.GetAtomWithIdx(0).GetPDBResidueInfo().GetTempFactor(), 13.79, 2)
    m = Chem.MolFromPDBBlock(Chem.MolToPDBBlock(m))
    self.assertEqual(m.GetNumAtoms(), 327)
    self.assertEqual(m.GetNumBonds(), 337)
    self.assertTrue(m.GetAtomWithIdx(0).GetPDBResidueInfo())
    self.assertEqual(m.GetAtomWithIdx(0).GetPDBResidueInfo().GetName(), " N  ")
    self.assertEqual(m.GetAtomWithIdx(0).GetPDBResidueInfo().GetResidueName(), "THR")
    self.assertAlmostEqual(m.GetAtomWithIdx(0).GetPDBResidueInfo().GetTempFactor(), 13.79, 2)
    # test multivalent Hs
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         '2c92_hypervalentH.pdb')
    mol = Chem.MolFromPDBFile(fileN, sanitize=False, removeHs=False)
    atom = mol.GetAtomWithIdx(84)
    self.assertEqual(atom.GetAtomicNum(), 1)  # is it H
    self.assertEqual(atom.GetDegree(), 1)  # H should have 1 bond
    for n in atom.GetNeighbors():  # Check if neighbor is from the same residue
      self.assertEqual(atom.GetPDBResidueInfo().GetResidueName(),
                       n.GetPDBResidueInfo().GetResidueName())
    # test unbinding metals (ZN)
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         '1ps3_zn.pdb')
    mol = Chem.MolFromPDBFile(fileN, sanitize=False, removeHs=False)
    atom = mol.GetAtomWithIdx(40)
    self.assertEqual(atom.GetAtomicNum(), 30)  # is it Zn
    self.assertEqual(atom.GetDegree(), 4)  # Zn should have 4 zero-order bonds
    self.assertEqual(atom.GetExplicitValence(), 0)
    bonds_order = [bond.GetBondType() for bond in atom.GetBonds()]
    self.assertEqual(bonds_order, [Chem.BondType.ZERO] * atom.GetDegree())

    # test metal bonds without proximity bonding
    mol = Chem.MolFromPDBFile(fileN, sanitize=False, removeHs=False, proximityBonding=False)
    atom = mol.GetAtomWithIdx(40)
    self.assertEqual(atom.GetAtomicNum(), 30)  # is it Zn
    self.assertEqual(atom.GetDegree(), 4)  # Zn should have 4 zero-order bonds
    self.assertEqual(atom.GetExplicitValence(), 0)
    bonds_order = [bond.GetBondType() for bond in atom.GetBonds()]
    self.assertEqual(bonds_order, [Chem.BondType.ZERO] * atom.GetDegree())
    # test unbinding HOHs
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         '2vnf_bindedHOH.pdb')
    mol = Chem.MolFromPDBFile(fileN, sanitize=False, removeHs=False)
    atom = mol.GetAtomWithIdx(10)
    self.assertEqual(atom.GetPDBResidueInfo().GetResidueName(), 'HOH')
    self.assertEqual(atom.GetDegree(), 0)  # HOH should have no bonds
    # test metal bonding in ligand
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         '2dej_APW.pdb')
    mol = Chem.MolFromPDBFile(fileN, sanitize=False, removeHs=False)
    atom = mol.GetAtomWithIdx(6)
    self.assertEqual(atom.GetAtomicNum(), 12)
    self.assertEqual(atom.GetDegree(), 2)
    atom = mol.GetAtomWithIdx(35)
    self.assertEqual(atom.GetPDBResidueInfo().GetResidueName(), 'HOH')
    self.assertEqual(atom.GetDegree(), 0)

  def test85AtomCopying(self):
    """Can a copied atom be added to a molecule?"""
    import copy
    m = Chem.MolFromSmiles('C1CC1')
    a = m.GetAtomWithIdx(0)
    a_copy1 = copy.copy(a)
    a_copy2 = Chem.Atom(a)
    m = None
    a = None

    def assert_is_valid_atom(a):
      new_m = Chem.RWMol()
      new_m.AddAtom(a)
      # This will not match if the owning mol is unset for a_copy,
      # or if there has been a clean up.
      self.assertEqual(new_m.GetAtomWithIdx(0).GetIdx(), 0)

    assert_is_valid_atom(a_copy1)
    assert_is_valid_atom(a_copy2)

  def test85MolCopying(self):
    m = Chem.MolFromSmiles('C1CC1[C@H](F)Cl')
    m.SetProp('foo', 'bar')
    m2 = Chem.Mol(m)
    self.assertEqual(Chem.MolToSmiles(m, True), Chem.MolToSmiles(m2, True))
    self.assertTrue(m2.HasProp('foo'))
    self.assertEqual(m2.GetProp('foo'), 'bar')
    ri = m2.GetRingInfo()
    self.assertTrue(ri)
    self.assertEqual(ri.NumRings(), 1)

  def test85MolCopying2(self):
    import copy
    m1 = Chem.MolFromSmiles('CC')
    m1.SetProp('Foo', 'bar')
    m1.foo = [1]
    m2 = copy.copy(m1)
    m3 = copy.copy(m2)
    m4 = copy.deepcopy(m1)
    m5 = copy.deepcopy(m2)
    m6 = copy.deepcopy(m4)

    self.assertEqual(m1.GetProp('Foo'), 'bar')
    self.assertEqual(m2.GetProp('Foo'), 'bar')
    self.assertEqual(m3.GetProp('Foo'), 'bar')
    self.assertEqual(m4.GetProp('Foo'), 'bar')
    self.assertEqual(m5.GetProp('Foo'), 'bar')
    self.assertEqual(m6.GetProp('Foo'), 'bar')

    m2.foo.append(4)
    self.assertEqual(m1.foo, [1, 4])
    self.assertEqual(m2.foo, [1, 4])
    self.assertEqual(m3.foo, [1, 4])
    self.assertEqual(m4.foo, [1])
    self.assertEqual(m5.foo, [1])
    self.assertEqual(m6.foo, [1])

    m7 = Chem.RWMol(m1)
    self.assertFalse(hasattr(m7, 'foo'))
    m7.foo = [1]
    m8 = copy.copy(m7)
    m9 = copy.deepcopy(m7)
    m8.foo.append(4)
    self.assertEqual(m7.GetProp('Foo'), 'bar')
    self.assertEqual(m8.GetProp('Foo'), 'bar')
    self.assertEqual(m9.GetProp('Foo'), 'bar')
    self.assertEqual(m8.foo, [1, 4])
    self.assertEqual(m9.foo, [1])

  def test86MolRenumbering(self):
    import random
    m = Chem.MolFromSmiles('C[C@H]1CC[C@H](C/C=C/[C@H](F)Cl)CC1')
    cSmi = Chem.MolToSmiles(m, True)
    for i in range(m.GetNumAtoms()):
      ans = list(range(m.GetNumAtoms()))
      random.shuffle(ans)
      m2 = Chem.RenumberAtoms(m, ans)
      nSmi = Chem.MolToSmiles(m2, True)
      self.assertEqual(cSmi, nSmi)

  def test87FragmentOnBonds(self):
    m = Chem.MolFromSmiles('CC1CC(O)C1CCC1CC1')
    bis = m.GetSubstructMatches(Chem.MolFromSmarts('[!R][R]'))
    bs = []
    labels = []
    for bi in bis:
      b = m.GetBondBetweenAtoms(bi[0], bi[1])
      if b.GetBeginAtomIdx() == bi[0]:
        labels.append((10, 1))
      else:
        labels.append((1, 10))
      bs.append(b.GetIdx())
    nm = Chem.FragmentOnBonds(m, bs)
    frags = Chem.GetMolFrags(nm)
    self.assertEqual(len(frags), 5)
    self.assertEqual(frags,
                     ((0, 12), (1, 2, 3, 5, 11, 14, 16), (4, 13), (6, 7, 15, 18), (8, 9, 10, 17)))
    smi = Chem.MolToSmiles(nm, True)
    self.assertEqual(smi, '*C1CC([4*])C1[6*].[1*]C.[3*]O.[5*]CC[8*].[7*]C1CC1')

    nm = Chem.FragmentOnBonds(m, bs, dummyLabels=labels)
    frags = Chem.GetMolFrags(nm)
    self.assertEqual(len(frags), 5)
    self.assertEqual(frags,
                     ((0, 12), (1, 2, 3, 5, 11, 14, 16), (4, 13), (6, 7, 15, 18), (8, 9, 10, 17)))
    smi = Chem.MolToSmiles(nm, True)
    self.assertEqual(smi, '[1*]C.[1*]CC[1*].[1*]O.[10*]C1CC([10*])C1[10*].[10*]C1CC1')

    m = Chem.MolFromSmiles('CCC(=O)CC(=O)C')
    bis = m.GetSubstructMatches(Chem.MolFromSmarts('C=O'))
    bs = []
    for bi in bis:
      b = m.GetBondBetweenAtoms(bi[0], bi[1])
      bs.append(b.GetIdx())
    bts = [Chem.BondType.DOUBLE] * len(bs)
    nm = Chem.FragmentOnBonds(m, bs, bondTypes=bts)
    frags = Chem.GetMolFrags(nm)
    self.assertEqual(len(frags), 3)
    smi = Chem.MolToSmiles(nm, True)
    self.assertEqual(smi, '[2*]=O.[3*]=C(CC)CC(=[6*])C.[5*]=O')

    # github issue 430:
    m = Chem.MolFromSmiles('OCCCCN')
    self.assertRaises(ValueError, lambda: Chem.FragmentOnBonds(m, ()))

  def test88QueryAtoms(self):
    from rdkit.Chem import rdqueries
    m = Chem.MolFromSmiles('c1nc(C)n(CC)c1')

    qa = rdqueries.ExplicitDegreeEqualsQueryAtom(3)
    l = tuple([x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)])
    self.assertEqual(l, (2, 4))

    qa.ExpandQuery(rdqueries.AtomNumEqualsQueryAtom(6, negate=True))
    l = tuple([x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)])
    self.assertEqual(l, (4, ))

    qa = rdqueries.ExplicitDegreeEqualsQueryAtom(3)
    qa.ExpandQuery(rdqueries.AtomNumEqualsQueryAtom(6, negate=True),
                   how=Chem.CompositeQueryType.COMPOSITE_OR)
    l = tuple([x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)])
    self.assertEqual(l, (1, 2, 4))

    qa = rdqueries.ExplicitDegreeEqualsQueryAtom(3)
    qa.ExpandQuery(rdqueries.AtomNumEqualsQueryAtom(6, negate=True),
                   how=Chem.CompositeQueryType.COMPOSITE_XOR)
    l = tuple([x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)])
    self.assertEqual(l, (1, 2))

    qa = rdqueries.ExplicitDegreeGreaterQueryAtom(2)
    l = tuple([x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)])
    self.assertEqual(l, (2, 4))

    qa = rdqueries.ExplicitDegreeLessQueryAtom(2)
    l = tuple([x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)])
    self.assertEqual(l, (3, 6))

    m = Chem.MolFromSmiles('N[CH][CH]')
    qa = rdqueries.NumRadicalElectronsGreaterQueryAtom(0)
    l = tuple([x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)])
    self.assertEqual(l, (1, 2))
    qa = rdqueries.NumRadicalElectronsGreaterQueryAtom(1)
    l = tuple([x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)])
    self.assertEqual(l, (2, ))

    m = Chem.MolFromSmiles('F[C@H](Cl)C')
    qa = rdqueries.HasChiralTagQueryAtom()
    l = tuple([x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)])
    self.assertEqual(l, (1, ))
    qa = rdqueries.MissingChiralTagQueryAtom()
    l = tuple([x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)])
    self.assertEqual(l, ())

    m = Chem.MolFromSmiles('F[CH](Cl)C')
    qa = rdqueries.HasChiralTagQueryAtom()
    l = tuple([x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)])
    self.assertEqual(l, ())
    qa = rdqueries.MissingChiralTagQueryAtom()
    l = tuple([x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)])
    self.assertEqual(l, (1, ))

    m = Chem.MolFromSmiles('CNCON')
    qa = rdqueries.NumHeteroatomNeighborsEqualsQueryAtom(2)
    l = tuple([x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)])
    self.assertEqual(l, (2, ))
    qa = rdqueries.NumHeteroatomNeighborsGreaterQueryAtom(0)
    l = tuple([x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)])
    self.assertEqual(l, (0, 2, 3, 4))

    m = Chem.MolFromSmiles('CC12CCN(CC1)C2')
    qa = rdqueries.IsBridgeheadQueryAtom()
    l = tuple([x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)])
    self.assertEqual(l, (1, 4))

    m = Chem.MolFromSmiles('OCCOC')
    qa = rdqueries.NonHydrogenDegreeEqualsQueryAtom(2)
    l = tuple([x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)])
    self.assertEqual(l, (1, 2, 3))

  def test89UnicodeInput(self):
    m = Chem.MolFromSmiles(u'c1ccccc1')
    self.assertTrue(m is not None)
    self.assertEqual(m.GetNumAtoms(), 6)
    m = Chem.MolFromSmarts(u'c1ccccc1')
    self.assertTrue(m is not None)
    self.assertEqual(m.GetNumAtoms(), 6)

  def test90FragmentOnSomeBonds(self):
    m = Chem.MolFromSmiles('OCCCCN')
    pieces = Chem.FragmentOnSomeBonds(m, (0, 2, 4), 2)
    self.assertEqual(len(pieces), 3)

    frags = Chem.GetMolFrags(pieces[0])
    self.assertEqual(len(frags), 3)
    self.assertEqual(len(frags[0]), 2)
    self.assertEqual(len(frags[1]), 4)
    self.assertEqual(len(frags[2]), 4)

    frags = Chem.GetMolFrags(pieces[1])
    self.assertEqual(len(frags), 3)
    self.assertEqual(len(frags[0]), 2)
    self.assertEqual(len(frags[1]), 6)
    self.assertEqual(len(frags[2]), 2)

    frags = Chem.GetMolFrags(pieces[2])
    self.assertEqual(len(frags), 3)
    self.assertEqual(len(frags[0]), 4)
    self.assertEqual(len(frags[1]), 4)
    self.assertEqual(len(frags[2]), 2)

    pieces, cpa = Chem.FragmentOnSomeBonds(m, (0, 2, 4), 2, returnCutsPerAtom=True)
    self.assertEqual(len(pieces), 3)
    self.assertEqual(len(cpa), 3)
    self.assertEqual(len(cpa[0]), m.GetNumAtoms())

    # github issue 430:
    m = Chem.MolFromSmiles('OCCCCN')
    self.assertRaises(ValueError, lambda: Chem.FragmentOnSomeBonds(m, ()))

    pieces = Chem.FragmentOnSomeBonds(m, (0, 2, 4), 0)
    self.assertEqual(len(pieces), 0)

  def test91RankAtoms(self):
    m = Chem.MolFromSmiles('ONCS.ONCS')
    ranks = Chem.CanonicalRankAtoms(m, breakTies=False)
    self.assertEqual(list(ranks[0:4]), list(ranks[4:]))

    m = Chem.MolFromSmiles("c1ccccc1")
    ranks = Chem.CanonicalRankAtoms(m, breakTies=False)
    for x in ranks:
      self.assertEqual(x, 0)

    m = Chem.MolFromSmiles("C1NCN1")
    ranks = Chem.CanonicalRankAtoms(m, breakTies=False)
    self.assertEqual(ranks[0], ranks[2])
    self.assertEqual(ranks[1], ranks[3])

  def test92RankAtomsInFragment(self):
    m = Chem.MolFromSmiles('ONCS.ONCS')
    ranks = Chem.CanonicalRankAtomsInFragment(m, [0, 1, 2, 3], [0, 1, 2])

    ranks2 = Chem.CanonicalRankAtomsInFragment(m, [4, 5, 6, 7], [3, 4, 5])
    self.assertEqual(list(ranks[0:4]), list(ranks2[4:]))
    self.assertEqual(list(ranks[4:]), [-1] * 4)
    self.assertEqual(list(ranks2[0:4]), [-1] * 4)

    # doc tests
    mol = Chem.MolFromSmiles('C1NCN1.C1NCN1')
    self.assertEqual(
      list(Chem.CanonicalRankAtomsInFragment(mol, atomsToUse=range(0, 4), breakTies=False)),
      [4, 6, 4, 6, -1, -1, -1, -1])
    self.assertNotEqual(
      list(Chem.CanonicalRankAtomsInFragment(mol, atomsToUse=range(0, 4), breakTies=True)),
      [4, 6, 4, 6, -1, -1, -1, -1])
    self.assertEqual(
      list(Chem.CanonicalRankAtomsInFragment(mol, atomsToUse=range(4, 8), breakTies=False)),
      [-1, -1, -1, -1, 4, 6, 4, 6])
    self.assertNotEqual(
      list(Chem.CanonicalRankAtomsInFragment(mol, atomsToUse=range(4, 8), breakTies=True)),
      [-1, -1, -1, -1, 4, 6, 4, 6])

  def test93RWMolsAsROMol(self):
    """ test the RWMol class as a proper ROMol

    """
    mol = Chem.MolFromSmiles('C1CCC1')
    self.assertTrue(type(mol) == Chem.Mol)
    rwmol = Chem.RWMol(mol)
    self.assertEqual(Chem.MolToSmiles(rwmol, True), Chem.MolToSmiles(rwmol.GetMol()))
    newAt = Chem.Atom(8)
    rwmol.ReplaceAtom(0, newAt)
    self.assertEqual(Chem.MolToSmiles(rwmol, True), Chem.MolToSmiles(rwmol.GetMol()))

  def test94CopyWithConfs(self):
    """ test copying Mols with some conformers

    """
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'cmpd2.tpl')
    m1 = Chem.MolFromTPLFile(fileN)
    self.assertTrue(m1 is not None)
    self.assertEqual(m1.GetNumAtoms(), 12)
    self.assertEqual(m1.GetNumConformers(), 2)
    self.assertEqual(m1.GetConformer(0).GetNumAtoms(), 12)
    self.assertEqual(m1.GetConformer(1).GetNumAtoms(), 12)

    m2 = Chem.Mol(m1)
    self.assertEqual(m2.GetNumAtoms(), 12)
    self.assertEqual(m2.GetNumConformers(), 2)
    self.assertEqual(m2.GetConformer(0).GetNumAtoms(), 12)
    self.assertEqual(m2.GetConformer(1).GetNumAtoms(), 12)

    m2 = Chem.Mol(m1, False, 0)
    self.assertEqual(m2.GetNumAtoms(), 12)
    self.assertEqual(m2.GetNumConformers(), 1)
    self.assertEqual(m2.GetConformer(0).GetNumAtoms(), 12)

    m2 = Chem.Mol(m1, False, 1)
    self.assertEqual(m2.GetNumAtoms(), 12)
    self.assertEqual(m2.GetNumConformers(), 1)
    self.assertEqual(m2.GetConformer(1).GetNumAtoms(), 12)

    m2 = Chem.Mol(m1, True)
    self.assertTrue(m2.GetNumAtoms() == 12)
    self.assertTrue(m2.GetNumConformers() == 0)

    m2 = Chem.RWMol(m1)
    self.assertEqual(m2.GetNumAtoms(), 12)
    self.assertEqual(m2.GetNumConformers(), 2)
    self.assertEqual(m2.GetConformer(0).GetNumAtoms(), 12)
    self.assertEqual(m2.GetConformer(1).GetNumAtoms(), 12)

    m2 = Chem.RWMol(m1, False, 0)
    self.assertEqual(m2.GetNumAtoms(), 12)
    self.assertEqual(m2.GetNumConformers(), 1)
    self.assertEqual(m2.GetConformer(0).GetNumAtoms(), 12)

    m2 = Chem.RWMol(m1, False, 1)
    self.assertEqual(m2.GetNumAtoms(), 12)
    self.assertEqual(m2.GetNumConformers(), 1)
    self.assertEqual(m2.GetConformer(1).GetNumAtoms(), 12)

    m2 = Chem.RWMol(m1, True)
    self.assertTrue(m2.GetNumAtoms() == 12)
    self.assertTrue(m2.GetNumConformers() == 0)

  def testAtomPropQueries(self):
    """ test the property queries
    """
    from rdkit.Chem import rdqueries

    m = Chem.MolFromSmiles("C" * 14)
    atoms = m.GetAtoms()
    atoms[0].SetProp("hah", "hah")
    atoms[1].SetIntProp("bar", 1)
    atoms[2].SetIntProp("bar", 2)
    atoms[3].SetBoolProp("baz", True)
    atoms[4].SetBoolProp("baz", False)
    atoms[5].SetProp("boo", "hoo")
    atoms[6].SetProp("boo", "-urns")
    atoms[7].SetDoubleProp("boot", 1.0)
    atoms[8].SetDoubleProp("boot", 4.0)
    atoms[9].SetDoubleProp("number", 4.0)
    atoms[10].SetIntProp("number", 4)

    tests = ((rdqueries.HasIntPropWithValueQueryAtom, "bar", {
      1: [1],
      2: [2]
    }), (rdqueries.HasBoolPropWithValueQueryAtom, "baz", {
      True: [3],
      False: [4]
    }), (rdqueries.HasStringPropWithValueQueryAtom, "boo", {
      "hoo": [5],
      "-urns": [6]
    }), (rdqueries.HasDoublePropWithValueQueryAtom, "boot", {
      1.0: [7],
      4.0: [8]
    }))

    for query, name, lookups in tests:
      for t, v in lookups.items():
        q = query(name, t)
        self.assertEqual(v, [x.GetIdx() for x in m.GetAtomsMatchingQuery(q)])
        q = query(name, t, negate=True)
        self.assertEqual(sorted(set(range(14)) - set(v)),
                         [x.GetIdx() for x in m.GetAtomsMatchingQuery(q)])

    # check tolerances
    self.assertEqual([
      x.GetIdx() for x in m.GetAtomsMatchingQuery(
        rdqueries.HasDoublePropWithValueQueryAtom("boot", 1.0, tolerance=3.))
    ], [7, 8])

    # numbers are numbers?, i.e. int!=double
    self.assertEqual([
      x.GetIdx()
      for x in m.GetAtomsMatchingQuery(rdqueries.HasIntPropWithValueQueryAtom("number", 4))
    ], [10])

  def testBondPropQueries(self):
    """ test the property queries
    """
    from rdkit.Chem import rdqueries

    m = Chem.MolFromSmiles("C" * 14)
    bonds = m.GetBonds()
    bonds[0].SetProp("hah", "hah")
    bonds[1].SetIntProp("bar", 1)
    bonds[2].SetIntProp("bar", 2)
    bonds[3].SetBoolProp("baz", True)
    bonds[4].SetBoolProp("baz", False)
    bonds[5].SetProp("boo", "hoo")
    bonds[6].SetProp("boo", "-urns")
    bonds[7].SetDoubleProp("boot", 1.0)
    bonds[8].SetDoubleProp("boot", 4.0)
    bonds[9].SetDoubleProp("number", 4.0)
    bonds[10].SetIntProp("number", 4)

    tests = ((rdqueries.HasIntPropWithValueQueryBond, "bar", {
      1: [1],
      2: [2]
    }), (rdqueries.HasBoolPropWithValueQueryBond, "baz", {
      True: [3],
      False: [4]
    }), (rdqueries.HasStringPropWithValueQueryBond, "boo", {
      "hoo": [5],
      "-urns": [6]
    }), (rdqueries.HasDoublePropWithValueQueryBond, "boot", {
      1.0: [7],
      4.0: [8]
    }))

    for query, name, lookups in tests:
      for t, v in lookups.items():
        q = query(name, t)
        self.assertEqual(v, [x.GetIdx() for x in m.GetBonds() if q.Match(x)])
        q = query(name, t, negate=True)
        self.assertEqual(sorted(set(range(13)) - set(v)),
                         [x.GetIdx() for x in m.GetBonds() if q.Match(x)])

    # check tolerances
    q = rdqueries.HasDoublePropWithValueQueryBond("boot", 1.0, tolerance=3.)
    self.assertEqual([x.GetIdx() for x in m.GetBonds() if q.Match(x)], [7, 8])

    # numbers are numbers?, i.e. int!=double
    q = rdqueries.HasIntPropWithValueQueryBond("number", 4)
    self.assertEqual([x.GetIdx() for x in m.GetBonds() if q.Match(x)], [10])

  def testGetShortestPath(self):
    """ test the GetShortestPath() wrapper
    """
    smi = "CC(OC1C(CCCC3)C3C(CCCC2)C2C1OC(C)=O)=O"
    m = Chem.MolFromSmiles(smi)
    path = Chem.GetShortestPath(m, 1, 20)
    self.assertEqual(path, (1, 2, 3, 16, 17, 18, 20))

  def testGithub497(self):
    with tempfile.TemporaryFile() as tmp, gzip.open(tmp) as outf:
      with self.assertRaises(ValueError):
        w = Chem.SDWriter(outf)

  def testGithub498(self):
    if (sys.version_info < (3, 0)):
      mode = 'w+'
    else:
      mode = 'wt+'
    m = Chem.MolFromSmiles('C')
    with tempfile.NamedTemporaryFile() as tmp, gzip.open(tmp, mode) as outf:
      w = Chem.SDWriter(outf)
      w.write(m)
      w.close()

  def testReplaceBond(self):
    origmol = Chem.RWMol(Chem.MolFromSmiles("CC"))
    bonds = list(origmol.GetBonds())
    self.assertEqual(len(bonds), 1)
    singlebond = bonds[0]
    self.assertEqual(singlebond.GetBondType(), Chem.BondType.SINGLE)

    # this is the only way we create a bond, is take it from another molecule
    doublebonded = Chem.MolFromSmiles("C=C")
    doublebond = list(doublebonded.GetBonds())[0]

    # make sure replacing the bond changes the smiles
    self.assertEqual(Chem.MolToSmiles(origmol), "CC")
    origmol.ReplaceBond(singlebond.GetIdx(), doublebond)
    Chem.SanitizeMol(origmol)

    self.assertEqual(Chem.MolToSmiles(origmol), "C=C")

  def testAdjustQueryProperties(self):
    m = Chem.MolFromSmarts('C1CCC1*')
    am = Chem.AdjustQueryProperties(m)
    self.assertTrue(Chem.MolFromSmiles('C1CCC1C').HasSubstructMatch(m))
    self.assertTrue(Chem.MolFromSmiles('C1CCC1C').HasSubstructMatch(am))
    self.assertTrue(Chem.MolFromSmiles('C1CC(C)C1C').HasSubstructMatch(m))
    self.assertFalse(Chem.MolFromSmiles('C1CC(C)C1C').HasSubstructMatch(am))

    m = Chem.MolFromSmiles('C1CCC1*')
    am = Chem.AdjustQueryProperties(m)
    self.assertFalse(Chem.MolFromSmiles('C1CCC1C').HasSubstructMatch(m))
    self.assertTrue(Chem.MolFromSmiles('C1CCC1C').HasSubstructMatch(am))
    qps = Chem.AdjustQueryParameters()
    qps.makeDummiesQueries = False
    am = Chem.AdjustQueryProperties(m, qps)
    self.assertFalse(Chem.MolFromSmiles('C1CCC1C').HasSubstructMatch(am))

    m = Chem.MolFromSmiles('C1=CC=CC=C1', sanitize=False)
    am = Chem.AdjustQueryProperties(m)
    self.assertTrue(Chem.MolFromSmiles('c1ccccc1').HasSubstructMatch(am))
    qp = Chem.AdjustQueryParameters()
    qp.aromatizeIfPossible = False
    am = Chem.AdjustQueryProperties(m, qp)
    self.assertFalse(Chem.MolFromSmiles('c1ccccc1').HasSubstructMatch(am))

    m = Chem.MolFromSmiles('C1CCC1OC')
    qps = Chem.AdjustQueryParameters()
    qps.makeAtomsGeneric = True
    am = Chem.AdjustQueryProperties(m, qps)
    self.assertEqual(Chem.MolToSmarts(am), '*1-*-*-*-1-*-*')
    qps.makeAtomsGenericFlags = Chem.ADJUST_IGNORERINGS
    am = Chem.AdjustQueryProperties(m, qps)
    self.assertEqual(Chem.MolToSmarts(am), '[#6&D2]1-[#6&D2]-[#6&D2]-[#6&D3]-1-*-*')

    qps = Chem.AdjustQueryParameters()
    qps.makeBondsGeneric = True
    am = Chem.AdjustQueryProperties(m, qps)
    self.assertEqual(Chem.MolToSmarts(am), '[#6&D2]1~[#6&D2]~[#6&D2]~[#6&D3]~1~[#8]~[#6]')
    qps.makeBondsGenericFlags = Chem.ADJUST_IGNORERINGS
    am = Chem.AdjustQueryProperties(m, qps)
    self.assertEqual(Chem.MolToSmarts(am), '[#6&D2]1-[#6&D2]-[#6&D2]-[#6&D3]-1~[#8]~[#6]')

  def testMolFragmentSmarts(self):
    m = Chem.MolFromSmiles('C1CCC1OC')
    self.assertEqual(Chem.MolFragmentToSmarts(m, [0, 1, 2]), '[#6]-[#6]-[#6]')
    # if bondsToUse is honored, the ring won't show up
    self.assertEqual(Chem.MolFragmentToSmarts(m, [0, 1, 2, 3], bondsToUse=[0, 1, 2, 3]),
                     '[#6]-[#6]-[#6]-[#6]')

    # Does MolFragmentToSmarts accept output of AdjustQueryProperties?
    qps = Chem.AdjustQueryParameters()
    qps.makeAtomsGeneric = True
    am = Chem.AdjustQueryProperties(m, qps)
    self.assertEqual(Chem.MolFragmentToSmarts(am, [0, 1, 2]), '*-*-*')

  def testAdjustQueryPropertiesgithubIssue1474(self):
    core = Chem.MolFromSmiles('[*:1]C1N([*:2])C([*:3])O1')
    core.GetAtomWithIdx(0).SetProp('foo', 'bar')
    core.GetAtomWithIdx(1).SetProp('foo', 'bar')

    ap = Chem.AdjustQueryProperties(core)
    self.assertEqual(ap.GetAtomWithIdx(0).GetPropsAsDict()["foo"], "bar")
    self.assertEqual(ap.GetAtomWithIdx(1).GetPropsAsDict()["foo"], "bar")

  def testReplaceAtomWithQueryAtom(self):
    mol = Chem.MolFromSmiles("CC(C)C")
    qmol = Chem.MolFromSmiles("C")
    matches = mol.GetSubstructMatches(qmol)
    self.assertEqual(((0,), (1,), (2,), (3,)), matches)

    atom = qmol.GetAtomWithIdx(0)
    natom = rdqueries.ReplaceAtomWithQueryAtom(qmol, atom)
    qa = rdqueries.ExplicitDegreeEqualsQueryAtom(3)
    natom.ExpandQuery(qa, Chem.CompositeQueryType.COMPOSITE_AND)
    matches = mol.GetSubstructMatches(qmol)
    self.assertEqual(((1,),), matches)

  def testGithubIssue579(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'NCI_aids_few.sdf.gz')
    inf = gzip.open(fileN)
    suppl = Chem.ForwardSDMolSupplier(inf)
    m0 = next(suppl)
    self.assertIsNot(m0, None)
    inf.close()
    del suppl

  def testSequenceBasics(self):
    " very basic round-tripping of the sequence reader/writer support "
    helm = 'PEPTIDE1{C.Y.I.Q.N.C.P.L.G}$$$$'
    seq = 'CYIQNCPLG'
    fasta = '>\nCYIQNCPLG\n'
    smi = 'CC[C@H](C)[C@H](NC(=O)[C@H](Cc1ccc(O)cc1)NC(=O)[C@@H](N)CS)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H](CS)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CC(C)C)C(=O)NCC(=O)O'

    m = Chem.MolFromSequence(seq)
    self.assertTrue(m is not None)
    self.assertEqual(Chem.MolToSequence(m), seq)
    self.assertEqual(Chem.MolToHELM(m), helm)
    self.assertEqual(Chem.MolToFASTA(m), fasta)
    self.assertEqual(Chem.MolToSmiles(m, isomericSmiles=True), smi)

    m = Chem.MolFromHELM(helm)
    self.assertTrue(m is not None)
    self.assertEqual(Chem.MolToSequence(m), seq)
    self.assertEqual(Chem.MolToHELM(m), helm)
    self.assertEqual(Chem.MolToFASTA(m), fasta)
    self.assertEqual(Chem.MolToSmiles(m, isomericSmiles=True), smi)

    m = Chem.MolFromFASTA(fasta)
    self.assertTrue(m is not None)
    self.assertEqual(Chem.MolToSequence(m), seq)
    self.assertEqual(Chem.MolToHELM(m), helm)
    self.assertEqual(Chem.MolToFASTA(m), fasta)
    self.assertEqual(Chem.MolToSmiles(m, isomericSmiles=True), smi)

    seq = "CGCGAATTACCGCG"
    m = Chem.MolFromSequence(seq, flavor=6)  # DNA
    self.assertEqual(Chem.MolToSequence(m), 'CGCGAATTACCGCG')
    self.assertEqual(
      Chem.MolToHELM(m),
      'RNA1{[dR](C)P.[dR](G)P.[dR](C)P.[dR](G)P.[dR](A)P.[dR](A)P.[dR](T)P.[dR](T)P.[dR](A)P.[dR](C)P.[dR](C)P.[dR](G)P.[dR](C)P.[dR](G)}$$$$'
    )
    seq = "CGCGAAUUACCGCG"
    m = Chem.MolFromSequence(seq, flavor=2)  # RNA
    self.assertEqual(Chem.MolToSequence(m), 'CGCGAAUUACCGCG')
    self.assertEqual(
      Chem.MolToHELM(m),
      'RNA1{R(C)P.R(G)P.R(C)P.R(G)P.R(A)P.R(A)P.R(U)P.R(U)P.R(A)P.R(C)P.R(C)P.R(G)P.R(C)P.R(G)}$$$$'
    )
    m = Chem.MolFromSequence(seq, flavor=3)  # RNA - 5' cap
    self.assertEqual(Chem.MolToSequence(m), 'CGCGAAUUACCGCG')
    self.assertEqual(
      Chem.MolToHELM(m),
      'RNA1{P.R(C)P.R(G)P.R(C)P.R(G)P.R(A)P.R(A)P.R(U)P.R(U)P.R(A)P.R(C)P.R(C)P.R(G)P.R(C)P.R(G)}$$$$'
    )

  def testResMolSupplier(self):
    mol = Chem.MolFromSmiles('CC')
    resMolSuppl = Chem.ResonanceMolSupplier(mol)
    del resMolSuppl
    resMolSuppl = Chem.ResonanceMolSupplier(mol)
    self.assertEqual(resMolSuppl.GetNumConjGrps(), 0)
    self.assertEqual(len(resMolSuppl), 1)
    self.assertEqual(resMolSuppl.GetNumConjGrps(), 0)

    mol = Chem.MolFromSmiles('NC(=[NH2+])c1ccc(cc1)C(=O)[O-]')
    totalFormalCharge = getTotalFormalCharge(mol)

    resMolSuppl = Chem.ResonanceMolSupplier(mol)
    self.assertFalse(resMolSuppl.GetIsEnumerated())
    self.assertEqual(len(resMolSuppl), 4)
    self.assertTrue(resMolSuppl.GetIsEnumerated())

    resMolSuppl = Chem.ResonanceMolSupplier(mol)
    self.assertFalse(resMolSuppl.GetIsEnumerated())
    resMolSuppl.Enumerate()
    self.assertTrue(resMolSuppl.GetIsEnumerated())
    self.assertTrue(
      (resMolSuppl[0].GetBondBetweenAtoms(0, 1).GetBondType() != resMolSuppl[1].GetBondBetweenAtoms(
        0, 1).GetBondType()) or (resMolSuppl[0].GetBondBetweenAtoms(9, 10).GetBondType()
                                 != resMolSuppl[1].GetBondBetweenAtoms(9, 10).GetBondType()))

    resMolSuppl = Chem.ResonanceMolSupplier(mol, Chem.KEKULE_ALL)
    self.assertEqual(len(resMolSuppl), 8)
    bondTypeSet = set()
    # check that we actually have two alternate Kekule structures
    bondTypeSet.add(resMolSuppl[0].GetBondBetweenAtoms(3, 4).GetBondType())
    bondTypeSet.add(resMolSuppl[1].GetBondBetweenAtoms(3, 4).GetBondType())
    self.assertEqual(len(bondTypeSet), 2)

    bondTypeDict = {}
    resMolSuppl = Chem.ResonanceMolSupplier(
      mol, Chem.ALLOW_INCOMPLETE_OCTETS
      | Chem.UNCONSTRAINED_CATIONS
      | Chem.UNCONSTRAINED_ANIONS)
    self.assertEqual(len(resMolSuppl), 32)
    for i in range(len(resMolSuppl)):
      resMol = resMolSuppl[i]
      self.assertEqual(getTotalFormalCharge(resMol), totalFormalCharge)
    while (not resMolSuppl.atEnd()):
      resMol = next(resMolSuppl)
      self.assertEqual(getTotalFormalCharge(resMol), totalFormalCharge)
    resMolSuppl.reset()
    cmpFormalChargeBondOrder(self, resMolSuppl[0], next(resMolSuppl))

    resMolSuppl = Chem.ResonanceMolSupplier(
      mol, Chem.ALLOW_INCOMPLETE_OCTETS
      | Chem.UNCONSTRAINED_CATIONS
      | Chem.UNCONSTRAINED_ANIONS, 10)
    self.assertEqual(len(resMolSuppl), 10)

    crambinPdb = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                              '1CRN.pdb')
    mol = Chem.MolFromPDBFile(crambinPdb)
    resMolSuppl = Chem.ResonanceMolSupplier(mol)
    self.assertEqual(len(resMolSuppl), 1)
    resMolSuppl = Chem.ResonanceMolSupplier(mol, Chem.KEKULE_ALL)
    self.assertEqual(len(resMolSuppl), 8)

  def testGithub4884(self):
    # test that we don't hang
    mol = Chem.MolFromSmiles('O=[N+][O-]')
    supp = Chem.ResonanceMolSupplier(mol)
    supp.atEnd()

  def testSubstructMatchAcetate(self):
    mol = Chem.MolFromSmiles('CC(=O)[O-]')
    query = Chem.MolFromSmarts('C(=O)[O-]')

    resMolSuppl = Chem.ResonanceMolSupplier(mol)
    matches = mol.GetSubstructMatches(query)
    self.assertEqual(len(matches), 1)
    self.assertEqual(matches, ((1, 2, 3), ))
    matches = mol.GetSubstructMatches(query, uniquify=True)
    self.assertEqual(len(matches), 1)
    self.assertEqual(matches, ((1, 2, 3), ))
    matches = mol.GetSubstructMatches(query, uniquify=False)
    self.assertEqual(len(matches), 1)
    self.assertEqual(matches, ((1, 2, 3), ))
    matches = resMolSuppl.GetSubstructMatches(query)
    self.assertEqual(len(matches), 2)
    self.assertEqual(matches, ((1, 2, 3), (1, 3, 2)))
    matches = resMolSuppl.GetSubstructMatches(query, uniquify=True)
    self.assertEqual(len(matches), 1)
    self.assertEqual(matches, ((1, 2, 3), ))
    matches = resMolSuppl.GetSubstructMatches(query, uniquify=False)
    self.assertEqual(len(matches), 2)
    self.assertEqual(matches, ((1, 2, 3), (1, 3, 2)))
    query = Chem.MolFromSmarts('C(~O)~O')
    matches = mol.GetSubstructMatches(query, uniquify=False)
    self.assertEqual(len(matches), 2)
    self.assertEqual(matches, ((1, 2, 3), (1, 3, 2)))
    matches = mol.GetSubstructMatches(query, uniquify=True)
    self.assertEqual(len(matches), 1)
    self.assertEqual(matches, ((1, 2, 3), ))
    matches = resMolSuppl.GetSubstructMatches(query, uniquify=False)
    self.assertEqual(len(matches), 2)
    self.assertEqual(matches, ((1, 2, 3), (1, 3, 2)))
    matches = resMolSuppl.GetSubstructMatches(query, uniquify=True)
    self.assertEqual(len(matches), 1)
    self.assertEqual(matches, ((1, 2, 3), ))

  def testSubstructMatchDMAP(self):
    mol = Chem.MolFromSmiles('C(C)Nc1cc[nH+]cc1')
    query = Chem.MolFromSmarts('[#7+]')

    resMolSuppl = Chem.ResonanceMolSupplier(mol)
    matches = mol.GetSubstructMatches(query, False, False, False)
    self.assertEqual(len(matches), 1)
    p = matches[0]
    self.assertEqual(p[0], 6)
    matches = resMolSuppl.GetSubstructMatches(query, False, False, False)
    self.assertEqual(len(matches), 2)
    v = []
    p = matches[0]
    v.append(p[0])
    p = matches[1]
    v.append(p[0])
    v.sort()
    self.assertEqual(v[0], 2)
    self.assertEqual(v[1], 6)

  def testCrambin(self):
    crambinPdb = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                              '1CRN.pdb')
    crambin = Chem.MolFromPDBFile(crambinPdb)
    res = []
    # protonate NH2
    res.append(Chem.MolFromSmarts('[Nh2][Ch;Ch2]'))
    # protonate Arg
    res.append(Chem.MolFromSmarts('[Nh][C]([Nh2])=[Nh]'))
    setResidueFormalCharge(crambin, res, 1)
    res = []
    # deprotonate COOH
    res.append(Chem.MolFromSmarts('C(=O)[Oh]'))
    setResidueFormalCharge(crambin, res, -1)
    res = []
    resMolSupplST = Chem.ResonanceMolSupplier(crambin)
    # crambin has 2 Arg (3 resonance structures each); 1 Asp, 1 Glu
    # and 1 terminal COO- (2 resonance structures each)
    # so possible resonance structures are 3^2 * 2^3 = 72
    self.assertEqual(len(resMolSupplST), 72)
    self.assertEqual(resMolSupplST.GetNumConjGrps(), 56)
    carboxylateQuery = Chem.MolFromSmarts('C(=O)[O-]')
    guanidiniumQuery = Chem.MolFromSmarts('NC(=[NH2+])N')
    matches = crambin.GetSubstructMatches(carboxylateQuery)
    self.assertEqual(len(matches), 3)
    matches = crambin.GetSubstructMatches(carboxylateQuery, uniquify=False)
    self.assertEqual(len(matches), 3)
    matches = crambin.GetSubstructMatches(guanidiniumQuery)
    self.assertEqual(len(matches), 0)
    matches = crambin.GetSubstructMatches(guanidiniumQuery, uniquify=False)
    self.assertEqual(len(matches), 0)
    matches = resMolSupplST.GetSubstructMatches(carboxylateQuery)
    self.assertEqual(len(matches), 6)
    self.assertEqual(matches, ((166, 167, 168), (166, 168, 167), (298, 299, 300), (298, 300, 299),
                               (320, 321, 326), (320, 326, 321)))
    matches = resMolSupplST.GetSubstructMatches(carboxylateQuery, uniquify=True)
    self.assertEqual(len(matches), 3)
    self.assertEqual(matches, ((166, 167, 168), (298, 299, 300), (320, 321, 326)))
    matches = resMolSupplST.GetSubstructMatches(guanidiniumQuery)
    self.assertEqual(len(matches), 8)
    self.assertEqual(matches, ((66, 67, 68, 69), (66, 67, 69, 68), (68, 67, 69, 66),
                               (69, 67, 68, 66), (123, 124, 125, 126), (123, 124, 126, 125),
                               (125, 124, 126, 123), (126, 124, 125, 123)))
    matches = resMolSupplST.GetSubstructMatches(guanidiniumQuery, uniquify=True)
    self.assertEqual(len(matches), 2)
    self.assertEqual(matches, ((66, 67, 68, 69), (123, 124, 125, 126)))
    btList2ST = getBtList2(resMolSupplST)
    self.assertTrue(btList2ST)
    resMolSupplMT = Chem.ResonanceMolSupplier(crambin)
    resMolSupplMT.SetNumThreads(0)
    self.assertEqual(len(resMolSupplST), len(resMolSupplMT))
    btList2MT = getBtList2(resMolSupplMT)
    self.assertTrue(btList2MT)
    self.assertEqual(len(btList2ST), len(btList2MT))
    for i in range(len(btList2ST)):
      for j in range(len(btList2ST)):
        self.assertEqual(btList2ST[i][j], btList2MT[i][j])
    for suppl in [resMolSupplST, resMolSupplMT]:
      matches = suppl.GetSubstructMatches(carboxylateQuery, numThreads=0)
      self.assertEqual(len(matches), 6)
      self.assertEqual(matches, ((166, 167, 168), (166, 168, 167), (298, 299, 300), (298, 300, 299),
                                 (320, 321, 326), (320, 326, 321)))
      matches = suppl.GetSubstructMatches(carboxylateQuery, uniquify=True, numThreads=0)
      self.assertEqual(len(matches), 3)
      self.assertEqual(matches, ((166, 167, 168), (298, 299, 300), (320, 321, 326)))
      matches = suppl.GetSubstructMatches(guanidiniumQuery, numThreads=0)
      self.assertEqual(len(matches), 8)
      self.assertEqual(matches, ((66, 67, 68, 69), (66, 67, 69, 68), (68, 67, 69, 66),
                                 (69, 67, 68, 66), (123, 124, 125, 126), (123, 124, 126, 125),
                                 (125, 124, 126, 123), (126, 124, 125, 123)))
      matches = suppl.GetSubstructMatches(guanidiniumQuery, uniquify=True, numThreads=0)
      self.assertEqual(len(matches), 2)
      self.assertEqual(matches, ((66, 67, 68, 69), (123, 124, 125, 126)))

  def testGitHub1166(self):
    mol = Chem.MolFromSmiles('NC(=[NH2+])c1ccc(cc1)C(=O)[O-]')
    resMolSuppl = Chem.ResonanceMolSupplier(mol, Chem.KEKULE_ALL)
    self.assertEqual(len(resMolSuppl), 8)
    # check that formal charges on odd indices are in the same position
    # as on even indices
    for i in range(0, len(resMolSuppl), 2):
      self.assertEqual(resMolSuppl[i].GetNumAtoms(), resMolSuppl[i + 1].GetNumAtoms())
      for atomIdx in range(resMolSuppl[i].GetNumAtoms()):
        self.assertEqual(resMolSuppl[i].GetAtomWithIdx(atomIdx).GetFormalCharge(),
                         resMolSuppl[i + 1].GetAtomWithIdx(atomIdx).GetFormalCharge())
      # check that bond orders are alternate on aromatic bonds between
      # structures on odd indices and structures on even indices
      self.assertEqual(resMolSuppl[i].GetNumBonds(), resMolSuppl[i + 1].GetNumBonds())
      for bondIdx in range(resMolSuppl[i].GetNumBonds()):
        self.assertTrue(
          ((not resMolSuppl[i].GetBondWithIdx(bondIdx).GetIsAromatic()) and
           (not resMolSuppl[i + 1].GetBondWithIdx(bondIdx).GetIsAromatic()) and
           (resMolSuppl[i].GetBondWithIdx(bondIdx).GetBondType()
            == resMolSuppl[i + 1].GetBondWithIdx(bondIdx).GetBondType()))
          or (resMolSuppl[i].GetBondWithIdx(bondIdx).GetIsAromatic()
              and resMolSuppl[i + 1].GetBondWithIdx(bondIdx).GetIsAromatic() and (int(
                round(resMolSuppl[i].GetBondWithIdx(bondIdx).GetBondTypeAsDouble() +
                      resMolSuppl[i + 1].GetBondWithIdx(bondIdx).GetBondTypeAsDouble())) == 3)))

  def testConjGrpPerception(self):
    mol1 = Chem.MolFromMolBlock("""\

     RDKit          2D

 14 15  0  0  0  0  0  0  0  0999 V2000
    3.7539   -1.2744    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4317   -0.5660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1571   -1.3568    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1651   -0.6484    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4397   -1.4393    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3921   -2.9385    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7619   -0.7309    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8095    0.7684    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1316    1.4768    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5349    1.5592    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2127    0.8508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0619    1.6417    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    2.3841    0.9333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6587    1.7241    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  4  0
  3  4  4  0
  4  5  4  0
  5  6  1  0
  5  7  4  0
  7  8  4  0
  8  9  1  0
  8 10  4  0
 10 11  4  0
 11 12  4  0
 12 13  4  0
 13 14  1  0
 13  2  4  0
 11  4  4  0
M  END
$$$$
""")
    mol2 = Chem.MolFromMolBlock("""\

     RDKit          2D

 14 15  0  0  0  0  0  0  0  0999 V2000
    1.0619   -1.6417    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2127   -0.8508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5349   -1.5592    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8095   -0.7684    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7619    0.7309    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4397    1.4393    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1651    0.6484    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1571    1.3568    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4317    0.5660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7539    1.2744    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3841   -0.9333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6587   -1.7241    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1316   -1.4768    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3921    2.9385    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  4  0
  3  4  4  0
  4  5  4  0
  5  6  4  0
  2  3  4  0
  2  7  4  0
  7  8  4  0
  8  9  4  0
  9 10  1  0
  9 11  4  0
 11 12  1  0
 11  1  4  0
  6  7  4  0
  4 13  1  0
  6 14  1  0
M  END
$$$$
""")
    resMolSuppl1 = Chem.ResonanceMolSupplier(mol1, Chem.KEKULE_ALL)
    self.assertEqual(len(resMolSuppl1), 3)
    resMolSuppl2 = Chem.ResonanceMolSupplier(mol2, Chem.KEKULE_ALL)
    self.assertEqual(len(resMolSuppl2), 3)

  def testGitHub2597(self):

    class MyBrokenCallBack(Chem.ResonanceMolSupplier):

      def __call__(self):
        return True

    class MyBrokenCallBack2(Chem.ResonanceMolSupplierCallback):
      pass

    class ExceedNumStructures(Chem.ResonanceMolSupplierCallback):

      def __init__(self, parent):
        super().__init__()
        self._parent = parent

      def __call__(self):
        self._parent.assertEqual(self.GetNumConjGrps(), 1)
        return (self.GetNumStructures(0) < 12)

    class ExceedNumDiverseStructures(Chem.ResonanceMolSupplierCallback):

      def __init__(self, parent):
        super().__init__()
        self._parent = parent

      def __call__(self):
        self._parent.assertEqual(self.GetNumConjGrps(), 1)
        return (self.GetNumDiverseStructures(0) < 8)

    class ExceedTimeout(Chem.ResonanceMolSupplierCallback):

      def __init__(self, parent):
        super().__init__()
        self.start_time = None
        self.timeout = timedelta(seconds=3)
        self._parent = parent

      def __call__(self):
        if (self.start_time is None):
          self.start_time = datetime.now()
        return (datetime.now() - self.start_time < self.timeout)

    mol = Chem.MolFromSmiles(
      "ClC1=NC(NC2=CC=CC3=C2C(=O)C2=CC=CC=C2C3=O)=NC(NC2=CC=CC3=C2C(=O)C2=CC=CC=C2C3=O)=N1")
    resMolSuppl = Chem.ResonanceMolSupplier(mol)
    self.assertEqual(len(resMolSuppl), 1)
    resMolSuppl = Chem.ResonanceMolSupplier(mol, Chem.KEKULE_ALL)

    self.assertEqual(len(resMolSuppl), 32)
    resMolSuppl = Chem.ResonanceMolSupplier(mol, Chem.ALLOW_CHARGE_SEPARATION, 10)
    self.assertEqual(len(resMolSuppl), 10)
    self.assertFalse(resMolSuppl.WasCanceled())
    resMolSuppl = Chem.ResonanceMolSupplier(mol, Chem.ALLOW_CHARGE_SEPARATION)
    callback = resMolSuppl.GetProgressCallback()
    self.assertIsNone(callback)
    resMolSuppl.SetProgressCallback(ExceedNumStructures(self))
    callback = resMolSuppl.GetProgressCallback()
    self.assertTrue(isinstance(callback, ExceedNumStructures))
    resMolSuppl.SetProgressCallback(None)
    callback = resMolSuppl.GetProgressCallback()
    self.assertIsNone(callback)
    resMolSuppl.SetProgressCallback(ExceedNumStructures(self))
    self.assertEqual(len(resMolSuppl), 12)
    self.assertTrue(resMolSuppl.WasCanceled())
    resMolSuppl = Chem.ResonanceMolSupplier(mol, Chem.ALLOW_CHARGE_SEPARATION)
    with self.assertRaises(TypeError):
      resMolSuppl.SetProgressCallback(MyBrokenCallBack())
    with self.assertRaises(AttributeError):
      resMolSuppl.SetProgressCallback(MyBrokenCallBack2())
    resMolSuppl.SetProgressCallback(ExceedNumDiverseStructures(self))
    self.assertEqual(len(resMolSuppl), 9)
    self.assertTrue(resMolSuppl.WasCanceled())
    resMolSuppl = Chem.ResonanceMolSupplier(
      mol, Chem.UNCONSTRAINED_CATIONS | Chem.UNCONSTRAINED_ANIONS | Chem.KEKULE_ALL)
    resMolSuppl.SetProgressCallback(ExceedTimeout(self))
    resMolSuppl.Enumerate()
    print(len(resMolSuppl))
    self.assertTrue(resMolSuppl.WasCanceled())

  def testAtomBondProps(self):
    m = Chem.MolFromSmiles('c1ccccc1')
    for atom in m.GetAtoms():
      d = atom.GetPropsAsDict()
      self.assertEqual(set(d.keys()), set(['_CIPRank', '__computedProps']))
      self.assertEqual(d['_CIPRank'], 0)
      self.assertEqual(list(d['__computedProps']), ['_CIPRank'])

    for bond in m.GetBonds():
      self.assertEqual(bond.GetPropsAsDict(), {})

  def testSDProps(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'NCI_aids_few.sdf')
    #fileN = "../FileParsers/test_data/NCI_aids_few.sdf"
    sddata = [
      {
        '_MolFileInfo': 'BBtclserve11129916382D 0   0.00000     0.00000    48',
        'NSC': 48,
        'NCI_AIDS_Antiviral_Screen_IC50': '2.00E-04\tM\t=\t2.46E-05\t3',
        '_Name': 48,
        'CAS_RN': '15716-70-8',
        '_MolFileChiralFlag': 0,
        '_MolFileComments': '15716-70-8',
        'NCI_AIDS_Antiviral_Screen_EC50': '2.00E-04\tM\t>\t2.00E-04\t3',
        'NCI_AIDS_Antiviral_Screen_Conclusion': 'CI'
      },
      {
        '_MolFileInfo': 'BBtclserve11129916382D 0   0.00000     0.00000    78',
        'NSC': 78,
        'NCI_AIDS_Antiviral_Screen_IC50': '2.00E-04\tM\t=\t9.80E-05\t3',
        '_Name': 78,
        'CAS_RN': '6290-84-2',
        '_MolFileChiralFlag': 0,
        '_MolFileComments': '6290-84-2',
        'NCI_AIDS_Antiviral_Screen_EC50': '2.00E-04\tM\t>\t2.00E-04\t3',
        'NCI_AIDS_Antiviral_Screen_Conclusion': 'CI'
      },
      {
        '_MolFileInfo': 'BBtclserve11129916382D 0   0.00000     0.00000   128',
        'NSC': 128,
        'NCI_AIDS_Antiviral_Screen_IC50': '2.00E-04\tM\t=\t4.60E-05\t4',
        '_Name': 128,
        'CAS_RN': '5395-10-8',
        '_MolFileChiralFlag': 0,
        '_MolFileComments': '5395-10-8',
        'NCI_AIDS_Antiviral_Screen_EC50': '2.00E-04\tM\t>\t2.00E-04\t4',
        'NCI_AIDS_Antiviral_Screen_Conclusion': 'CI'
      },
      {
        '_MolFileInfo': 'BBtclserve11129916382D 0   0.00000     0.00000   163',
        'NSC': 163,
        'NCI_AIDS_Antiviral_Screen_IC50': '6.75E-04\tM\t>\t6.75E-04\t2',
        '_Name': 163,
        'CAS_RN': '81-11-8',
        '_MolFileChiralFlag': 0,
        '_MolFileComments': '81-11-8',
        'NCI_AIDS_Antiviral_Screen_EC50': '6.75E-04\tM\t>\t6.75E-04\t2',
        'NCI_AIDS_Antiviral_Screen_Conclusion': 'CI'
      },
      {
        '_MolFileInfo': 'BBtclserve11129916382D 0   0.00000     0.00000   164',
        'NSC': 164,
        'NCI_AIDS_Antiviral_Screen_IC50': '2.00E-04\tM\t>\t2.00E-04\t2',
        '_Name': 164,
        'CAS_RN': '5325-43-9',
        '_MolFileChiralFlag': 0,
        '_MolFileComments': '5325-43-9',
        'NCI_AIDS_Antiviral_Screen_EC50': '2.00E-04\tM\t>\t2.00E-04\t2',
        'NCI_AIDS_Antiviral_Screen_Conclusion': 'CI'
      },
      {
        '_MolFileInfo': 'BBtclserve11129916382D 0   0.00000     0.00000   170',
        'NSC': 170,
        '_Name': 170,
        'CAS_RN': '999-99-9',
        '_MolFileChiralFlag': 0,
        '_MolFileComments': '999-99-9',
        'NCI_AIDS_Antiviral_Screen_EC50': '9.47E-04\tM\t>\t9.47E-04\t1',
        'NCI_AIDS_Antiviral_Screen_Conclusion': 'CI'
      },
      {
        '_MolFileInfo': 'BBtclserve11129916382D 0   0.00000     0.00000   180',
        'NSC': 180,
        'NCI_AIDS_Antiviral_Screen_IC50':
        '6.46E-04\tM\t=\t5.80E-04\t2\n1.81E-03\tM\t=\t6.90E-04\t2',
        '_Name': 180,
        'CAS_RN': '69-72-7',
        '_MolFileChiralFlag': 0,
        '_MolFileComments': '69-72-7',
        'NCI_AIDS_Antiviral_Screen_EC50':
        '6.46E-04\tM\t>\t6.46E-04\t2\n1.81E-03\tM\t>\t1.81E-03\t2',
        'NCI_AIDS_Antiviral_Screen_Conclusion': 'CI'
      },
      {
        '_MolFileInfo': 'BBtclserve11129916382D 0   0.00000     0.00000   186',
        'NSC': 186,
        'NCI_AIDS_Antiviral_Screen_IC50': '1.44E-04\tM\t=\t2.49E-05\t2',
        '_Name': 186,
        'CAS_RN': '518-75-2',
        '_MolFileChiralFlag': 0,
        '_MolFileComments': '518-75-2',
        'NCI_AIDS_Antiviral_Screen_EC50': '1.44E-04\tM\t>\t1.44E-04\t2',
        'NCI_AIDS_Antiviral_Screen_Conclusion': 'CI'
      },
      {
        '_MolFileInfo': 'BBtclserve11129916382D 0   0.00000     0.00000   192',
        'NSC': 192,
        'NCI_AIDS_Antiviral_Screen_IC50': '2.00E-04\tM\t=\t3.38E-06\t2',
        '_Name': 192,
        'CAS_RN': '2217-55-2',
        '_MolFileChiralFlag': 0,
        '_MolFileComments': '2217-55-2',
        'NCI_AIDS_Antiviral_Screen_EC50': '2.00E-04\tM\t>\t2.00E-04\t2',
        'NCI_AIDS_Antiviral_Screen_Conclusion': 'CI'
      },
      {
        '_MolFileInfo': 'BBtclserve11129916382D 0   0.00000     0.00000   203',
        'NSC': 203,
        '_Name': 203,
        'CAS_RN': '1155-00-6',
        '_MolFileChiralFlag': 0,
        '_MolFileComments': '1155-00-6',
        'NCI_AIDS_Antiviral_Screen_Conclusion': 'CI'
      },
      {
        '_MolFileInfo': 'BBtclserve11129916382D 0   0.00000     0.00000   210',
        'NSC': 210,
        'NCI_AIDS_Antiviral_Screen_IC50': '1.33E-03\tM\t>\t1.33E-03\t2',
        '_Name': 210,
        'CAS_RN': '5325-75-7',
        '_MolFileChiralFlag': 0,
        '_MolFileComments': '5325-75-7',
        'NCI_AIDS_Antiviral_Screen_EC50': '1.33E-03\tM\t>\t1.33E-03\t2',
        'NCI_AIDS_Antiviral_Screen_Conclusion': 'CI'
      },
      {
        '_MolFileInfo': 'BBtclserve11129916382D 0   0.00000     0.00000   211',
        'NSC': 211,
        'NCI_AIDS_Antiviral_Screen_IC50':
        '2.00E-04\tM\t>\t2.00E-04\t8\n2.00E-03\tM\t=\t1.12E-03\t2',
        '_Name': 211,
        'CAS_RN': '5325-76-8',
        '_MolFileChiralFlag': 0,
        '_MolFileComments': '5325-76-8',
        'NCI_AIDS_Antiviral_Screen_EC50':
        '2.00E-04\tM\t>\t7.42E-05\t8\n2.00E-03\tM\t=\t6.35E-05\t2',
        'NCI_AIDS_Antiviral_Screen_Conclusion': 'CM'
      },
      {
        '_MolFileInfo': 'BBtclserve11129916382D 0   0.00000     0.00000   213',
        'NSC': 213,
        'NCI_AIDS_Antiviral_Screen_IC50': '2.00E-04\tM\t>\t2.00E-04\t4',
        '_Name': 213,
        'CAS_RN': '119-80-2',
        '_MolFileChiralFlag': 0,
        '_MolFileComments': '119-80-2',
        'NCI_AIDS_Antiviral_Screen_EC50': '2.00E-04\tM\t>\t2.00E-04\t4',
        'NCI_AIDS_Antiviral_Screen_Conclusion': 'CI'
      },
      {
        '_MolFileInfo': 'BBtclserve11129916382D 0   0.00000     0.00000   220',
        'NSC': 220,
        'NCI_AIDS_Antiviral_Screen_IC50': '2.00E-04\tM\t>\t2.00E-04\t4',
        '_Name': 220,
        'CAS_RN': '5325-83-7',
        '_MolFileChiralFlag': 0,
        '_MolFileComments': '5325-83-7',
        'NCI_AIDS_Antiviral_Screen_EC50': '2.00E-04\tM\t>\t2.00E-04\t4',
        'NCI_AIDS_Antiviral_Screen_Conclusion': 'CI'
      },
      {
        '_MolFileInfo': 'BBtclserve11129916382D 0   0.00000     0.00000   229',
        'NSC': 229,
        'NCI_AIDS_Antiviral_Screen_IC50': '2.00E-04\tM\t>\t2.00E-04\t2',
        '_Name': 229,
        'CAS_RN': '5325-88-2',
        '_MolFileChiralFlag': 0,
        '_MolFileComments': '5325-88-2',
        'NCI_AIDS_Antiviral_Screen_EC50': '2.00E-04\tM\t>\t2.00E-04\t2',
        'NCI_AIDS_Antiviral_Screen_Conclusion': 'CI'
      },
      {
        '_MolFileInfo': 'BBtclserve11129916382D 0   0.00000     0.00000   256',
        'NSC': 256,
        'NCI_AIDS_Antiviral_Screen_IC50': '2.00E-04\tM\t>\t2.00E-04\t4',
        '_Name': 256,
        'CAS_RN': '5326-06-7',
        '_MolFileChiralFlag': 0,
        '_MolFileComments': '5326-06-7',
        'NCI_AIDS_Antiviral_Screen_EC50': '2.00E-04\tM\t>\t2.00E-04\t4',
        'NCI_AIDS_Antiviral_Screen_Conclusion': 'CI'
      },
    ]
    sdSup = Chem.SDMolSupplier(fileN)
    for i, mol in enumerate(sdSup):
      self.assertEqual(mol.GetPropsAsDict(includePrivate=True), sddata[i])

  def testGetSetProps(self):
    m = Chem.MolFromSmiles("CC")
    errors = {
      "int":
      r"key `foo` exists but does not result in an integer value reason: [B,b]ad any[\ ,_]cast",
      "uint overflow":
      "key `foo` exists but does not result in an unsigned integer value reason: bad numeric conversion: negative overflow",
      "int overflow":
      "key `foo` exists but does not result in an integer value reason: bad numeric conversion: positive overflow",
      "double":
      r"key `foo` exists but does not result in a double value reason: [B,b]ad any[\ ,_]cast",
      "bool":
      r"key `foo` exists but does not result in a True or False value reason: [B,b]ad any[\ ,_]cast"
    }

    for ob in [m, list(m.GetAtoms())[0], list(m.GetBonds())[0]]:
      ob.SetDoubleProp("foo", 2.0)
      with self.assertRaises(ValueError) as e:
        ob.GetBoolProp("foo")
      self.assertRegex(str(e.exception), errors["bool"])

      with self.assertRaises(ValueError) as e:
        ob.GetIntProp("foo")
      self.assertRegex(str(e.exception), errors["int"])

      ob.SetBoolProp("foo", True)
      with self.assertRaises(ValueError) as e:
        ob.GetDoubleProp("foo")
      self.assertRegex(str(e.exception), errors["double"])

      with self.assertRaises(ValueError) as e:
        ob.GetIntProp("foo")
      self.assertRegex(str(e.exception), errors["int"])

      ob.SetIntProp("foo", -1)
      with self.assertRaises(ValueError) as e:
        ob.GetUnsignedProp("foo")
      self.assertEqual(str(e.exception), errors["uint overflow"])

      ob.SetUnsignedProp("foo", 4294967295)
      self.assertEqual(ob.GetUnsignedProp("foo"), 4294967295)
      with self.assertRaises(ValueError) as e:
        ob.GetIntProp("foo")
      self.assertEqual(str(e.exception), errors["int overflow"])

  def testInvariantException(self):
    m = Chem.MolFromSmiles("C")
    try:
      m.GetAtomWithIdx(3)
    except RuntimeError as e:
      import platform
      details = str(e)
      if platform.system() == 'Windows':
        details = details.replace('\\', '/')
      self.assertTrue("Code/GraphMol/ROMol.cpp".lower() in details.lower())
      self.assertTrue("Failed Expression: 3 < 1" in details)
      self.assertTrue("RDKIT:" in details)
      self.assertTrue(__version__ in details)

  def testGetSDText(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'NCI_aids_few.sdf')
    #fileN = "../FileParsers/test_data/NCI_aids_few.sdf"
    sdSup = Chem.SDMolSupplier(fileN)
    for m in sdSup:
      sdt = Chem.SDWriter.GetText(m)
      ts = Chem.SDMolSupplier()
      ts.SetData(sdt)
      nm = next(ts)
      self.assertEqual(Chem.MolToSmiles(m, True), Chem.MolToSmiles(nm, True))
      for pn in m.GetPropNames():
        self.assertTrue(nm.HasProp(pn))
        self.assertEqual(m.GetProp(pn), nm.GetProp(pn))

  def testUnfoldedRDKFingerprint(self):
    from rdkit.Chem import AllChem

    m = Chem.MolFromSmiles('c1ccccc1N')
    fp = AllChem.UnfoldedRDKFingerprintCountBased(m)
    fpDict = fp.GetNonzeroElements()
    self.assertEqual(len(fpDict.items()), 19)
    self.assertTrue(374073638 in fpDict)
    self.assertEqual(fpDict[374073638], 6)
    self.assertTrue(464351883 in fpDict)
    self.assertEqual(fpDict[464351883], 2)
    self.assertTrue(1949583554 in fpDict)
    self.assertEqual(fpDict[1949583554], 6)
    self.assertTrue(4105342207 in fpDict)
    self.assertEqual(fpDict[4105342207], 1)
    self.assertTrue(794080973 in fpDict)
    self.assertEqual(fpDict[794080973], 1)
    self.assertTrue(3826517238 in fpDict)
    self.assertEqual(fpDict[3826517238], 2)

    m = Chem.MolFromSmiles('Cl')
    fp = AllChem.UnfoldedRDKFingerprintCountBased(m)
    fpDict = fp.GetNonzeroElements()
    self.assertEqual(len(fpDict.items()), 0)

    m = Chem.MolFromSmiles('CCCO')
    aBits = {}
    fp = AllChem.UnfoldedRDKFingerprintCountBased(m, bitInfo=aBits)
    fpDict = fp.GetNonzeroElements()
    self.assertEqual(len(fpDict.items()), 5)
    self.assertTrue(1524090560 in fpDict)
    self.assertEqual(fpDict[1524090560], 1)
    self.assertTrue(1940446997 in fpDict)
    self.assertEqual(fpDict[1940446997], 1)
    self.assertTrue(3977409745 in fpDict)
    self.assertEqual(fpDict[3977409745], 1)
    self.assertTrue(4274652475 in fpDict)
    self.assertEqual(fpDict[4274652475], 1)
    self.assertTrue(4275705116 in fpDict)
    self.assertEqual(fpDict[4275705116], 2)

    self.assertTrue(1524090560 in aBits)
    self.assertEqual(aBits[1524090560], [[1, 2]])
    self.assertTrue(1940446997 in aBits)
    self.assertEqual(aBits[1940446997], [[0, 1]])
    self.assertTrue(3977409745 in aBits)
    self.assertEqual(aBits[3977409745], [[0, 1, 2]])
    self.assertTrue(4274652475 in aBits)
    self.assertEqual(aBits[4274652475], [[2]])
    self.assertTrue(4275705116 in aBits)
    self.assertEqual(aBits[4275705116], [[0], [1]])

  def testRDKFingerprintBitInfo(self):

    m = Chem.MolFromSmiles('CCCO')
    aBits = {}
    fp1 = Chem.RDKFingerprint(m, bitInfo=aBits)
    self.assertTrue(1183 in aBits)
    self.assertEqual(aBits[1183], [[1, 2]])
    self.assertTrue(709 in aBits)
    self.assertEqual(aBits[709], [[0, 1]])
    self.assertTrue(1118 in aBits)
    self.assertEqual(aBits[1118], [[0, 1, 2]])
    self.assertTrue(562 in aBits)
    self.assertEqual(aBits[562], [[2]])
    self.assertTrue(1772 in aBits)
    self.assertEqual(aBits[1772], [[0], [1]])

  def testSimpleAromaticity(self):
    m = Chem.MolFromSmiles('c1ccccc1')
    self.assertTrue(m.GetBondWithIdx(0).GetIsAromatic())
    self.assertTrue(m.GetAtomWithIdx(0).GetIsAromatic())
    Chem.Kekulize(m, True)
    self.assertFalse(m.GetBondWithIdx(0).GetIsAromatic())
    self.assertFalse(m.GetAtomWithIdx(0).GetIsAromatic())
    Chem.SetAromaticity(m, Chem.AROMATICITY_SIMPLE)
    self.assertTrue(m.GetBondWithIdx(0).GetIsAromatic())
    self.assertTrue(m.GetAtomWithIdx(0).GetIsAromatic())

    m = Chem.MolFromSmiles('c1c[nH]cc1')
    self.assertTrue(m.GetBondWithIdx(0).GetIsAromatic())
    self.assertTrue(m.GetAtomWithIdx(0).GetIsAromatic())
    Chem.Kekulize(m, True)
    self.assertFalse(m.GetBondWithIdx(0).GetIsAromatic())
    self.assertFalse(m.GetAtomWithIdx(0).GetIsAromatic())
    Chem.SetAromaticity(m, Chem.AROMATICITY_SIMPLE)
    self.assertTrue(m.GetBondWithIdx(0).GetIsAromatic())
    self.assertTrue(m.GetAtomWithIdx(0).GetIsAromatic())

    m = Chem.MolFromSmiles('c1cccoocc1')
    self.assertTrue(m.GetBondWithIdx(0).GetIsAromatic())
    self.assertTrue(m.GetAtomWithIdx(0).GetIsAromatic())
    Chem.Kekulize(m, True)
    self.assertFalse(m.GetBondWithIdx(0).GetIsAromatic())
    self.assertFalse(m.GetAtomWithIdx(0).GetIsAromatic())
    Chem.SetAromaticity(m, Chem.AROMATICITY_SIMPLE)
    self.assertFalse(m.GetBondWithIdx(0).GetIsAromatic())
    self.assertFalse(m.GetAtomWithIdx(0).GetIsAromatic())

    m = Chem.MolFromSmiles('c1ooc1')
    self.assertTrue(m.GetBondWithIdx(0).GetIsAromatic())
    self.assertTrue(m.GetAtomWithIdx(0).GetIsAromatic())
    Chem.Kekulize(m, True)
    self.assertFalse(m.GetBondWithIdx(0).GetIsAromatic())
    self.assertFalse(m.GetAtomWithIdx(0).GetIsAromatic())
    Chem.SetAromaticity(m, Chem.AROMATICITY_SIMPLE)
    self.assertFalse(m.GetBondWithIdx(0).GetIsAromatic())
    self.assertFalse(m.GetAtomWithIdx(0).GetIsAromatic())

    m = Chem.MolFromSmiles('C1=CC2=CC=CC=CC2=C1')
    self.assertTrue(m.GetBondWithIdx(0).GetIsAromatic())
    self.assertTrue(m.GetAtomWithIdx(0).GetIsAromatic())
    Chem.Kekulize(m, True)
    self.assertFalse(m.GetBondWithIdx(0).GetIsAromatic())
    self.assertFalse(m.GetAtomWithIdx(0).GetIsAromatic())
    Chem.SetAromaticity(m, Chem.AROMATICITY_SIMPLE)
    self.assertFalse(m.GetBondWithIdx(0).GetIsAromatic())
    self.assertFalse(m.GetAtomWithIdx(0).GetIsAromatic())

  def testGithub955(self):
    m = Chem.MolFromSmiles("CCC")
    m.GetAtomWithIdx(0).SetProp("foo", "1")
    self.assertEqual(list(m.GetAtomWithIdx(0).GetPropNames()), ["foo"])
    m.GetBondWithIdx(0).SetProp("foo", "1")
    self.assertEqual(list(m.GetBondWithIdx(0).GetPropNames()), ["foo"])

  def testMDLProps(self):
    m = Chem.MolFromSmiles("CCC")
    m.GetAtomWithIdx(0).SetAtomMapNum(1)
    Chem.SetAtomAlias(m.GetAtomWithIdx(1), "foo")
    Chem.SetAtomValue(m.GetAtomWithIdx(1), "bar")

    m = Chem.MolFromMolBlock(Chem.MolToMolBlock(m))
    self.assertEqual(m.GetAtomWithIdx(0).GetAtomMapNum(), 1)
    self.assertEqual(Chem.GetAtomAlias(m.GetAtomWithIdx(1)), "foo")
    self.assertEqual(Chem.GetAtomValue(m.GetAtomWithIdx(1)), "bar")

  def testSmilesProps(self):
    m = Chem.MolFromSmiles("C")
    Chem.SetSupplementalSmilesLabel(m.GetAtomWithIdx(0), 'xxx')
    self.assertEqual(Chem.MolToSmiles(m), "Cxxx")

  def testGithub1051(self):
    # just need to test that this exists:
    self.assertTrue(Chem.BondDir.EITHERDOUBLE)

  def testGithub1041(self):
    a = Chem.Atom(6)
    self.assertRaises(RuntimeError, lambda: a.GetOwningMol())
    self.assertRaises(RuntimeError, lambda: a.GetNeighbors())
    self.assertRaises(RuntimeError, lambda: a.GetBonds())
    self.assertRaises(RuntimeError, lambda: a.IsInRing())
    self.assertRaises(RuntimeError, lambda: a.IsInRingSize(4))

  def testSmilesParseParams(self):
    smi = "CCC |$foo;;bar$| ourname"
    m = Chem.MolFromSmiles(smi)
    self.assertTrue(m is not None)
    ps = Chem.SmilesParserParams()
    ps.allowCXSMILES = False
    ps.parseName = False
    m = Chem.MolFromSmiles(smi, ps)
    self.assertTrue(m is None)
    ps.allowCXSMILES = True
    ps.parseName = True
    m = Chem.MolFromSmiles(smi, ps)
    self.assertTrue(m is not None)
    self.assertTrue(m.GetAtomWithIdx(0).HasProp('atomLabel'))
    self.assertEqual(m.GetAtomWithIdx(0).GetProp('atomLabel'), "foo")
    self.assertTrue(m.HasProp('_Name'))
    self.assertEqual(m.GetProp('_Name'), "ourname")
    self.assertEqual(m.GetProp("_CXSMILES_Data"), "|$foo;;bar$|")

  def testWriteCXSmiles(self):
    smi = "CCC |$foo;;bar$|"
    ps = Chem.SmilesParserParams()
    ps.allowCXSMILES = True
    m = Chem.MolFromSmiles(smi, ps)
    self.assertTrue(m is not None)
    self.assertTrue(m.GetAtomWithIdx(0).HasProp('atomLabel'))
    self.assertEqual(m.GetAtomWithIdx(0).GetProp('atomLabel'), "foo")
    self.assertEqual(Chem.MolToCXSmiles(m), 'CCC |$foo;;bar$|')

    smi = "Cl.CCC |$;foo;;bar$|"
    m = Chem.MolFromSmiles(smi, ps)
    self.assertTrue(m is not None)
    self.assertTrue(m.GetAtomWithIdx(1).HasProp('atomLabel'))
    self.assertEqual(m.GetAtomWithIdx(1).GetProp('atomLabel'), "foo")
    self.assertEqual(Chem.MolFragmentToCXSmiles(m, atomsToUse=(1, 2, 3)), 'CCC |$foo;;bar$|')

  def testPickleProps(self):
    import pickle
    m = Chem.MolFromSmiles('C1=CN=CC=C1')
    m.SetProp("_Name", "Name")
    for atom in m.GetAtoms():
      atom.SetProp("_foo", "bar" + str(atom.GetIdx()))
      atom.SetProp("foo", "baz" + str(atom.GetIdx()))

    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
    pkl = pickle.dumps(m)
    m2 = pickle.loads(pkl)
    smi1 = Chem.MolToSmiles(m)
    smi2 = Chem.MolToSmiles(m2)
    self.assertTrue(smi1 == smi2)
    self.assertEqual(m2.GetProp("_Name"), "Name")
    for atom in m2.GetAtoms():
      self.assertEqual(atom.GetProp("_foo"), "bar" + str(atom.GetIdx()))
      self.assertEqual(atom.GetProp("foo"), "baz" + str(atom.GetIdx()))

    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AtomProps)
    pkl = pickle.dumps(m)
    m2 = pickle.loads(pkl)
    smi1 = Chem.MolToSmiles(m)
    smi2 = Chem.MolToSmiles(m2)
    self.assertTrue(smi1 == smi2)
    self.assertFalse(m2.HasProp("_Name"))
    for atom in m2.GetAtoms():
      self.assertFalse(atom.HasProp("_foo"))
      self.assertEqual(atom.GetProp("foo"), "baz" + str(atom.GetIdx()))

    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.NoProps)
    pkl = pickle.dumps(m)
    m2 = pickle.loads(pkl)
    smi1 = Chem.MolToSmiles(m)
    smi2 = Chem.MolToSmiles(m2)
    self.assertTrue(smi1 == smi2)
    self.assertFalse(m2.HasProp("_Name"))
    for atom in m2.GetAtoms():
      self.assertFalse(atom.HasProp("_foo"))
      self.assertFalse(atom.HasProp("foo"))

    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.MolProps
                                    | Chem.PropertyPickleOptions.PrivateProps)
    pkl = pickle.dumps(m)
    m2 = pickle.loads(pkl)
    smi1 = Chem.MolToSmiles(m)
    smi2 = Chem.MolToSmiles(m2)
    self.assertTrue(smi1 == smi2)
    self.assertEqual(m2.GetProp("_Name"), "Name")
    for atom in m2.GetAtoms():
      self.assertFalse(atom.HasProp("_foo"))
      self.assertFalse(atom.HasProp("foo"))

  def testGithub1352(self):
    self.assertTrue('SP' in Chem.HybridizationType.names)
    self.assertTrue('S' in Chem.HybridizationType.names)
    m = Chem.MolFromSmiles('CC(=O)O.[Na]')
    self.assertEqual(m.GetAtomWithIdx(0).GetHybridization().name, 'SP3')
    self.assertEqual(m.GetAtomWithIdx(4).GetHybridization().name, 'S')

  def testGithub1366(self):
    mol = Chem.MolFromSmiles('*C*')
    mol = Chem.RWMol(mol)
    ats = iter(mol.GetAtoms())
    atom = next(ats)
    mol.RemoveAtom(atom.GetIdx())
    self.assertRaises(RuntimeError, next, ats)

    mol = Chem.MolFromSmiles('*C*')
    mol = Chem.RWMol(mol)
    bonds = iter(mol.GetBonds())
    bond = next(bonds)
    mol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    self.assertRaises(RuntimeError, next, bonds)

  def testGithub1478(self):
    data = """
  MJ150720

  8  8  0  0  0  0  0  0  0  0999 V2000
   -0.4242   -1.4883    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2901   -1.0758    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0046    0.9865    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
    1.0046    0.1614    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
    0.2901   -0.2508    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4243    0.1614    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4243    0.9865    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
    0.2901    1.3990    0.0000 A   0  0  0  0  0  0  0  0  0  0  0  0
  7  6  4  0  0  0  0
  8  7  4  0  0  0  0
  6  5  4  0  0  0  0
  5  4  4  0  0  0  0
  5  2  1  0  0  0  0
  4  3  4  0  0  0  0
  8  3  4  0  0  0  0
  2  1  2  0  0  0  0
M  END
"""
    pattern = Chem.MolFromMolBlock(data)
    m = Chem.MolFromSmiles("c1ccccc1C=O")
    self.assertTrue(m.HasSubstructMatch(pattern))

  def testGithub1320(self):
    import pickle
    mol = Chem.MolFromSmiles('N[C@@H](C)O')
    mol2 = pickle.loads(pickle.dumps(mol))
    self.assertEqual(Chem.MolToSmiles(mol, isomericSmiles=True),
                     Chem.MolToSmiles(mol2, isomericSmiles=True))
    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AtomProps
                                    | Chem.PropertyPickleOptions.BondProps
                                    | Chem.PropertyPickleOptions.MolProps
                                    | Chem.PropertyPickleOptions.PrivateProps
                                    | Chem.PropertyPickleOptions.ComputedProps)
    mol3 = pickle.loads(pickle.dumps(mol))

    for a1, a2 in zip(mol.GetAtoms(), mol3.GetAtoms()):
      d1 = a1.GetPropsAsDict()
      d2 = a2.GetPropsAsDict()
      if "__computedProps" in d1:
        c1 = list(d1["__computedProps"])
        c2 = list(d2["__computedProps"])
        del d1["__computedProps"]
        del d2["__computedProps"]
        self.assertEqual(c1, c2)

      assert d1 == d2

    for a1, a2 in zip(mol.GetBonds(), mol3.GetBonds()):
      d1 = a1.GetPropsAsDict()
      d2 = a2.GetPropsAsDict()
      if "__computedProps" in d1:
        c1 = list(d1["__computedProps"])
        c2 = list(d2["__computedProps"])
        del d1["__computedProps"]
        del d2["__computedProps"]
        self.assertEqual(c1, c2)

      assert d1 == d2

    self.assertEqual(Chem.MolToSmiles(mol, isomericSmiles=True),
                     Chem.MolToSmiles(mol3, isomericSmiles=True))

  def testOldPropPickles(self):
    data = 'crdkit.Chem.rdchem\nMol\np0\n(S\'\\xef\\xbe\\xad\\xde\\x00\\x00\\x00\\x00\\x08\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00\\x00)\\x00\\x00\\x00-\\x00\\x00\\x00\\x80\\x01\\x06\\x00`\\x00\\x00\\x00\\x01\\x03\\x07\\x00`\\x00\\x00\\x00\\x02\\x01\\x06 4\\x00\\x00\\x00\\x01\\x01\\x04\\x06\\x00`\\x00\\x00\\x00\\x01\\x03\\x06\\x00(\\x00\\x00\\x00\\x03\\x04\\x08\\x00(\\x00\\x00\\x00\\x03\\x02\\x07\\x00h\\x00\\x00\\x00\\x03\\x02\\x01\\x06 4\\x00\\x00\\x00\\x02\\x01\\x04\\x06\\x00(\\x00\\x00\\x00\\x03\\x04\\x08\\x00(\\x00\\x00\\x00\\x03\\x02\\x07\\x00(\\x00\\x00\\x00\\x03\\x03\\x06\\x00`\\x00\\x00\\x00\\x02\\x02\\x06 4\\x00\\x00\\x00\\x01\\x01\\x04\\x08\\x00(\\x00\\x00\\x00\\x03\\x02\\x06@(\\x00\\x00\\x00\\x03\\x04\\x06@h\\x00\\x00\\x00\\x03\\x03\\x01\\x06@h\\x00\\x00\\x00\\x03\\x03\\x01\\x06@h\\x00\\x00\\x00\\x03\\x03\\x01\\x06@h\\x00\\x00\\x00\\x03\\x03\\x01\\x06@h\\x00\\x00\\x00\\x03\\x03\\x01\\x06\\x00`\\x00\\x00\\x00\\x02\\x02\\x06 4\\x00\\x00\\x00\\x01\\x01\\x04\\x06\\x00(\\x00\\x00\\x00\\x03\\x04\\x08\\x00(\\x00\\x00\\x00\\x03\\x02\\x07\\x00h\\x00\\x00\\x00\\x03\\x02\\x01\\x06 4\\x00\\x00\\x00\\x02\\x01\\x04\\x06\\x00`\\x00\\x00\\x00\\x02\\x02\\x06\\x00`\\x00\\x00\\x00\\x02\\x02\\x06\\x00`\\x00\\x00\\x00\\x02\\x02\\x06@(\\x00\\x00\\x00\\x03\\x04\\x06@h\\x00\\x00\\x00\\x03\\x03\\x01\\x06@h\\x00\\x00\\x00\\x03\\x03\\x01\\x06@h\\x00\\x00\\x00\\x03\\x03\\x01\\x06@h\\x00\\x00\\x00\\x03\\x03\\x01\\x06@(\\x00\\x00\\x00\\x03\\x04\\x06\\x00`\\x00\\x00\\x00\\x03\\x01\\x06\\x00`\\x00\\x00\\x00\\x02\\x02\\x06\\x00`\\x00\\x00\\x00\\x02\\x02\\x06\\x00`\\x00\\x00\\x00\\x02\\x02\\x06\\x00`\\x00\\x00\\x00\\x02\\x02\\x06\\x00`\\x00\\x00\\x00\\x02\\x02\\x0b\\x00\\x01\\x00\\x01\\x02\\x00\\x02\\x03\\x00\\x02\\x04\\x00\\x04\\x05(\\x02\\x04\\x06 \\x06\\x07\\x00\\x07\\x08\\x00\\x08\\t(\\x02\\x08\\n \\n\\x0b\\x00\\x0b\\x0c\\x00\\x0c\\r\\x00\\r\\x0e \\x0e\\x0fh\\x0c\\x0f\\x10h\\x0c\\x10\\x11h\\x0c\\x11\\x12h\\x0c\\x12\\x13h\\x0c\\x0c\\x14\\x00\\x14\\x15\\x00\\x15\\x16\\x00\\x16\\x17(\\x02\\x16\\x18 \\x18\\x19\\x00\\x19\\x1a\\x00\\x1a\\x1b\\x00\\x1b\\x1c\\x00\\x1c\\x1d\\x00\\x1d\\x1eh\\x0c\\x1e\\x1fh\\x0c\\x1f h\\x0c !h\\x0c!"h\\x0c\\x07#\\x00#$\\x00$%\\x00%&\\x00&\\\'\\x00\\\'(\\x00\\x15\\n\\x00"\\x19\\x00(#\\x00\\x13\\x0eh\\x0c"\\x1dh\\x0c\\x14\\x05\\x05\\x0b\\n\\x15\\x14\\x0c\\x06\\x0f\\x10\\x11\\x12\\x13\\x0e\\x06\\x1a\\x1b\\x1c\\x1d"\\x19\\x06\\x1e\\x1f !"\\x1d\\x06$%&\\\'(#\\x17\\x00\\x00\\x00\\x00\\x12\\x03\\x00\\x00\\x00\\x07\\x00\\x00\\x00numArom\\x01\\x02\\x00\\x00\\x00\\x0f\\x00\\x00\\x00_StereochemDone\\x01\\x01\\x00\\x00\\x00\\x03\\x00\\x00\\x00foo\\x00\\x03\\x00\\x00\\x00bar\\x13:\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x12\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x01\\x00\\x00\\x000\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x1d\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x01\\x00\\x00\\x001\\x04\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x15\\x00\\x00\\x00\\x12\\x00\\x00\\x00_ChiralityPossible\\x01\\x01\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPCode\\x00\\x01\\x00\\x00\\x00S\\x05\\x00\\x00\\x00myidx\\x00\\x01\\x00\\x00\\x002\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x00\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x01\\x00\\x00\\x003\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x1a\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x01\\x00\\x00\\x004\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02"\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x01\\x00\\x00\\x005\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x1f\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x01\\x00\\x00\\x006\\x04\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x16\\x00\\x00\\x00\\x12\\x00\\x00\\x00_ChiralityPossible\\x01\\x01\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPCode\\x00\\x01\\x00\\x00\\x00S\\x05\\x00\\x00\\x00myidx\\x00\\x01\\x00\\x00\\x007\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x1c\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x01\\x00\\x00\\x008\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02$\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x01\\x00\\x00\\x009\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02 \\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0010\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x13\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0011\\x04\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x18\\x00\\x00\\x00\\x12\\x00\\x00\\x00_ChiralityPossible\\x01\\x01\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPCode\\x00\\x01\\x00\\x00\\x00S\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0012\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02!\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0013\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x19\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0014\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x0f\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0015\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x0b\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0016\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x08\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0017\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x0b\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0018\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x0f\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0019\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x07\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0020\\x04\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x17\\x00\\x00\\x00\\x12\\x00\\x00\\x00_ChiralityPossible\\x01\\x01\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPCode\\x00\\x01\\x00\\x00\\x00S\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0021\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x1b\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0022\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02#\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0023\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x1e\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0024\\x04\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x14\\x00\\x00\\x00\\x12\\x00\\x00\\x00_ChiralityPossible\\x01\\x01\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPCode\\x00\\x01\\x00\\x00\\x00R\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0025\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x06\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0026\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x03\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0027\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x05\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0028\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x10\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0029\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x0c\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0030\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\t\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0031\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\n\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0032\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\r\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0033\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x11\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0034\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x0e\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0035\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x04\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0036\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x02\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0037\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x01\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0038\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x02\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0039\\x02\\x00\\x00\\x00\\x08\\x00\\x00\\x00_CIPRank\\x02\\x04\\x00\\x00\\x00\\x05\\x00\\x00\\x00myidx\\x00\\x02\\x00\\x00\\x0040\\x13\\x16\'\np1\ntp2\nRp3\n.'
    import pickle

    # bonds were broken in v1
    m2 = pickle.loads(data.encode("utf-8"), encoding='bytes')

    self.assertEqual(m2.GetProp("foo"), "bar")
    for atom in m2.GetAtoms():
      self.assertEqual(atom.GetProp("myidx"), str(atom.GetIdx()))

    self.assertEqual(
      Chem.MolToSmiles(m2, True),
      Chem.MolToSmiles(
        Chem.MolFromSmiles(
          "CN[C@@H](C)C(=O)N[C@H](C(=O)N1C[C@@H](Oc2ccccc2)C[C@H]1C(=O)N[C@@H]1CCCc2ccccc21)C1CCCCC1"
        ), True))

  def testGithub1461(self):
    # this is simple, it should throw a precondition and not seg fault
    m = Chem.RWMol()
    try:
      m.AddBond(0, 1, Chem.BondType.SINGLE)
      self.assertFalse(True)  # shouldn't get here
    except RuntimeError:
      pass

  def testMolBundles1(self):
    b = Chem.MolBundle()
    smis = ('CC(Cl)(F)CC(F)(Br)', 'C[C@](Cl)(F)C[C@H](F)(Br)', 'C[C@](Cl)(F)C[C@@H](F)(Br)')
    for smi in smis:
      b.AddMol(Chem.MolFromSmiles(smi))
    self.assertEqual(len(b), 3)
    self.assertEqual(b.Size(), 3)
    self.assertRaises(IndexError, lambda: b[4])
    self.assertEqual(Chem.MolToSmiles(b[1], isomericSmiles=True),
                     Chem.MolToSmiles(Chem.MolFromSmiles(smis[1]), isomericSmiles=True))
    self.assertTrue(b.HasSubstructMatch(Chem.MolFromSmiles('CC(Cl)(F)CC(F)(Br)'),
                                        useChirality=True))
    self.assertTrue(
      b.HasSubstructMatch(Chem.MolFromSmiles('C[C@](Cl)(F)C[C@@H](F)(Br)'), useChirality=True))
    self.assertTrue(
      b.HasSubstructMatch(Chem.MolFromSmiles('C[C@@](Cl)(F)C[C@@H](F)(Br)'), useChirality=False))
    self.assertFalse(
      b.HasSubstructMatch(Chem.MolFromSmiles('C[C@@](Cl)(F)C[C@@H](F)(Br)'), useChirality=True))

    self.assertEqual(
      len(b.GetSubstructMatch(Chem.MolFromSmiles('CC(Cl)(F)CC(F)(Br)'), useChirality=True)), 8)
    self.assertEqual(
      len(b.GetSubstructMatch(Chem.MolFromSmiles('C[C@](Cl)(F)C[C@@H](F)(Br)'), useChirality=True)),
      8)
    self.assertEqual(
      len(b.GetSubstructMatch(Chem.MolFromSmiles('C[C@@](Cl)(F)C[C@@H](F)(Br)'),
                              useChirality=False)), 8)
    self.assertEqual(
      len(b.GetSubstructMatch(Chem.MolFromSmiles('C[C@@](Cl)(F)C[C@@H](F)(Br)'),
                              useChirality=True)), 0)

    self.assertEqual(
      len(b.GetSubstructMatches(Chem.MolFromSmiles('CC(Cl)(F)CC(F)(Br)'), useChirality=True)), 1)
    self.assertEqual(
      len(b.GetSubstructMatches(Chem.MolFromSmiles('C[C@](Cl)(F)C[C@@H](F)(Br)'),
                                useChirality=True)), 1)
    self.assertEqual(
      len(
        b.GetSubstructMatches(Chem.MolFromSmiles('C[C@@](Cl)(F)C[C@@H](F)(Br)'),
                              useChirality=False)), 1)
    self.assertEqual(
      len(
        b.GetSubstructMatches(Chem.MolFromSmiles('C[C@@](Cl)(F)C[C@@H](F)(Br)'),
                              useChirality=True)), 0)
    self.assertEqual(
      len(b.GetSubstructMatches(Chem.MolFromSmiles('CC(Cl)(F)CC(F)(Br)'), useChirality=True)[0]), 8)
    self.assertEqual(
      len(
        b.GetSubstructMatches(Chem.MolFromSmiles('C[C@](Cl)(F)C[C@@H](F)(Br)'),
                              useChirality=True)[0]), 8)
    self.assertEqual(
      len(
        b.GetSubstructMatches(Chem.MolFromSmiles('C[C@@](Cl)(F)C[C@@H](F)(Br)'),
                              useChirality=False)[0]), 8)

    if Chem.MolBundleCanSerialize():
      for b2 in (pickle.loads(pickle.dumps(b)), Chem.MolBundle(b.ToBinary())):
        self.assertEqual(len(b2), len(b))
        for m, m2 in zip(b, b2):
          self.assertEqual(Chem.MolToSmiles(m), Chem.MolToSmiles(m2))

  def testMolBundles2(self):
    b = Chem.MolBundle()
    smis = ('Fc1c(Cl)cccc1', 'Fc1cc(Cl)ccc1', 'Fc1ccc(Cl)cc1')
    for smi in smis:
      b.AddMol(Chem.MolFromSmiles(smi))
    self.assertEqual(len(b), 3)
    self.assertEqual(b.Size(), 3)
    self.assertTrue(Chem.MolFromSmiles('Fc1c(Cl)cccc1').HasSubstructMatch(b))
    self.assertTrue(Chem.MolFromSmiles('Fc1cc(Cl)ccc1').HasSubstructMatch(b))
    self.assertTrue(Chem.MolFromSmiles('Fc1c(Cl)cccc1C').HasSubstructMatch(b))
    self.assertTrue(Chem.MolFromSmiles('Fc1cc(Cl)ccc1C').HasSubstructMatch(b))
    self.assertFalse(Chem.MolFromSmiles('Fc1c(Br)cccc1').HasSubstructMatch(b))

    self.assertEqual(len(Chem.MolFromSmiles('Fc1c(Cl)cccc1').GetSubstructMatch(b)), 8)
    self.assertEqual(len(Chem.MolFromSmiles('Fc1c(Cl)cccc1').GetSubstructMatches(b)), 1)
    self.assertEqual(len(Chem.MolFromSmiles('Fc1c(Cl)cccc1').GetSubstructMatches(b)[0]), 8)
    self.assertEqual(len(Chem.MolFromSmiles('Fc1ccc(Cl)cc1').GetSubstructMatches(b)), 1)
    self.assertEqual(
      len(Chem.MolFromSmiles('Fc1ccc(Cl)cc1').GetSubstructMatches(b, uniquify=False)), 2)

    self.assertEqual(len(Chem.MolFromSmiles('Fc1c(C)cccc1').GetSubstructMatch(b)), 0)
    self.assertEqual(len(Chem.MolFromSmiles('Fc1c(C)cccc1').GetSubstructMatches(b)), 0)

  def testMolBundles3(self):
    smis = ('CCC', 'CCO', 'CCN')
    b = Chem.FixedMolSizeMolBundle()
    for smi in smis:
      b.AddMol(Chem.MolFromSmiles(smi))
    self.assertEqual(len(b), 3)
    self.assertEqual(b.Size(), 3)
    with self.assertRaises(ValueError):
      b.AddMol(Chem.MolFromSmiles('CCCC'))

    b = Chem.MolBundle()
    for smi in smis:
      b.AddMol(Chem.MolFromSmiles(smi))
    self.assertEqual(len(b), 3)
    self.assertEqual(b.Size(), 3)
    b.AddMol(Chem.MolFromSmiles('CCCC'))
    self.assertEqual(b.Size(), 4)

  def testGithub1622(self):
    nonaromatics = (
      "C1=C[N]C=C1",  # radicals are not two electron donors
      "O=C1C=CNC=C1",  # exocyclic double bonds don't steal electrons
      "C1=CS(=O)C=C1",  # not sure how to classify this example from the
      # OEChem docs
      "C1#CC=CC=C1"  # benzyne
      # 5-membered heterocycles
      "C1=COC=C1",  # furan
      "C1=CSC=C1",  # thiophene
      "C1=CNC=C1",  #pyrrole
      "C1=COC=N1",  # oxazole
      "C1=CSC=N1",  # thiazole
      "C1=CNC=N1",  # imidazole
      "C1=CNN=C1",  # pyrazole
      "C1=CON=C1",  # isoxazole
      "C1=CSN=C1",  # isothiazole
      "C1=CON=N1",  # 1,2,3-oxadiazole
      "C1=CNN=N1",  # 1,2,3-triazole
      "N1=CSC=N1",  # 1,3,4-thiadiazole
      # not outside the second rows
      "C1=CC=C[Si]=C1",
      "C1=CC=CC=P1",
      # 5-membered heterocycles outside the second row
      "C1=C[Se]C=C1",
      'C1=C[Te]C=C1')
    for smi in nonaromatics:
      m = Chem.MolFromSmiles(smi, sanitize=False)
      Chem.SanitizeMol(m, Chem.SANITIZE_ALL ^ Chem.SANITIZE_SETAROMATICITY)
      Chem.SetAromaticity(m, Chem.AROMATICITY_MDL)
      self.assertFalse(m.GetAtomWithIdx(0).GetIsAromatic())
    aromatics = (
      "C1=CC=CC=C1",  # benzene, of course
      # hetrocyclics
      "N1=CC=CC=C1",  # pyridine
      "N1=CC=CC=N1",  # pyridazine
      "N1=CC=CN=C1",  # pyrimidine
      "N1=CC=NC=C1",  # pyrazine
      "N1=CN=CN=C1",  # 1,3,5-triazine
      # polycyclic aromatics
      "C1=CC2=CC=CC=CC2=C1",  # azulene
      "C1=CC=CC2=CC=CC=C12",
      "C1=CC2=CC=CC=CC=C12",
      "C1=CC=C2C(=C1)N=CC=N2",
      "C1=CN=CC2C=CC=CC1=2",
      "C1=CC=C2C(=C1)N=C3C=CC=CC3=N2",
      "C1=CN=NC2C=CC=CC1=2",
      # macrocycle aromatics
      "C1=CC=CC=CC=CC=C1",
      "C1=CC=CC=CC=CC=CC=CC=CC=CC=C1",
      "N1=CN=NC=CC=CC=CC=CC=CC=CC=CC=CC=CC=CC=CC=CC=C1")
    for smi in aromatics:
      m = Chem.MolFromSmiles(smi, sanitize=False)
      Chem.SanitizeMol(m, Chem.SANITIZE_ALL ^ Chem.SANITIZE_SETAROMATICITY)
      Chem.SetAromaticity(m, Chem.AROMATICITY_MDL)
      self.assertTrue(m.GetAtomWithIdx(0).GetIsAromatic())

  def testMolBlockChirality(self):
    m = Chem.MolFromSmiles('C[C@H](Cl)Br')
    mb = Chem.MolToMolBlock(m)
    m2 = Chem.MolFromMolBlock(mb)
    csmi1 = Chem.MolToSmiles(m, isomericSmiles=True)
    csmi2 = Chem.MolToSmiles(m2, isomericSmiles=True)
    self.assertEqual(csmi1, csmi2)

  def testIssue1735(self):
    # this shouldn't seg fault...
    m = Chem.RWMol()
    ranks = Chem.CanonicalRankAtoms(m, breakTies=False)
    ranks = Chem.CanonicalRankAtoms(m, breakTies=True)

  def testGithub1615(self):
    mb = """Issue399a.mol
  ChemDraw04050615582D

  4  4  0  0  0  0  0  0  0  0999 V2000
   -0.7697    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0553    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7697    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7697   -0.4125    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  0
  2  3  1  0
  3  4  1  0
  2  4  1  0
M  END"""
    m = Chem.MolFromMolBlock(mb)
    self.assertFalse(m.GetAtomWithIdx(1).HasProp("_CIPCode"))
    self.assertEqual(m.GetBondWithIdx(0).GetBondDir(), Chem.BondDir.NONE)
    self.assertEqual(m.GetAtomWithIdx(1).GetChiralTag(), Chem.ChiralType.CHI_UNSPECIFIED)
    m.GetAtomWithIdx(1).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    Chem.AssignStereochemistry(m, force=True)
    self.assertTrue(m.GetAtomWithIdx(1).HasProp("_CIPCode"))
    self.assertEqual(m.GetAtomWithIdx(1).GetProp("_CIPCode"), "S")
    self.assertEqual(m.GetBondWithIdx(0).GetBondDir(), Chem.BondDir.NONE)
    Chem.WedgeBond(m.GetBondWithIdx(0), 1, m.GetConformer())
    self.assertEqual(m.GetBondWithIdx(0).GetBondDir(), Chem.BondDir.BEGINWEDGE)

  def testSmilesToAtom(self):
    a = Chem.AtomFromSmiles("C")
    self.assertEqual(a.GetAtomicNum(), 6)
    b = Chem.BondFromSmiles("=")
    self.assertEqual(b.GetBondType(), Chem.BondType.DOUBLE)
    a = Chem.AtomFromSmiles("error")
    self.assertIs(a, None)
    b = Chem.BondFromSmiles("d")
    self.assertIs(b, None)

    a = Chem.AtomFromSmarts("C")
    self.assertEqual(a.GetAtomicNum(), 6)
    b = Chem.BondFromSmarts("=")
    self.assertEqual(b.GetBondType(), Chem.BondType.DOUBLE)
    a = Chem.AtomFromSmarts("error")
    self.assertIs(a, None)
    b = Chem.BondFromSmarts("d")
    self.assertIs(b, None)

  def testSVGParsing(self):
    svg = """<?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='200px' height='200px' >
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='200' height='200' x='0' y='0'> </rect>
<path d='M 9.09091,89.4974 24.2916,84.7462' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 24.2916,84.7462 39.4923,79.9949' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 86.2908,106.814 75.1709,93.4683 72.0765,96.8285 86.2908,106.814' style='fill:#000000;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 75.1709,93.4683 57.8622,86.8431 64.051,80.1229 75.1709,93.4683' style='fill:#0000FF;fill-rule:evenodd;stroke:#0000FF;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 75.1709,93.4683 72.0765,96.8285 57.8622,86.8431 75.1709,93.4683' style='fill:#0000FF;fill-rule:evenodd;stroke:#0000FF;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 86.2908,106.814 82.1459,125.293' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 82.1459,125.293 78.0009,143.772' style='fill:none;fill-rule:evenodd;stroke:#00CC00;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 86.2908,106.814 129.89,93.1862' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 134.347,94.186 138.492,75.7069' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 138.492,75.7069 142.637,57.2277' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 125.432,92.1865 129.577,73.7074' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 129.577,73.7074 133.722,55.2282' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 129.89,93.1862 142.557,104.852' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 142.557,104.852 155.224,116.517' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<text x='39.4923' y='83.483' style='font-size:15px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#0000FF' ><tspan>NH</tspan></text>
<text x='67.6656' y='158.998' style='font-size:15px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#00CC00' ><tspan>Cl</tspan></text>
<text x='132.777' y='56.228' style='font-size:15px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#FF0000' ><tspan>O</tspan></text>
<text x='149.782' y='131.743' style='font-size:15px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#FF0000' ><tspan>OH</tspan></text>
<text x='89.9952' y='194' style='font-size:12px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' ><tspan>m1</tspan></text>
<metadata>
<rdkit:mol xmlns:rdkit = "http://www.rdkit.org/xml" version="0.9">
<rdkit:atom idx="1" atom-smiles="[CH3]" drawing-x="9.09091" drawing-y="89.4974" x="-2.78651" y="0.295614" z="0" />
<rdkit:atom idx="2" atom-smiles="[NH]" drawing-x="52.6897" drawing-y="75.8699" x="-1.35482" y="0.743114" z="0" />
<rdkit:atom idx="3" atom-smiles="[C@H]" drawing-x="86.2908" drawing-y="106.814" x="-0.251428" y="-0.273019" z="0" />
<rdkit:atom idx="4" atom-smiles="[Cl]" drawing-x="76.2932" drawing-y="151.385" x="-0.579728" y="-1.73665" z="0" />
<rdkit:atom idx="5" atom-smiles="[C]" drawing-x="129.89" drawing-y="93.1862" x="1.18027" y="0.174481" z="0" />
<rdkit:atom idx="6" atom-smiles="[O]" drawing-x="139.887" drawing-y="48.6148" x="1.50857" y="1.63811" z="0" />
<rdkit:atom idx="7" atom-smiles="[OH]" drawing-x="163.491" drawing-y="124.13" x="2.28366" y="-0.841652" z="0" />
<rdkit:bond idx="1" begin-atom-idx="1" end-atom-idx="2" bond-smiles="-" />
<rdkit:bond idx="2" begin-atom-idx="2" end-atom-idx="3" bond-smiles="-" />
<rdkit:bond idx="3" begin-atom-idx="3" end-atom-idx="4" bond-smiles="-" />
<rdkit:bond idx="4" begin-atom-idx="3" end-atom-idx="5" bond-smiles="-" />
<rdkit:bond idx="5" begin-atom-idx="5" end-atom-idx="6" bond-smiles="=" />
<rdkit:bond idx="6" begin-atom-idx="5" end-atom-idx="7" bond-smiles="-" />
</rdkit:mol></metadata>
</svg>"""
    mol = Chem.MolFromRDKitSVG(svg)
    self.assertEqual(mol.GetNumAtoms(), 7)
    self.assertEqual(Chem.MolToSmiles(mol), 'CN[C@H](Cl)C(=O)O')

    svg2 = """<?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
          xmlns='http://www.w3.org/2000/svg'
                  xmlns:rdkit='http://www.rdkit.org/xml'
                  xmlns:xlink='http://www.w3.org/1999/xlink'
              xml:space='preserve'
width='200px' height='200px' >
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='200' height='200' x='0' y='0'> </rect>
<path d='M 9.09091,89.4974 24.2916,84.7462' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 24.2916,84.7462 39.4923,79.9949' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 86.2908,106.814 75.1709,93.4683 72.0765,96.8285 86.2908,106.814' style='fill:#000000;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 75.1709,93.4683 57.8622,86.8431 64.051,80.1229 75.1709,93.4683' style='fill:#0000FF;fill-rule:evenodd;stroke:#0000FF;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 75.1709,93.4683 72.0765,96.8285 57.8622,86.8431 75.1709,93.4683' style='fill:#0000FF;fill-rule:evenodd;stroke:#0000FF;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 86.2908,106.814 82.1459,125.293' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 82.1459,125.293 78.0009,143.772' style='fill:none;fill-rule:evenodd;stroke:#00CC00;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 86.2908,106.814 129.89,93.1862' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 134.347,94.186 138.492,75.7069' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 138.492,75.7069 142.637,57.2277' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 125.432,92.1865 129.577,73.7074' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 129.577,73.7074 133.722,55.2282' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 129.89,93.1862 142.557,104.852' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 142.557,104.852 155.224,116.517' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<text x='39.4923' y='83.483' style='font-size:15px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#0000FF' ><tspan>NH</tspan></text>
<text x='67.6656' y='158.998' style='font-size:15px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#00CC00' ><tspan>Cl</tspan></text>
<text x='132.777' y='56.228' style='font-size:15px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#FF0000' ><tspan>O</tspan></text>
<text x='149.782' y='131.743' style='font-size:15px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#FF0000' ><tspan>OH</tspan></text>
<text x='89.9952' y='194' style='font-size:12px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' ><tspan>m1</tspan></text>
</svg>"""
    mol = Chem.MolFromRDKitSVG(svg2)
    self.assertTrue(mol is None)

    with self.assertRaises(RuntimeError):
      mol = Chem.MolFromRDKitSVG("bad svg")

  def testAssignChiralTypesFromBondDirs(self):
    """
    Just check to see that AssignChiralTypesFromBondDirs is wrapped.
    Critical tests of the underlying C++ function already exist
    in SD file reader tests.
    """
    mol = Chem.MolFromSmiles('C(F)(Cl)Br')
    rdkit.Chem.rdDepictor.Compute2DCoords(mol)
    atom0 = mol.GetAtomWithIdx(0)
    self.assertEqual(atom0.GetChiralTag(), Chem.rdchem.ChiralType.CHI_UNSPECIFIED)
    bond = mol.GetBondBetweenAtoms(0, 1)
    bond.SetBondDir(Chem.rdchem.BondDir.BEGINWEDGE)
    Chem.AssignChiralTypesFromBondDirs(mol)
    self.assertEqual(atom0.GetChiralTag(), Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)

  def testAssignStereochemistryFrom3D(self):

    def _stereoTester(mol, expectedCIP, expectedStereo):
      mol.UpdatePropertyCache()
      self.assertEqual(mol.GetNumAtoms(), 9)
      self.assertFalse(mol.GetAtomWithIdx(1).HasProp("_CIPCode"))
      self.assertEqual(mol.GetBondWithIdx(3).GetStereo(), Chem.BondStereo.STEREONONE)
      for bond in mol.GetBonds():
        bond.SetBondDir(Chem.BondDir.NONE)

      Chem.AssignStereochemistryFrom3D(mol)
      self.assertTrue(mol.GetAtomWithIdx(1).HasProp("_CIPCode"))
      self.assertEqual(mol.GetAtomWithIdx(1).GetProp("_CIPCode"), expectedCIP)
      self.assertEqual(mol.GetBondWithIdx(3).GetStereo(), expectedStereo)

    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'test_data', 'stereochem.sdf')
    suppl = Chem.SDMolSupplier(fileN, sanitize=False)
    expected = (
      ("R", Chem.BondStereo.STEREOZ),
      ("R", Chem.BondStereo.STEREOE),
      ("S", Chem.BondStereo.STEREOZ),
      ("S", Chem.BondStereo.STEREOE),
    )
    for i, mol in enumerate(suppl):
      cip, stereo = expected[i]
      _stereoTester(mol, cip, stereo)

  def testGitHub2082(self):
    ctab = """
  MJ150720

  9  9  0  0  0  0  0  0  0  0999 V2000
    2.5687   -0.7144    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1562    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5687    0.7144    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.3312    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9187   -0.7144    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0937   -0.7144    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3187    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0937    0.7144    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.9187    0.7144    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  6
  2  3  1  0
  2  4  1  0
  4  5  2  0
  5  6  1  0
  6  7  2  0
  7  8  1  0
  8  9  2  0
  9  4  1  0
M  END
"""
    mol = Chem.MolFromMolBlock(ctab)
    self.assertFalse(mol.GetConformer().Is3D())
    self.assertTrue("@" in Chem.MolToSmiles(mol, True))

  def testGitHub2082_2(self):
    # test a mol block that lies is 3D but labelled 2D
    ofile = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Wrap', 'test_data',
                         'issue2082.mol')
    with open(ofile) as inf:
      ctab = inf.read()
    m = Chem.MolFromMolBlock(ctab)
    self.assertTrue(m.GetConformer().Is3D())

  def testSetQuery(self):
    from rdkit.Chem import rdqueries
    pat = Chem.MolFromSmarts("[C]")
    self.assertFalse(Chem.MolFromSmiles("c1ccccc1").HasSubstructMatch(pat))

    q = rdqueries.AtomNumEqualsQueryAtom(6)
    for atom in pat.GetAtoms():
      atom.SetQuery(q)

    self.assertTrue(Chem.MolFromSmiles("c1ccccc1").HasSubstructMatch(pat))

  def testGetQueryType(self):
    query_a = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                           'query_A.mol')
    m = next(Chem.SDMolSupplier(query_a))
    self.assertTrue(m.GetAtomWithIdx(6).HasQuery())
    self.assertTrue(m.GetAtomWithIdx(6).GetQueryType() == "A")

    query_a_v3k = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                               'query_A.v3k.mol')
    m = next(Chem.SDMolSupplier(query_a_v3k))
    self.assertTrue(m.GetAtomWithIdx(6).HasQuery())
    self.assertTrue(m.GetAtomWithIdx(6).GetQueryType() == "A")

    query_q = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                           'query_Q.mol')
    m = next(Chem.SDMolSupplier(query_q))
    self.assertTrue(m.GetAtomWithIdx(6).HasQuery())
    self.assertTrue(m.GetAtomWithIdx(6).GetQueryType() == "Q")

    query_q_v3k = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                               'query_Q.v3k.mol')
    m = next(Chem.SDMolSupplier(query_q_v3k))
    self.assertTrue(m.GetAtomWithIdx(6).HasQuery())
    self.assertTrue(m.GetAtomWithIdx(6).GetQueryType() == "Q")

    m = Chem.MolFromSmiles("*CC")
    params = Chem.rdmolops.AdjustQueryParameters.NoAdjustments()
    params.makeDummiesQueries = True
    m = Chem.rdmolops.AdjustQueryProperties(m, params)
    self.assertTrue(m.GetAtomWithIdx(0).HasQuery())
    self.assertTrue(m.GetAtomWithIdx(0).GetQueryType() == "")

  def testBondSetQuery(self):
    pat = Chem.MolFromSmarts('[#6]=[#6]')
    mol = Chem.MolFromSmiles("c1ccccc1")
    self.assertFalse(mol.HasSubstructMatch(pat))
    pat2 = Chem.MolFromSmarts('C:C')
    for bond in pat.GetBonds():
      bond.SetQuery(pat2.GetBondWithIdx(0))
    self.assertTrue(mol.HasSubstructMatch(pat))

  def testBondExpandQuery(self):
    pat = Chem.MolFromSmarts('C-C')
    mol = Chem.MolFromSmiles("C=C-C")
    self.assertEqual(len(mol.GetSubstructMatches(pat)), 1)
    pat2 = Chem.MolFromSmarts('C=C')
    for bond in pat.GetBonds():
      bond.ExpandQuery(pat2.GetBondWithIdx(0), Chem.CompositeQueryType.COMPOSITE_OR)
    self.assertEqual(len(mol.GetSubstructMatches(pat)), 2)

  def testGitHub1985(self):
    # simple check, this used to throw an exception
    try:
      Chem.MolToSmarts(Chem.MolFromSmarts("[C@]"))
    except Exception:
      self.fail("[C@] caused an exception when roundtripping smarts")

  def testGetEnhancedStereo(self):

    rdbase = os.environ['RDBASE']
    filename = os.path.join(rdbase, 'Code/GraphMol/FileParsers/test_data/two_centers_or.mol')
    m = Chem.MolFromMolFile(filename)

    sg = m.GetStereoGroups()
    self.assertEqual(len(sg), 2)
    group1 = sg[1]
    self.assertEqual(group1.GetGroupType(), Chem.StereoGroupType.STEREO_OR)
    stereo_atoms = group1.GetAtoms()
    self.assertEqual(len(stereo_atoms), 2)
    # file is 1 indexed and says 5
    self.assertEqual(stereo_atoms[1].GetIdx(), 4)

    # make sure the atoms are connected to the parent molecule
    stereo_atoms[1].SetProp("foo", "bar")
    self.assertTrue(m.GetAtomWithIdx(4).HasProp("foo"))

    # make sure that we can iterate over the atoms:
    for at in stereo_atoms:
      at.SetProp("foo2", "bar2")
      self.assertTrue(m.GetAtomWithIdx(at.GetIdx()).HasProp("foo2"))

  def testEnhancedStereoPreservesMol(self):
    """
    Check that the stereo group (and the atoms therein) preserve the lifetime
    of the associated mol.
    """
    rdbase = os.environ['RDBASE']
    filename = os.path.join(rdbase, 'Code/GraphMol/FileParsers/test_data/two_centers_or.mol')
    m = Chem.MolFromMolFile(filename)

    sg = m.GetStereoGroups()
    m = None
    gc.collect()
    self.assertEqual(len(sg), 2)
    group1 = sg[1]
    stereo_atoms = group1.GetAtoms()
    sg = None
    gc.collect()
    self.assertEqual(stereo_atoms[1].GetIdx(), 4)
    self.assertEqual(stereo_atoms[1].GetOwningMol().GetNumAtoms(), 8)

  def testSetEnhancedStereoGroup(self):
    m = Chem.MolFromSmiles('F[C@@H](Br)[C@H](F)Cl |o1:1|')
    m2 = Chem.RWMol(m)

    groups = m2.GetStereoGroups()
    self.assertEqual(len(groups), 1)
    # Can clear the StereoGroups by setting to an empty list
    m2.SetStereoGroups([])
    self.assertEqual(len(m2.GetStereoGroups()), 0)

    # Can add new StereoGroups
    group1 = Chem.rdchem.CreateStereoGroup(Chem.rdchem.StereoGroupType.STEREO_OR, m2, [1])
    m2.SetStereoGroups([group1])
    self.assertEqual(len(m2.GetStereoGroups()), 1)

  def testSetEnhancedStereoGroupOwnershipCheck(self):
    # make sure that the object returned by CreateStereoGroup()
    # preserves the owning molecule:
    m = Chem.RWMol(Chem.MolFromSmiles('F[C@@H](Br)[C@H](F)Cl'))
    group1 = Chem.rdchem.CreateStereoGroup(Chem.rdchem.StereoGroupType.STEREO_OR, m, [1])
    m.SetStereoGroups([group1])
    self.assertEqual(len(m.GetStereoGroups()), 1)

    m = None
    gc.collect()
    stereo_atoms = group1.GetAtoms()
    self.assertEqual(stereo_atoms[0].GetIdx(), 1)
    return
    self.assertEqual(stereo_atoms[0].GetOwningMol().GetNumAtoms(), 6)

    # make sure we can't add StereoGroups constructed from one molecule
    # to a different one:
    m2 = Chem.RWMol(Chem.MolFromSmiles('F[C@@H](Br)[C@H](F)Cl'))
    with self.assertRaises(ValueError):
      m2.SetStereoGroups([group1])

  def testSetEnhancedStereoTypeChecking(self):
    m = Chem.RWMol(Chem.MolFromSmiles('F[C@@H](Br)[C@H](F)Cl'))

    # List or tuple should be allowed:
    group = Chem.rdchem.CreateStereoGroup(Chem.rdchem.StereoGroupType.STEREO_OR, m, [1, 3])
    group = Chem.rdchem.CreateStereoGroup(Chem.rdchem.StereoGroupType.STEREO_OR, m, (1, 3))

    # Python ValueError (range error) with index past the end
    with self.assertRaises(ValueError):
      group = Chem.rdchem.CreateStereoGroup(Chem.rdchem.StereoGroupType.STEREO_OR, m, [100])

    # Mol is None
    with self.assertRaises(TypeError):
      group = Chem.rdchem.CreateStereoGroup(Chem.rdchem.StereoGroupType.STEREO_OR, None, [1])

    # Atom indices must be numbers
    with self.assertRaises(TypeError):
      group = Chem.rdchem.CreateStereoGroup(Chem.rdchem.StereoGroupType.STEREO_OR, m, [1, 'text'])

  def testSubstructParameters(self):
    m = Chem.MolFromSmiles('C[C@](F)(Cl)OCC')
    p1 = Chem.MolFromSmiles('C[C@](F)(Cl)O')
    p2 = Chem.MolFromSmiles('C[C@@](F)(Cl)O')
    p3 = Chem.MolFromSmiles('CC(F)(Cl)O')

    ps = Chem.SubstructMatchParameters()
    self.assertTrue(m.HasSubstructMatch(p1, ps))
    self.assertTrue(m.HasSubstructMatch(p2, ps))
    self.assertTrue(m.HasSubstructMatch(p3, ps))
    self.assertEqual(m.GetSubstructMatch(p1, ps), (0, 1, 2, 3, 4))
    self.assertEqual(m.GetSubstructMatch(p2, ps), (0, 1, 2, 3, 4))
    self.assertEqual(m.GetSubstructMatch(p3, ps), (0, 1, 2, 3, 4))
    self.assertEqual(m.GetSubstructMatches(p1, ps), ((0, 1, 2, 3, 4), ))
    self.assertEqual(m.GetSubstructMatches(p2, ps), ((0, 1, 2, 3, 4), ))
    self.assertEqual(m.GetSubstructMatches(p3, ps), ((0, 1, 2, 3, 4), ))
    ps.useChirality = True
    self.assertTrue(m.HasSubstructMatch(p1, ps))
    self.assertFalse(m.HasSubstructMatch(p2, ps))
    self.assertTrue(m.HasSubstructMatch(p3, ps))
    self.assertEqual(m.GetSubstructMatch(p1, ps), (0, 1, 2, 3, 4))
    self.assertEqual(m.GetSubstructMatch(p2, ps), ())
    self.assertEqual(m.GetSubstructMatch(p3, ps), (0, 1, 2, 3, 4))
    self.assertEqual(m.GetSubstructMatches(p1, ps), ((0, 1, 2, 3, 4), ))
    self.assertEqual(m.GetSubstructMatches(p2, ps), ())
    self.assertEqual(m.GetSubstructMatches(p3, ps), ((0, 1, 2, 3, 4), ))

  def testForwardEnhancedStereoGroupIds(self):
    m = Chem.MolFromSmiles('C[C@H](O)Cl |o5:1|')
    self.assertIsNotNone(m)

    # StereoGroup id is read, but not forwarded to "write id"
    stgs = m.GetStereoGroups()
    self.assertEqual(len(stgs), 1)
    self.assertEqual(stgs[0].GetGroupType(), Chem.StereoGroupType.STEREO_OR)
    self.assertEqual(stgs[0].GetReadId(), 5)
    self.assertEqual(stgs[0].GetWriteId(), 0)

    self.assertEqual(Chem.MolToCXSmiles(m), 'C[C@H](O)Cl |o1:1|')

    stgs[0].SetWriteId(7)
    self.assertEqual(stgs[0].GetWriteId(), 7)
    self.assertEqual(Chem.MolToCXSmiles(m), 'C[C@H](O)Cl |o7:1|')

    # ids are forwarded to copies of the mol
    m2 = Chem.RWMol(m)
    self.assertIsNotNone(m2)

    stgs2 = m2.GetStereoGroups()
    self.assertEqual(len(stgs), 1)
    self.assertEqual(stgs2[0].GetGroupType(), Chem.StereoGroupType.STEREO_OR)
    self.assertEqual(stgs2[0].GetReadId(), 5)
    self.assertEqual(stgs2[0].GetWriteId(), 7)

    self.assertEqual(Chem.MolToCXSmiles(m2), 'C[C@H](O)Cl |o7:1|')

    # Forwardings the ids overrides the WriteId
    Chem.ForwardStereoGroupIds(m)
    self.assertEqual(stgs[0].GetWriteId(), 5)
    self.assertEqual(Chem.MolToCXSmiles(m), 'C[C@H](O)Cl |o5:1|')

    # the copy mol is not affected
    self.assertEqual(stgs2[0].GetWriteId(), 7)

  def testSubstructParametersBundles(self):
    b = Chem.MolBundle()
    smis = ('C[C@](F)(Cl)O', 'C[C@](Br)(Cl)O', 'C[C@](I)(Cl)O')
    for smi in smis:
      b.AddMol(Chem.MolFromSmiles(smi))
    self.assertEqual(len(b), 3)
    self.assertEqual(b.Size(), 3)
    ps = Chem.SubstructMatchParameters()
    ps.useChirality = True
    self.assertTrue(Chem.MolFromSmiles('C[C@](F)(Cl)OCC').HasSubstructMatch(b, ps))
    self.assertFalse(Chem.MolFromSmiles('C[C@@](F)(Cl)OCC').HasSubstructMatch(b, ps))
    self.assertTrue(Chem.MolFromSmiles('C[C@](I)(Cl)OCC').HasSubstructMatch(b, ps))
    self.assertFalse(Chem.MolFromSmiles('C[C@@](I)(Cl)OCC').HasSubstructMatch(b, ps))

    self.assertEqual(
      Chem.MolFromSmiles('C[C@](F)(Cl)OCC').GetSubstructMatch(b, ps), (0, 1, 2, 3, 4))
    self.assertEqual(Chem.MolFromSmiles('C[C@@](F)(Cl)OCC').GetSubstructMatch(b, ps), ())
    self.assertEqual(
      Chem.MolFromSmiles('C[C@](I)(Cl)OCC').GetSubstructMatch(b, ps), (0, 1, 2, 3, 4))
    self.assertEqual(Chem.MolFromSmiles('C[C@@](I)(Cl)OCC').GetSubstructMatch(b, ps), ())

    self.assertEqual(
      Chem.MolFromSmiles('C[C@](F)(Cl)OCC').GetSubstructMatches(b, ps), ((0, 1, 2, 3, 4), ))
    self.assertEqual(Chem.MolFromSmiles('C[C@@](F)(Cl)OCC').GetSubstructMatches(b, ps), ())
    self.assertEqual(
      Chem.MolFromSmiles('C[C@](I)(Cl)OCC').GetSubstructMatches(b, ps), ((0, 1, 2, 3, 4), ))
    self.assertEqual(Chem.MolFromSmiles('C[C@@](I)(Cl)OCC').GetSubstructMatches(b, ps), ())

  def testSubstructParametersBundles2(self):
    b = Chem.MolBundle()
    smis = ('C[C@](F)(Cl)O', 'C[C@](Br)(Cl)O', 'C[C@](I)(Cl)O')
    for smi in smis:
      b.AddMol(Chem.MolFromSmiles(smi))
    self.assertEqual(len(b), 3)
    b2 = Chem.MolBundle()
    smis = ('C[C@@](F)(Cl)O', 'C[C@@](Br)(Cl)O', 'C[C@@](I)(Cl)O')
    for smi in smis:
      b2.AddMol(Chem.MolFromSmiles(smi))
    self.assertEqual(len(b2), 3)
    ps = Chem.SubstructMatchParameters()
    ps.useChirality = True
    self.assertTrue(b.HasSubstructMatch(b, ps))
    self.assertFalse(b.HasSubstructMatch(b2, ps))
    self.assertFalse(b2.HasSubstructMatch(b, ps))

    self.assertEqual(b.GetSubstructMatch(b, ps), (0, 1, 2, 3, 4))
    self.assertEqual(b.GetSubstructMatch(b2, ps), ())
    self.assertEqual(b2.GetSubstructMatch(b, ps), ())

    self.assertEqual(b.GetSubstructMatches(b, ps), ((0, 1, 2, 3, 4), ))
    self.assertEqual(b.GetSubstructMatches(b2, ps), ())
    self.assertEqual(b2.GetSubstructMatches(b, ps), ())

  def testSubstructMatchAtomProperties(self):
    m = Chem.MolFromSmiles("CCCCCCCCC")
    query = Chem.MolFromSmiles("CCC")
    m.GetAtomWithIdx(0).SetProp("test_prop", "1")
    query.GetAtomWithIdx(0).SetProp("test_prop", "1")
    ps = Chem.SubstructMatchParameters()
    ps.atomProperties = ["test_prop"]

    self.assertEqual(len(m.GetSubstructMatches(query)), 7)
    self.assertEqual(len(m.GetSubstructMatches(query, ps)), 1)

    # more than one property works as well
    m.GetAtomWithIdx(1).SetProp("test_prop2", "1")
    query.GetAtomWithIdx(1).SetProp("test_prop2", "1")
    ps.atomProperties = ["test_prop", "test_prop2"]
    self.assertEqual(len(m.GetSubstructMatches(query, ps)), 1)

  def testSubstructMatchBondProperties(self):
    m = Chem.MolFromSmiles("CCCCCCCCC")
    query = Chem.MolFromSmiles("CCC")
    m.GetBondWithIdx(0).SetProp("test_prop", "1")
    query.GetBondWithIdx(0).SetProp("test_prop", "1")
    ps = Chem.SubstructMatchParameters()
    ps.bondProperties = ["test_prop"]

    self.assertEqual(len(m.GetSubstructMatches(query)), 7)
    self.assertEqual(len(m.GetSubstructMatches(query, ps)), 1)

    # more than one property works as well
    m.GetBondWithIdx(1).SetProp("test_prop2", "1")
    query.GetBondWithIdx(1).SetProp("test_prop2", "1")
    ps.bondProperties = ["test_prop", "test_prop2"]
    self.assertEqual(len(m.GetSubstructMatches(query, ps)), 1)

    # atom and bond properties work together
    m.GetAtomWithIdx(0).SetProp("test_prop", "1")
    query.GetAtomWithIdx(0).SetProp("test_prop", "1")
    ps.atomProperties = ["test_prop"]
    self.assertEqual(len(m.GetSubstructMatches(query, ps)), 1)

  def testGithub2285(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'github2285.sdf')

    supp = Chem.ForwardSDMolSupplier(fileN, removeHs=False)
    if hasattr(supp, "__next__"):
      self.assertTrue(supp.__next__() is not None)
    else:
      self.assertTrue(supp.next() is not None)

  def testBitVectProp(self):
    bv = DataStructs.ExplicitBitVect(100)
    m = Chem.MolFromSmiles("CC")
    for atom in m.GetAtoms():
      bv.SetBit(atom.GetIdx())
      atom.SetExplicitBitVectProp("prop", bv)

    for atom in m.GetAtoms():
      bv = atom.GetExplicitBitVectProp("prop")
      self.assertTrue(bv.GetBit(atom.GetIdx()))

  def testBitVectQuery(self):
    bv = DataStructs.ExplicitBitVect(4)
    bv.SetBit(0)
    bv.SetBit(2)

    # wow, what a mouthfull..
    qa = rdqueries.HasBitVectPropWithValueQueryAtom("prop", bv, tolerance=0.0)

    m = Chem.MolFromSmiles("CC")
    for atom in m.GetAtoms():
      if atom.GetIdx() == 0:
        atom.SetExplicitBitVectProp("prop", bv)

    l = tuple([x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)])
    self.assertEqual(l, (0, ))

    m = Chem.MolFromSmiles("CC")
    for atom in m.GetAtoms():
      bv = DataStructs.ExplicitBitVect(4)
      bv.SetBit(atom.GetIdx())
      atom.SetExplicitBitVectProp("prop", bv)

    sma = Chem.MolFromSmarts("C")
    for atom in sma.GetAtoms():
      bv = DataStructs.ExplicitBitVect(4)
      bv.SetBit(1)
      qa = rdqueries.HasBitVectPropWithValueQueryAtom("prop", bv, tolerance=0.0)
      atom.ExpandQuery(qa)

    res = m.GetSubstructMatches(sma)
    self.assertEqual(res, ((1, ), ))

    sma = Chem.MolFromSmarts("C")
    for atom in sma.GetAtoms():
      bv = DataStructs.ExplicitBitVect(4)
      bv.SetBit(0)
      qa = rdqueries.HasBitVectPropWithValueQueryAtom("prop", bv, tolerance=0.0)
      atom.ExpandQuery(qa)

    res = m.GetSubstructMatches(sma)
    self.assertEqual(res, ((0, ), ))

    sma = Chem.MolFromSmarts("C")
    for atom in sma.GetAtoms():
      bv = DataStructs.ExplicitBitVect(4)
      bv.SetBit(0)
      qa = rdqueries.HasBitVectPropWithValueQueryAtom("prop", bv, tolerance=1.0)
      atom.ExpandQuery(qa)

    res = m.GetSubstructMatches(sma)
    self.assertEqual(res, ((0, ), (1, )))

  def testGithub2441(self):
    m = Chem.MolFromSmiles("CC")
    conf = Chem.Conformer(2)
    m.AddConformer(conf, assignId=False)
    m.GetConformer().SetIntProp("foo", 1)
    m.GetConformer().SetProp("bar", "foo")
    self.assertTrue(m.GetConformer().HasProp("foo"))
    self.assertFalse(m.GetConformer().HasProp("food"))
    d = m.GetConformer().GetPropsAsDict()
    self.assertTrue('foo' in d)
    self.assertTrue('bar' in d)
    self.assertEqual(d['bar'], 'foo')
    self.assertEqual(m.GetConformer().GetProp("bar"), "foo")
    self.assertEqual(m.GetConformer().GetIntProp("foo"), 1)

  def testGithub2479(self):
    # Chemistry failure in last entry
    smi2 = '''c1ccccc  duff
c1ccccc1 ok
c1ccncc1 pyridine
C(C garbage
C1CC1 ok2
C1C(Cl)C1 ok3
CC(C)(C)(C)C duff2
'''
    suppl2 = Chem.SmilesMolSupplier()
    suppl2.SetData(smi2, titleLine=False, nameColumn=1)
    l = [x for x in suppl2]
    self.assertEqual(len(l), 7)
    self.assertTrue(l[6] is None)

    # SMILES failure in last entry
    smi2 = '''c1ccccc  duff
c1ccccc1 ok
c1ccncc1 pyridine
C(C garbage
C1CC1 ok2
C1C(Cl)C1 ok3
C1C(Cl)CCCC duff2
'''
    suppl2 = Chem.SmilesMolSupplier()
    suppl2.SetData(smi2, titleLine=False, nameColumn=1)
    l = [x for x in suppl2]
    self.assertEqual(len(l), 7)
    self.assertTrue(l[6] is None)

    sdf = b"""
  Mrv1810 06051911332D

  3  2  0  0  0  0            999 V2000
  -13.3985    4.9850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -12.7066    5.4343    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -12.0654    4.9151    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
$$$$

  Mrv1810 06051911332D

  3  2  0  0  0  0            999 V2000
  -10.3083    4.8496    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.6408    5.3345    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -9.0277    4.7825    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
$$$$

  Mrv1810 06051911332D

  3  2  0  0  0  0            999 V2000
  -10.3083    4.8496    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.6"""
    suppl3 = Chem.SDMolSupplier()
    suppl3.SetData(sdf)
    l = [x for x in suppl3]
    self.assertEqual(len(l), 3)
    self.assertTrue(l[1] is None)
    self.assertTrue(l[2] is None)

    sio = BytesIO(sdf)
    suppl3 = Chem.ForwardSDMolSupplier(sio)
    l = [x for x in suppl3]
    self.assertEqual(len(l), 3)
    self.assertTrue(l[1] is None)
    self.assertTrue(l[2] is None)

    sdf = b"""
  Mrv1810 06051911332D

  3  2  0  0  0  0            999 V2000
  -13.3985    4.9850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  -12.7066    5.4343    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  -12.0654    4.9151    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
>  <pval>  (1)
[1,2,]

$$$$

  Mrv1810 06051911332D

  3  2  0  0  0  0            999 V2000
  -10.3083    4.8496    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -9.6408    5.3345    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -9.0277    4.7825    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
>  <pval>  (1)
[1,2,]
"""
    suppl3 = Chem.SDMolSupplier()
    suppl3.SetData(sdf)
    l = [x for x in suppl3]
    self.assertEqual(len(l), 2)
    self.assertTrue(l[0] is not None)
    self.assertTrue(l[1] is not None)

    sio = BytesIO(sdf)
    suppl3 = Chem.ForwardSDMolSupplier(sio)
    l = [x for x in suppl3]
    self.assertEqual(len(l), 2)
    self.assertTrue(l[0] is not None)
    self.assertTrue(l[1] is not None)

  def testCMLWriter(self):
    self.maxDiff = None  # XXX
    conf = Chem.Conformer(11)

    conf.SetAtomPosition(0, [-0.95330, 0.60416, 1.01609])
    conf.SetAtomPosition(1, [-1.00832, 1.68746, 0.83520])
    conf.SetAtomPosition(2, [-1.96274, 0.16103, 0.94471])
    conf.SetAtomPosition(3, [-0.57701, 0.44737, 2.04167])
    conf.SetAtomPosition(4, [0.00000, 0.00000, 0.00000])
    conf.SetAtomPosition(5, [-0.43038, 0.18596, -1.01377])
    conf.SetAtomPosition(6, [0.22538, -1.36531, 0.19373])
    conf.SetAtomPosition(7, [1.21993, -1.33937, 0.14580])
    conf.SetAtomPosition(8, [1.38490, 0.73003, 0.00000])
    conf.SetAtomPosition(9, [1.38490, 1.96795, 0.00000])
    conf.SetAtomPosition(10, [2.35253, -0.07700, 0.00000])

    emol = Chem.EditableMol(Chem.Mol())
    for z in [6, 1, 1, 1, 6, 1, 8, 1, 6, 8, 8]:
      emol.AddAtom(Chem.Atom(z))

    emol.AddBond(0, 1, Chem.BondType.SINGLE)
    emol.AddBond(0, 2, Chem.BondType.SINGLE)
    emol.AddBond(0, 3, Chem.BondType.SINGLE)
    emol.AddBond(0, 4, Chem.BondType.SINGLE)
    emol.AddBond(4, 5, Chem.BondType.SINGLE)
    emol.AddBond(4, 6, Chem.BondType.SINGLE)
    emol.AddBond(4, 8, Chem.BondType.SINGLE)
    emol.AddBond(6, 7, Chem.BondType.SINGLE)
    emol.AddBond(8, 9, Chem.BondType.DOUBLE)
    emol.AddBond(8, 10, Chem.BondType.SINGLE)

    mol = emol.GetMol()
    mol.SetProp('_Name', 'S-lactic acid')
    mol.AddConformer(conf)

    mol.GetAtomWithIdx(7).SetIsotope(2)
    mol.GetAtomWithIdx(10).SetFormalCharge(-1)

    mol.GetAtomWithIdx(4).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)

    cmlblock_expected = """<?xml version="1.0" encoding="utf-8"?>
<cml xmlns="http://www.xml-cml.org/schema" xmlns:convention="http://www.xml-cml.org/convention/" convention="convention:molecular">
  <molecule id="m-1" formalCharge="-1" spinMultiplicity="1">
    <name>S-lactic acid</name>
    <atomArray>
      <atom id="a0" elementType="C" formalCharge="0" hydrogenCount="3" x3="-0.953300" y3="0.604160" z3="1.016090"/>
      <atom id="a1" elementType="H" formalCharge="0" hydrogenCount="0" x3="-1.008320" y3="1.687460" z3="0.835200"/>
      <atom id="a2" elementType="H" formalCharge="0" hydrogenCount="0" x3="-1.962740" y3="0.161030" z3="0.944710"/>
      <atom id="a3" elementType="H" formalCharge="0" hydrogenCount="0" x3="-0.577010" y3="0.447370" z3="2.041670"/>
      <atom id="a4" elementType="C" formalCharge="0" hydrogenCount="1" x3="0.000000" y3="0.000000" z3="0.000000">
        <atomParity atomRefs4="a0 a5 a6 a8">1</atomParity>
      </atom>
      <atom id="a5" elementType="H" formalCharge="0" hydrogenCount="0" x3="-0.430380" y3="0.185960" z3="-1.013770"/>
      <atom id="a6" elementType="O" formalCharge="0" hydrogenCount="1" x3="0.225380" y3="-1.365310" z3="0.193730"/>
      <atom id="a7" elementType="H" formalCharge="0" hydrogenCount="0" isotopeNumber="2" x3="1.219930" y3="-1.339370" z3="0.145800"/>
      <atom id="a8" elementType="C" formalCharge="0" hydrogenCount="0" x3="1.384900" y3="0.730030" z3="0.000000"/>
      <atom id="a9" elementType="O" formalCharge="0" hydrogenCount="0" x3="1.384900" y3="1.967950" z3="0.000000"/>
      <atom id="a10" elementType="O" formalCharge="-1" hydrogenCount="0" x3="2.352530" y3="-0.077000" z3="0.000000"/>
    </atomArray>
    <bondArray>
      <bond atomRefs2="a0 a1" id="b0" order="S"/>
      <bond atomRefs2="a0 a2" id="b1" order="S"/>
      <bond atomRefs2="a0 a3" id="b2" order="S"/>
      <bond atomRefs2="a0 a4" id="b3" order="S"/>
      <bond atomRefs2="a4 a5" id="b4" order="S" bondStereo="H"/>
      <bond atomRefs2="a4 a6" id="b5" order="S"/>
      <bond atomRefs2="a4 a8" id="b6" order="S"/>
      <bond atomRefs2="a6 a7" id="b7" order="S"/>
      <bond atomRefs2="a8 a9" id="b8" order="D"/>
      <bond atomRefs2="a8 a10" id="b9" order="S"/>
    </bondArray>
  </molecule>
</cml>
"""

    self.assertEqual(Chem.MolToCMLBlock(mol), cmlblock_expected)

  def testXYZ(self):
    conf = Chem.Conformer(5)
    conf.SetAtomPosition(0, [0.000, 0.000, 0.000])
    conf.SetAtomPosition(1, [-0.635, -0.635, 0.635])
    conf.SetAtomPosition(2, [-0.635, 0.635, -0.635])
    conf.SetAtomPosition(3, [0.635, -0.635, -0.635])
    conf.SetAtomPosition(4, [0.635, 0.635, 0.635])

    emol = Chem.EditableMol(Chem.Mol())
    for z in [6, 1, 1, 1, 1]:
      emol.AddAtom(Chem.Atom(z))
    mol = emol.GetMol()
    mol.SetProp('_Name', 'methane\nthis part should not be output')
    mol.AddConformer(conf)

    xyzblock_expected = """5
methane
C      0.000000    0.000000    0.000000
H     -0.635000   -0.635000    0.635000
H     -0.635000    0.635000   -0.635000
H      0.635000   -0.635000   -0.635000
H      0.635000    0.635000    0.635000
"""

    self.assertEqual(Chem.MolToXYZBlock(mol), xyzblock_expected)

  def testSanitizationExceptionBasics(self):
    try:
      Chem.SanitizeMol(Chem.MolFromSmiles('CFC', sanitize=False))
    except Chem.AtomValenceException as exc:
      self.assertEqual(exc.cause.GetAtomIdx(), 1)
    else:
      self.assertFalse(True)

    try:
      Chem.SanitizeMol(Chem.MolFromSmiles('c1cc1', sanitize=False))
    except Chem.KekulizeException as exc:
      self.assertEqual(exc.cause.GetAtomIndices(), (0, 1, 2))
    else:
      self.assertFalse(True)

  def testSanitizationExceptionHierarchy(self):
    with self.assertRaises(Chem.AtomValenceException):
      Chem.SanitizeMol(Chem.MolFromSmiles('CFC', sanitize=False))
    with self.assertRaises(Chem.AtomSanitizeException):
      Chem.SanitizeMol(Chem.MolFromSmiles('CFC', sanitize=False))
    with self.assertRaises(Chem.MolSanitizeException):
      Chem.SanitizeMol(Chem.MolFromSmiles('CFC', sanitize=False))
    with self.assertRaises(ValueError):
      Chem.SanitizeMol(Chem.MolFromSmiles('CFC', sanitize=False))

    with self.assertRaises(Chem.KekulizeException):
      Chem.SanitizeMol(Chem.MolFromSmiles('c1cc1', sanitize=False))
    with self.assertRaises(Chem.MolSanitizeException):
      Chem.SanitizeMol(Chem.MolFromSmiles('c1cc1', sanitize=False))
    with self.assertRaises(ValueError):
      Chem.SanitizeMol(Chem.MolFromSmiles('c1cc1', sanitize=False))

  def testNoExceptionSmilesParserParams(self):
    """
    MolFromSmiles should catch exceptions even when SmilesParserParams
    is provided.
    """
    smiles_params = Chem.SmilesParserParams()
    mol = Chem.MolFromSmiles("C1CC", smiles_params)
    self.assertIsNone(mol)

  def testDetectChemistryProblems(self):
    m = Chem.MolFromSmiles('CFCc1cc1FC', sanitize=False)
    ps = Chem.DetectChemistryProblems(m)
    self.assertEqual(len(ps), 3)
    self.assertEqual([x.GetType() for x in ps],
                     ['AtomValenceException', 'AtomValenceException', 'KekulizeException'])
    self.assertEqual(ps[0].GetAtomIdx(), 1)
    self.assertEqual(ps[1].GetAtomIdx(), 6)
    self.assertEqual(ps[2].GetAtomIndices(), (3, 4, 5))

  def testGithub2611(self):
    mol = Chem.MolFromSmiles('ONCS.ONCS')
    for atom in mol.GetAtoms():
      atom.SetIsotope(atom.GetIdx())

    order1 = list(
      Chem.CanonicalRankAtomsInFragment(mol, atomsToUse=range(0, 4), breakTies=False,
                                        includeIsotopes=True))
    order2 = list(
      Chem.CanonicalRankAtomsInFragment(mol, atomsToUse=range(0, 8), breakTies=False,
                                        includeIsotopes=False))
    self.assertNotEqual(order1[:4], order2[4:])
    # ensure that the orders are ignored in the second batch
    self.assertEqual(order2[:4], order2[4:])

    for smi in ['ONCS.ONCS', 'F[C@@H](Br)[C@H](F)Cl']:
      mol = Chem.MolFromSmiles(smi)
      for atom in mol.GetAtoms():
        atom.SetIsotope(atom.GetIdx())

        for iso, chiral in [(True, True), (True, False), (False, True), (False, False)]:
          order1 = list(
            Chem.CanonicalRankAtomsInFragment(mol, atomsToUse=range(0, mol.GetNumAtoms()),
                                              bondsToUse=range(0,
                                                               mol.GetNumBonds()), breakTies=False,
                                              includeIsotopes=iso, includeChirality=chiral))
          order2 = list(
            Chem.CanonicalRankAtomsInFragment(mol, atomsToUse=range(0, mol.GetNumAtoms()),
                                              bondsToUse=range(0,
                                                               mol.GetNumBonds()), breakTies=True,
                                              includeIsotopes=iso, includeChirality=chiral))
          order3 = list(
            Chem.CanonicalRankAtoms(mol, breakTies=False, includeIsotopes=iso,
                                    includeChirality=chiral))
          order4 = list(
            Chem.CanonicalRankAtoms(mol, breakTies=True, includeIsotopes=iso,
                                    includeChirality=chiral))
          self.assertEqual(order1, order3)
          self.assertEqual(order2, order4)

  def testSetBondStereoFromDirections(self):
    m1 = Chem.MolFromMolBlock(
      '''
  Mrv1810 10141909482D

  4  3  0  0  0  0            999 V2000
    3.3412   -2.9968    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5162   -2.9968    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1037   -3.7112    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7537   -2.2823    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  1  4  1  0  0  0  0
M  END
''', sanitize=False)
    self.assertEqual(m1.GetBondBetweenAtoms(0, 1).GetBondType(), Chem.BondType.DOUBLE)
    self.assertEqual(m1.GetBondBetweenAtoms(0, 1).GetStereo(), Chem.BondStereo.STEREONONE)
    Chem.SetBondStereoFromDirections(m1)
    self.assertEqual(m1.GetBondBetweenAtoms(0, 1).GetStereo(), Chem.BondStereo.STEREOTRANS)

    m2 = Chem.MolFromMolBlock(
      '''
  Mrv1810 10141909542D

  4  3  0  0  0  0            999 V2000
    3.4745   -5.2424    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6495   -5.2424    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.2370   -5.9569    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8870   -5.9569    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  1  4  1  0  0  0  0
M  END
''', sanitize=False)
    self.assertEqual(m2.GetBondBetweenAtoms(0, 1).GetBondType(), Chem.BondType.DOUBLE)
    self.assertEqual(m2.GetBondBetweenAtoms(0, 1).GetStereo(), Chem.BondStereo.STEREONONE)
    Chem.SetBondStereoFromDirections(m2)
    self.assertEqual(m2.GetBondBetweenAtoms(0, 1).GetStereo(), Chem.BondStereo.STEREOCIS)

  def testSetBondDirFromStereo(self):
    m1 = Chem.MolFromSmiles('CC=CC')
    m1.GetBondWithIdx(1).SetStereoAtoms(0, 3)
    m1.GetBondWithIdx(1).SetStereo(Chem.BondStereo.STEREOCIS)
    Chem.SetDoubleBondNeighborDirections(m1)
    self.assertEqual(Chem.MolToSmiles(m1), r"C/C=C\C")
    self.assertEqual(m1.GetBondWithIdx(0).GetBondDir(), Chem.BondDir.ENDUPRIGHT)
    self.assertEqual(m1.GetBondWithIdx(2).GetBondDir(), Chem.BondDir.ENDDOWNRIGHT)

  def testAssignChiralTypesFromMolParity(self):

    class TestAssignChiralTypesFromMolParity:

      class BondDef:

        def __init__(self, bi, ei, t):
          self.beginIdx = bi
          self.endIdx = ei
          self.type = t

      def __init__(self, mol, parent):
        self.parent = parent
        self.parityMap = {
          Chem.ChiralType.CHI_TETRAHEDRAL_CW: 1,
          Chem.ChiralType.CHI_TETRAHEDRAL_CCW: 2,
          Chem.ChiralType.CHI_UNSPECIFIED: 0,
          Chem.ChiralType.CHI_OTHER: 0
        }
        self.d_rwMol = Chem.RWMol(mol)
        self.assignMolParity()
        self.fillBondDefVect()
        Chem.AssignAtomChiralTagsFromMolParity(self.d_rwMol)
        self.d_refSmiles = Chem.MolToSmiles(self.d_rwMol)
        self.heapPermutation()

      def assignMolParity(self):
        Chem.AssignAtomChiralTagsFromStructure(self.d_rwMol)
        for a in self.d_rwMol.GetAtoms():
          parity = self.parityMap[a.GetChiralTag()]
          a.SetIntProp("molParity", parity)
          a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)

      def fillBondDefVect(self):
        self.d_bondDefVect = [
          self.BondDef(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondType())
          for b in self.d_rwMol.GetBonds()
        ]

      def stripBonds(self):
        for i in reversed(range(self.d_rwMol.GetNumBonds())):
          b = self.d_rwMol.GetBondWithIdx(i)
          self.d_rwMol.RemoveBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx())

      def addBonds(self):
        [
          self.d_rwMol.AddBond(bondDef.beginIdx, bondDef.endIdx, bondDef.type)
          for bondDef in self.d_bondDefVect
        ]

      def checkBondPermutation(self):
        self.stripBonds()
        self.addBonds()
        Chem.SanitizeMol(self.d_rwMol)
        Chem.AssignAtomChiralTagsFromMolParity(self.d_rwMol)
        self.parent.assertEqual(Chem.MolToSmiles(self.d_rwMol), self.d_refSmiles)

      def heapPermutation(self, s=0):
        # if size becomes 1 the permutation is ready to use
        if (s == 0):
          s = len(self.d_bondDefVect)
        if (s == 1):
          self.checkBondPermutation()
          return
        for i in range(s):
          self.heapPermutation(s - 1)
          # if size is odd, swap first and last element
          j = 0 if (s % 2 == 1) else i
          self.d_bondDefVect[j], self.d_bondDefVect[s - 1] = \
            self.d_bondDefVect[s - 1], self.d_bondDefVect[j]

    molb = """
     RDKit          3D

  6  5  0  0  1  0  0  0  0  0999 V2000
   -2.9747    1.7234    0.0753 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4586    1.4435    0.1253 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5885    2.6215    1.4893 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -3.7306    0.3885   -0.0148 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3395    3.0471    0.1580 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -1.1574    0.7125    1.2684 F   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  2  5  1  0
  2  6  1  0
M  END
"""
    m = Chem.RWMol(Chem.MolFromMolBlock(molb, sanitize=True, removeHs=False))
    self.assertIsNotNone(m)
    TestAssignChiralTypesFromMolParity(m, self)

  def testCXSMILESErrors(self):
    smi = "CCC |FAILURE|"
    ps = Chem.SmilesParserParams()
    ps.strictCXSMILES = False
    m = Chem.MolFromSmiles(smi, ps)
    self.assertTrue(m is not None)
    self.assertEqual(m.GetNumAtoms(), 3)

  def testRemoveHsParams(self):
    smips = Chem.SmilesParserParams()
    smips.removeHs = False

    m = Chem.MolFromSmiles('F.[H]', smips)
    ps = Chem.RemoveHsParameters()
    m = Chem.RemoveHs(m, ps)
    self.assertEqual(m.GetNumAtoms(), 2)
    ps.removeDegreeZero = True
    m = Chem.RemoveHs(m, ps)
    self.assertEqual(m.GetNumAtoms(), 1)

    m = Chem.MolFromSmiles('F[H-]F', smips)
    ps = Chem.RemoveHsParameters()
    m = Chem.RemoveHs(m, ps)
    self.assertEqual(m.GetNumAtoms(), 3)
    m = Chem.MolFromSmiles('F[H-]F', smips)
    ps.removeHigherDegrees = True
    m = Chem.RemoveHs(m, ps)
    self.assertEqual(m.GetNumAtoms(), 2)

    m = Chem.MolFromSmiles('[H][H]', smips)
    ps = Chem.RemoveHsParameters()
    m = Chem.RemoveHs(m, ps)
    self.assertEqual(m.GetNumAtoms(), 2)
    m = Chem.MolFromSmiles('[H][H]', smips)
    ps.removeOnlyHNeighbors = True
    m = Chem.RemoveHs(m, ps)
    self.assertEqual(m.GetNumAtoms(), 0)

    m = Chem.MolFromSmiles('F[2H]', smips)
    ps = Chem.RemoveHsParameters()
    m = Chem.RemoveHs(m, ps)
    self.assertEqual(m.GetNumAtoms(), 2)
    m = Chem.MolFromSmiles('F[2H]', smips)
    ps.removeIsotopes = True
    m = Chem.RemoveHs(m, ps)
    self.assertEqual(m.GetNumAtoms(), 1)

    m = Chem.MolFromSmiles('c1c(C([2H])([2H])F)cccc1', smips)
    ps = Chem.RemoveHsParameters()
    m_noh = Chem.RemoveHs(m, ps)
    self.assertEqual(m_noh.GetNumAtoms(), m.GetNumAtoms())
    ps.removeAndTrackIsotopes = True
    m_noh = Chem.RemoveHs(m, ps)
    self.assertEqual(m_noh.GetNumAtoms(), m.GetNumAtoms() - 2)
    self.assertTrue(m_noh.GetAtomWithIdx(2).HasProp("_isotopicHs"))
    self.assertEqual(tuple(m_noh.GetAtomWithIdx(2).GetPropsAsDict().get("_isotopicHs")), (2, 2))
    m_h = Chem.AddHs(m_noh)
    self.assertFalse(m_h.GetAtomWithIdx(2).HasProp("_isotopicHs"))
    self.assertEqual(
      sum([
        1 for nbr in m_h.GetAtomWithIdx(2).GetNeighbors()
        if (nbr.GetAtomicNum() == 1 and nbr.GetIsotope())
      ]), 2)

    m = Chem.MolFromSmiles('*[H]', smips)
    ps = Chem.RemoveHsParameters()
    m = Chem.RemoveHs(m, ps)
    self.assertEqual(m.GetNumAtoms(), 2)
    m = Chem.MolFromSmiles('*[H]', smips)
    ps.removeDummyNeighbors = True
    m = Chem.RemoveHs(m, ps)
    self.assertEqual(m.GetNumAtoms(), 1)

    m = Chem.MolFromSmiles('F/C=N/[H]', smips)
    ps = Chem.RemoveHsParameters()
    m = Chem.RemoveHs(m, ps)
    self.assertEqual(m.GetNumAtoms(), 4)
    m = Chem.MolFromSmiles('F/C=N/[H]', smips)
    ps.removeDefiningBondStereo = True
    m = Chem.RemoveHs(m, ps)
    self.assertEqual(m.GetNumAtoms(), 3)

    m = Chem.MolFromSmiles('FC([H])(O)Cl', smips)
    m.GetBondBetweenAtoms(1, 2).SetBondDir(Chem.BondDir.BEGINWEDGE)
    ps = Chem.RemoveHsParameters()
    m = Chem.RemoveHs(m, ps)
    self.assertEqual(m.GetNumAtoms(), 4)
    m = Chem.MolFromSmiles('FC([H])(O)Cl', smips)
    m.GetBondBetweenAtoms(1, 2).SetBondDir(Chem.BondDir.BEGINWEDGE)
    ps.removeWithWedgedBond = False
    m = Chem.RemoveHs(m, ps)
    self.assertEqual(m.GetNumAtoms(), 5)

    m = Chem.MolFromSmarts('F[#1]')
    ps = Chem.RemoveHsParameters()
    m = Chem.RemoveHs(m, ps)
    self.assertEqual(m.GetNumAtoms(), 2)
    m = Chem.MolFromSmarts('F[#1]')
    ps.removeWithQuery = True
    m = Chem.RemoveHs(m, ps)
    self.assertEqual(m.GetNumAtoms(), 1)

    m = Chem.MolFromSmiles('[C@]12([H])CCC1CO2.[H+].F[H-]F.[H][H].[H]*.F/C=C/[H]')
    m = Chem.RemoveAllHs(m)
    for at in m.GetAtoms():
      self.assertNotEqual(at.GetAtomicNum(), 1)

  def testPickleCoordsAsDouble(self):
    import pickle
    m = Chem.MolFromSmiles('C')
    test_num = 1e50
    conf = Chem.Conformer(1)
    conf.SetAtomPosition(0, (test_num, 0.0, 0.0))
    m.AddConformer(conf)

    opts = Chem.GetDefaultPickleProperties()
    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.NoProps)
    self.assertNotEqual(pickle.loads(pickle.dumps(m)).GetConformer().GetAtomPosition(0).x, test_num)

    try:
      Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.CoordsAsDouble)
      self.assertEqual(pickle.loads(pickle.dumps(m)).GetConformer().GetAtomPosition(0).x, test_num)
    finally:
      Chem.SetDefaultPickleProperties(opts)

  def testCustomSubstructMatchCheck(self):

    def accept_none(mol, vect):
      return False

    def accept_all(mol, vect):
      return True

    def accept_large(mol, vect):
      return sum(vect) > 5

    m = Chem.MolFromSmiles('CCOCC')
    p = Chem.MolFromSmiles('CCO')
    ps = Chem.SubstructMatchParameters()
    self.assertEqual(len(m.GetSubstructMatches(p, ps)), 2)

    ps.setExtraFinalCheck(accept_none)
    self.assertEqual(len(m.GetSubstructMatches(p, ps)), 0)
    self.assertEqual(len(m.GetSubstructMatch(p, ps)), 0)
    self.assertFalse(m.HasSubstructMatch(p, ps))

    ps.setExtraFinalCheck(accept_all)
    self.assertEqual(len(m.GetSubstructMatches(p, ps)), 2)
    self.assertEqual(len(m.GetSubstructMatch(p, ps)), 3)
    self.assertTrue(m.HasSubstructMatch(p, ps))

    ps.setExtraFinalCheck(accept_large)
    self.assertEqual(len(m.GetSubstructMatches(p, ps)), 1)
    self.assertEqual(len(m.GetSubstructMatch(p, ps)), 3)
    self.assertTrue(m.HasSubstructMatch(p, ps))

  def testMostSubstitutedCoreMatch(self):
    core = Chem.MolFromSmarts("[*:1]c1cc([*:2])ccc1[*:3]")
    orthoMeta = Chem.MolFromSmiles("c1ccc(-c2ccc(-c3ccccc3)c(-c3ccccc3)c2)cc1")
    ortho = Chem.MolFromSmiles("c1ccc(-c2ccccc2-c2ccccc2)cc1")
    meta = Chem.MolFromSmiles("c1ccc(-c2cccc(-c3ccccc3)c2)cc1")
    biphenyl = Chem.MolFromSmiles("c1ccccc1-c1ccccc1")
    phenyl = Chem.MolFromSmiles("c1ccccc1")

    def numHsMatchingDummies(mol, core, match):
      return sum([
        1 for qi, ai in enumerate(match) if core.GetAtomWithIdx(qi).GetAtomicNum() == 0
        and mol.GetAtomWithIdx(ai).GetAtomicNum() == 1
      ])

    for mol, res in ((orthoMeta, 0), (ortho, 1), (meta, 1), (biphenyl, 2), (phenyl, 3)):
      mol = Chem.AddHs(mol)
      matches = mol.GetSubstructMatches(core)
      bestMatch = Chem.GetMostSubstitutedCoreMatch(mol, core, matches)
      self.assertEqual(numHsMatchingDummies(mol, core, bestMatch), res)
      ctrlCounts = sorted([numHsMatchingDummies(mol, core, match) for match in matches])
      sortedCounts = [
        numHsMatchingDummies(mol, core, match)
        for match in Chem.SortMatchesByDegreeOfCoreSubstitution(mol, core, matches)
      ]
      self.assertEqual(len(ctrlCounts), len(sortedCounts))
      self.assertTrue(all(ctrl == sortedCounts[i] for i, ctrl in enumerate(ctrlCounts)))
    with self.assertRaises(ValueError):
      Chem.GetMostSubstitutedCoreMatch(orthoMeta, core, [])
    with self.assertRaises(ValueError):
      Chem.SortMatchesByDegreeOfCoreSubstitution(orthoMeta, core, [])

  def testSetCoordsTerminalAtom(self):
    mol = Chem.MolFromMolBlock("""
     RDKit          2D

  6  6  0  0  0  0  0  0  0  0999 V2000
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
M  END
""")
    mol = Chem.RWMol(mol)
    atom = Chem.Atom(0)
    idx = mol.AddAtom(atom)
    mol.AddBond(idx, 0)
    Chem.SetTerminalAtomCoords(mol, idx, 0)
    coord = mol.GetConformer().GetAtomPosition(idx)
    self.assertAlmostEqual(coord.x, 2.5, 2)
    self.assertAlmostEqual(coord.y, 0, 2)
    self.assertAlmostEqual(coord.z, 0, 2)

  def testSuppliersReadingDirectories(self):
    # this is an odd one, basically we need to check that we don't hang
    #  which is pretty much a bad test in my opinion, but YMMV
    d = tempfile.mkdtemp()
    self.assertTrue(os.path.exists(d))

    for supplier in [
        Chem.SmilesMolSupplier,
        Chem.SDMolSupplier,
        Chem.TDTMolSupplier,
        #Chem.CompressedSDMolSupplier,
    ]:
      print("supplier:", supplier)
      with self.assertRaises(OSError):
        suppl = supplier(d)
    if hasattr(Chem, 'MaeMolSupplier'):
      with self.assertRaises(OSError):
        suppl = Chem.MaeMolSupplier(d)
    os.rmdir(d)

  def testRandomSmilesVect(self):
    m = Chem.MolFromSmiles("C1OCC1N(CO)(Cc1ccccc1NCCl)")
    v = Chem.MolToRandomSmilesVect(m, 5, randomSeed=0xf00d)
    self.assertEqual(v, [
      "c1cc(CN(C2COC2)CO)c(cc1)NCCl", "N(CCl)c1c(CN(C2COC2)CO)cccc1", "N(CCl)c1ccccc1CN(C1COC1)CO",
      "OCN(Cc1ccccc1NCCl)C1COC1", "C(N(C1COC1)Cc1c(cccc1)NCCl)O"
    ])

    v = Chem.MolToRandomSmilesVect(m, 5, randomSeed=0xf00d, allHsExplicit=True)
    self.assertEqual(v, [
      "[cH]1[cH][c]([CH2][N]([CH]2[CH2][O][CH2]2)[CH2][OH])[c]([cH][cH]1)[NH][CH2][Cl]",
      "[NH]([CH2][Cl])[c]1[c]([CH2][N]([CH]2[CH2][O][CH2]2)[CH2][OH])[cH][cH][cH][cH]1",
      "[NH]([CH2][Cl])[c]1[cH][cH][cH][cH][c]1[CH2][N]([CH]1[CH2][O][CH2]1)[CH2][OH]",
      "[OH][CH2][N]([CH2][c]1[cH][cH][cH][cH][c]1[NH][CH2][Cl])[CH]1[CH2][O][CH2]1",
      "[CH2]([N]([CH]1[CH2][O][CH2]1)[CH2][c]1[c]([cH][cH][cH][cH]1)[NH][CH2][Cl])[OH]"
    ])

  def testGithub3198(self):

    ps = Chem.SmilesParserParams()
    ps.removeHs = False
    m = Chem.MolFromSmiles('[H]OC(F).[Na]', ps)

    qa = rdqueries.AAtomQueryAtom()
    self.assertEqual(tuple(x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)), (1, 2, 3, 4))

    qa = rdqueries.AHAtomQueryAtom()
    self.assertEqual(tuple(x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)), (0, 1, 2, 3, 4))

    qa = rdqueries.QAtomQueryAtom()
    self.assertEqual(tuple(x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)), (1, 3, 4))

    qa = rdqueries.QHAtomQueryAtom()
    self.assertEqual(tuple(x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)), (0, 1, 3, 4))

    qa = rdqueries.XAtomQueryAtom()
    self.assertEqual(tuple(x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)), (3, ))

    qa = rdqueries.XHAtomQueryAtom()
    self.assertEqual(tuple(x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)), (0, 3))

    qa = rdqueries.MAtomQueryAtom()
    self.assertEqual(tuple(x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)), (4, ))

    qa = rdqueries.MHAtomQueryAtom()
    self.assertEqual(tuple(x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)), (0, 4))

  def testAdjustQueryPropertiesFiveRings(self):
    mol = Chem.MolFromSmiles('c1cc[nH]c1')
    nops = Chem.AdjustQueryParameters.NoAdjustments()
    nmol = Chem.AdjustQueryProperties(mol, nops)
    self.assertEqual(Chem.MolToSmarts(nmol), "[#6]1:[#6]:[#6]:[#7H]:[#6]:1")

    nops.adjustConjugatedFiveRings = True
    nmol = Chem.AdjustQueryProperties(mol, nops)
    self.assertEqual(Chem.MolToSmarts(nmol), "[#6]1-,=,:[#6]-,=,:[#6]-,=,:[#7H]-,=,:[#6]-,=,:1")

  def testFindPotentialStereo(self):
    mol = Chem.MolFromSmiles('C[C@H](F)C=CC')
    si = Chem.FindPotentialStereo(mol)
    self.assertEqual(len(si), 2)
    self.assertEqual(si[0].type, Chem.StereoType.Atom_Tetrahedral)
    self.assertEqual(si[0].specified, Chem.StereoSpecified.Specified)
    self.assertEqual(si[0].centeredOn, 1)
    self.assertEqual(si[0].descriptor, Chem.StereoDescriptor.Tet_CCW)
    self.assertEqual(list(si[0].controllingAtoms), [0, 2, 3])
    self.assertEqual(si[1].type, Chem.StereoType.Bond_Double)
    self.assertEqual(si[1].specified, Chem.StereoSpecified.Unspecified)
    self.assertEqual(si[1].centeredOn, 3)
    self.assertEqual(si[1].descriptor, Chem.StereoDescriptor.NoValue)
    self.assertEqual(list(si[1].controllingAtoms),
                     [1, Chem.StereoInfo.NOATOM, 5, Chem.StereoInfo.NOATOM])

  def testNewFindMolChiralCenters(self):
    mol = Chem.MolFromSmiles('C[C@H](F)C=CC(F)Cl')
    ctrs = Chem.FindMolChiralCenters(mol, useLegacyImplementation=False)
    self.assertEqual(len(ctrs), 1)
    self.assertEqual(ctrs, [(1, 'S')])
    ctrs = Chem.FindMolChiralCenters(mol, useLegacyImplementation=False, includeCIP=False)
    self.assertEqual(len(ctrs), 1)
    self.assertEqual(ctrs, [(1, 'Tet_CCW')])
    ctrs = Chem.FindMolChiralCenters(mol, useLegacyImplementation=False, includeUnassigned=True,
                                     includeCIP=False)
    self.assertEqual(len(ctrs), 2)
    self.assertEqual(ctrs, [(1, 'Tet_CCW'), (5, '?')])
    ctrs = Chem.FindMolChiralCenters(mol, useLegacyImplementation=False, includeUnassigned=True,
                                     includeCIP=True)
    self.assertEqual(len(ctrs), 2)
    self.assertEqual(ctrs, [(1, 'S'), (5, '?')])

  def testGithub6945(self):
    origVal = Chem.GetUseLegacyStereoPerception()
    tgt = [(1, '?'), (4, 'R')]
    try:
      for opt in (not origVal, origVal):
        Chem.SetUseLegacyStereoPerception(opt)
        # make sure calling with the default value works:
        self.assertEqual(
          tgt,
          Chem.FindMolChiralCenters(Chem.MolFromSmiles('FC(Cl)(Br)[C@H](F)Cl'),
                                    includeUnassigned=True))
        for useLegacy in (True, False):
          self.assertEqual(
            tgt,
            Chem.FindMolChiralCenters(Chem.MolFromSmiles('FC(Cl)(Br)[C@H](F)Cl'),
                                      includeUnassigned=True, useLegacyImplementation=useLegacy))
    except:
      raise
    finally:
      Chem.SetUseLegacyStereoPerception(origVal)

  @unittest.skipUnless(hasattr(Chem, 'MolFromPNGFile'), "RDKit not built with iostreams support")
  def testMolFromPNG(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'colchicine.png')
    mol = Chem.MolFromPNGFile(fileN)
    self.assertIsNotNone(mol)
    self.assertEqual(mol.GetNumAtoms(), 29)

    with open(fileN, 'rb') as inf:
      d = inf.read()
    mol = Chem.MolFromPNGString(d)
    self.assertIsNotNone(mol)
    self.assertEqual(mol.GetNumAtoms(), 29)

  @unittest.skipUnless(hasattr(Chem, 'MolFromPNGFile'), "RDKit not built with iostreams support")
  def testMolToPNG(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'colchicine.no_metadata.png')

    with open(fileN, 'rb') as inf:
      d = inf.read()
    mol = Chem.MolFromSmiles("COc1cc2c(c(OC)c1OC)-c1ccc(OC)c(=O)cc1[C@@H](NC(C)=O)CC2")
    self.assertIsNotNone(mol)
    self.assertEqual(mol.GetNumAtoms(), 29)

    nd = Chem.MolMetadataToPNGString(mol, d)
    mol = Chem.MolFromPNGString(nd)
    self.assertIsNotNone(mol)
    self.assertEqual(mol.GetNumAtoms(), 29)

    nd = Chem.MolMetadataToPNGFile(mol, fileN)
    mol = Chem.MolFromPNGString(nd)
    self.assertIsNotNone(mol)
    self.assertEqual(mol.GetNumAtoms(), 29)

  @unittest.skipUnless(hasattr(Chem, 'MolFromPNGFile'), "RDKit not built with iostreams support")
  def testMolsFromPNG(self):
    refMols = [Chem.MolFromSmiles(x) for x in ('c1ccccc1', 'CCO', 'CC(=O)O', 'c1ccccn1')]
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'multiple_mols.png')
    mols = Chem.MolsFromPNGFile(fileN)
    self.assertEqual(len(mols), len(refMols))
    for mol, refMol in zip(mols, refMols):
      self.assertEqual(Chem.MolToSmiles(mol), Chem.MolToSmiles(refMol))

  @unittest.skipUnless(hasattr(Chem, 'MolFromPNGFile'), "RDKit not built with iostreams support")
  def testMetadataToPNG(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'colchicine.png')

    with open(fileN, 'rb') as inf:
      d = inf.read()
    mol = Chem.MolFromPNGString(d)
    nd = Chem.MolMetadataToPNGString(mol, d)
    vals = {'foo': '1', 'bar': '2'}
    nd = Chem.AddMetadataToPNGString(vals, nd)
    nvals = Chem.MetadataFromPNGString(nd)
    self.assertTrue('foo' in nvals)
    self.assertEqual(nvals['foo'], b'1')
    self.assertTrue('bar' in nvals)
    self.assertEqual(nvals['bar'], b'2')

    nd = Chem.AddMetadataToPNGFile(vals, fileN)
    nvals = Chem.MetadataFromPNGString(nd)
    self.assertTrue('foo' in nvals)
    self.assertEqual(nvals['foo'], b'1')
    self.assertTrue('bar' in nvals)
    self.assertEqual(nvals['bar'], b'2')

    vals = {'foo': 1, 'bar': '2'}
    with self.assertRaises(TypeError):
      nd = Chem.AddMetadataToPNGString(vals, d)

  def test_github3403(self):
    core1 = "[$(C-!@[a])](=O)(Cl)"
    sma = Chem.MolFromSmarts(core1)

    m = Chem.MolFromSmiles("c1ccccc1C(=O)Cl")
    self.assertFalse(m.HasSubstructMatch(sma, recursionPossible=False))

    m = Chem.MolFromSmiles("c1ccccc1C(=O)Cl")
    self.assertTrue(m.HasSubstructMatch(sma))

    m = Chem.MolFromSmiles("c1ccccc1C(=O)Cl")
    self.assertFalse(m.HasSubstructMatch(sma, recursionPossible=False))

  def test_github3517(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'NCI_aids_few.sdf')
    sdSup = Chem.SDMolSupplier(fileN)
    mols_list = list(sdSup)
    mols_list_compr = [m for m in sdSup]
    self.assertEqual(len(mols_list), len(mols_list_compr))

  def test_github3492(self):

    def read_smile(s):
      m = Chem.MolFromSmiles(s)
      rdkit.Chem.rdDepictor.Compute2DCoords(m)
      return m

    def sq_dist(a, b):
      ab = [a[i] - b[i] for i, _ in enumerate(a)]
      return sum([d * d for d in ab])

    self.assertIsNotNone(Chem.MolFromSmiles("OCCN").GetAtoms()[0].GetOwningMol())
    self.assertEqual([Chem.MolFromSmiles("OCCN").GetAtoms()[i].GetAtomicNum() for i in range(4)],
                     [8, 6, 6, 7])
    self.assertIsNotNone(Chem.MolFromSmiles("O=CCC=N").GetBonds()[0].GetOwningMol())
    self.assertEqual(
      [Chem.MolFromSmiles("O=CCC=N").GetBonds()[i].GetBondType() for i in range(4)],
      [Chem.BondType.DOUBLE, Chem.BondType.SINGLE, Chem.BondType.SINGLE, Chem.BondType.DOUBLE])
    self.assertIsNotNone(read_smile("CCC").GetConformers()[0].GetOwningMol())
    pos = read_smile("CCC").GetConformers()[0].GetPositions()
    self.assertAlmostEqual(sq_dist(pos[0], pos[1]), sq_dist(pos[1], pos[2]))

  def test_github3553(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Wrap', 'test_data',
                         'github3553.sdf')
    sdSup = Chem.SDMolSupplier(fileN)
    for mol in sdSup:
      pval = mol.GetProp('boiling.point.predicted')
      sio = StringIO()
      w = Chem.SDWriter(sio)
      w.SetKekulize(True)
      w.SetForceV3000(True)
      w.write(mol)
      w.flush()
      txt = sio.getvalue()
      self.assertTrue(pval in txt)

  def test_github1631(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         '1CRN.pdb')

    m = Chem.MolFromPDBFile(fileN)
    info = m.GetAtomWithIdx(0).GetPDBResidueInfo()
    self.assertEqual(info.GetName(), " N  ")
    self.assertEqual(info.GetResidueName(), "THR")
    self.assertAlmostEqual(info.GetTempFactor(), 13.79, 2)

    m2 = Chem.MolFromSmiles('CC(C(C(=O)O)N)O')
    self.assertTrue(m2.GetAtomWithIdx(6).GetPDBResidueInfo() is None)
    m2.GetAtomWithIdx(6).SetPDBResidueInfo(info)
    info2 = m2.GetAtomWithIdx(6).GetPDBResidueInfo()
    self.assertEqual(info2.GetName(), " N  ")
    self.assertEqual(info2.GetResidueName(), "THR")
    self.assertAlmostEqual(info2.GetTempFactor(), 13.79, 2)

  def testMolzip(self):
    tests = [["C[*:1]", "N[*:1]", "CN", Chem.MolzipParams()]]
    for a, b, res, params in tests:
      self.assertEqual(
        Chem.CanonSmiles(res),
        Chem.MolToSmiles(Chem.molzip(Chem.MolFromSmiles(a), Chem.MolFromSmiles(b), params)))

    # multiple arg test
    a = Chem.MolFromSmiles('C=C[1*]')
    b = Chem.MolFromSmiles('O/C=N/[1*]')
    p = Chem.MolzipParams()
    p.label = Chem.MolzipLabel.Isotope
    c = Chem.molzip(a, b, p)
    self.assertEqual(Chem.MolToSmiles(c), 'C=C/N=C/O')

    # single argument test
    a = Chem.MolFromSmiles('C=C[1*].O/C=N/[1*]')
    p = Chem.MolzipParams()
    p.label = Chem.MolzipLabel.Isotope
    c = Chem.molzip(a, p)
    self.assertEqual(Chem.MolToSmiles(c), 'C=C/N=C/O')

    a = Chem.MolFromSmiles("[C@H]([Xe])(F)([V])")
    b = Chem.MolFromSmiles("[Xe]N.[V]I")
    p = Chem.MolzipParams()
    p.label = Chem.MolzipLabel.AtomType
    p.setAtomSymbols(["Xe", "V"])
    c = Chem.molzip(a, b, p)
    self.assertEqual(Chem.MolToSmiles(c), "N[C@@H](F)I")

    a = Chem.MolFromSmiles("C(=[*:1])N")
    b = Chem.MolFromSmiles("[*:1]-N=C")
    p = Chem.MolzipParams()
    p.enforceValenceRules = False
    c = Chem.molzip(a, b, p)

  def testContextManagers(self):
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Wrap', 'test_data',
                         'github3553.sdf')
    with Chem.SDMolSupplier(fileN) as suppl:
      mols = [x for x in suppl if x is not None]
    sio = StringIO()
    with Chem.SDWriter(sio) as w:
      for m in mols:
        w.write(m)
    txt = sio.getvalue()
    self.assertEqual(txt.count('$$$$'), len(mols))
    with self.assertRaises(RuntimeError):
      w.write(mols[0])

    with Chem.ForwardSDMolSupplier(fileN) as suppl:
      mols = [x for x in suppl if x is not None]

    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'first_200.tpsa.csv')
    with Chem.SmilesMolSupplier(fileN, ",", 0, -1) as suppl:
      ms = [x for x in suppl if x is not None]

    sio = StringIO()
    with Chem.SmilesWriter(sio) as w:
      for m in mols:
        w.write(m)
    txt = sio.getvalue()
    self.assertEqual(txt.count('\n'), len(mols) + 1)
    with self.assertRaises(RuntimeError):
      w.write(mols[0])

    data = """$SMI<Cc1nnc(N)nc1C>
CAS<17584-12-2>
|
$SMI<Cc1n[nH]c(=O)nc1N>
CAS<~>
|
$SMI<Cc1n[nH]c(=O)[nH]c1=O>
CAS<932-53-6>
|
$SMI<Cc1nnc(NN)nc1O>
CAS<~>
|"""
    with Chem.TDTMolSupplier() as suppl:
      suppl.SetData(data, "CAS")
      ms = [x for x in suppl if x is not None]
    sio = StringIO()
    with Chem.TDTWriter(sio) as w:
      for m in mols:
        w.write(m)
    txt = sio.getvalue()
    self.assertEqual(txt.count('$SMI'), len(mols))
    with self.assertRaises(RuntimeError):
      w.write(mols[0])

    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         '1CRN.pdb')
    m = Chem.MolFromPDBFile(fileN)
    mols = [m, m]
    sio = StringIO()
    with Chem.PDBWriter(sio) as w:
      for m in mols:
        w.write(m)
    txt = sio.getvalue()
    self.assertEqual(txt.count('COMPND    CRAMBIN'), len(mols))
    with self.assertRaises(RuntimeError):
      w.write(mols[0])

    if hasattr(Chem, 'MaeMolSupplier'):
      fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                           'NCI_aids_few.mae')
      with Chem.MaeMolSupplier(fileN) as suppl:
        ms = [x for x in suppl if x is not None]

    RDLogger.EnableLog('rdApp.*')

  def testInsertMol(self):
    m = Chem.MolFromSmiles("CNO")
    m2 = Chem.MolFromSmiles("c1ccccc1")
    m3 = Chem.MolFromSmiles("C1CC1")

    rwmol = Chem.RWMol(m)
    rwmol.InsertMol(m2)
    rwmol.InsertMol(m3)
    self.assertEqual(Chem.MolToSmiles(rwmol), Chem.CanonSmiles("CNO.c1ccccc1.C1CC1"))

  def testBatchEdits(self):
    mol = Chem.MolFromSmiles("C1CCCO1")

    for rwmol in [Chem.EditableMol(mol), Chem.RWMol(mol)]:
      rwmol.BeginBatchEdit()
      rwmol.RemoveAtom(2)
      rwmol.RemoveAtom(3)
      rwmol.CommitBatchEdit()
      nmol = rwmol.GetMol()
      self.assertEqual(Chem.MolToSmiles(nmol), "CCO")

    for rwmol in [Chem.EditableMol(mol), Chem.RWMol(mol)]:
      rwmol = Chem.EditableMol(mol)
      rwmol.BeginBatchEdit()
      rwmol.RemoveAtom(3)
      rwmol.RemoveBond(4, 0)
      rwmol.CommitBatchEdit()
      nmol = rwmol.GetMol()
      self.assertEqual(Chem.MolToSmiles(nmol), "CCC.O")

    for rwmol in [Chem.EditableMol(mol), Chem.RWMol(mol)]:
      rwmol.BeginBatchEdit()
      rwmol.RemoveAtom(2)
      rwmol.RemoveAtom(3)
      rwmol.RollbackBatchEdit()
      nmol = rwmol.GetMol()
      self.assertEqual(Chem.MolToSmiles(nmol), "C1CCOC1")

  def testBatchEditContextManager(self):
    mol = Chem.MolFromSmiles("C1CCCO1")
    with Chem.RWMol(mol) as rwmol:
      rwmol.RemoveAtom(2)
      rwmol.RemoveAtom(3)
      # make sure we haven't actually changed anything yet:
      self.assertEqual(rwmol.GetNumAtoms(), mol.GetNumAtoms())
    self.assertEqual(rwmol.GetNumAtoms(), mol.GetNumAtoms() - 2)

    # make sure no changes get made if we throw an exception
    try:
      with Chem.RWMol(mol) as rwmol:
        rwmol.RemoveAtom(2)
        rwmol.RemoveAtom(6)
    except Exception:
      pass
    self.assertEqual(rwmol.GetNumAtoms(), mol.GetNumAtoms())

  def testGithub4138(self):
    m = Chem.MolFromSmiles('C1CCCO1')
    q = Chem.MolFromSmarts('')
    self.assertFalse(m.HasSubstructMatch(q))
    self.assertEqual(m.GetSubstructMatch(q), ())
    self.assertEqual(m.GetSubstructMatches(q), ())

    m = Chem.MolFromSmiles('')
    q = Chem.MolFromSmarts('C')
    self.assertFalse(m.HasSubstructMatch(q))
    self.assertEqual(m.GetSubstructMatch(q), ())
    self.assertEqual(m.GetSubstructMatches(q), ())

  def testGithub4144(self):
    ''' the underlying problem with #4144 was that the
    includeRings argument could not be passed to ClearComputedProps()
    from Python. Make sure that's fixed
    '''
    m = Chem.MolFromSmiles('c1ccccc1')
    self.assertEqual(m.GetRingInfo().NumRings(), 1)
    m.ClearComputedProps(includeRings=False)
    self.assertEqual(m.GetRingInfo().NumRings(), 1)
    m.ClearComputedProps()
    with self.assertRaises(RuntimeError):
      m.GetRingInfo().NumRings()

  def testAddWavyBondsForStereoAny(self):
    m = Chem.MolFromSmiles('CC=CC')
    m.GetBondWithIdx(1).SetStereo(Chem.BondStereo.STEREOANY)
    m2 = Chem.Mol(m)
    Chem.AddWavyBondsForStereoAny(m2)
    self.assertEqual(m2.GetBondWithIdx(1).GetStereo(), Chem.BondStereo.STEREONONE)
    self.assertEqual(m2.GetBondWithIdx(0).GetBondDir(), Chem.BondDir.UNKNOWN)
    m2 = Chem.Mol(m)
    Chem.AddWavyBondsForStereoAny(m2, clearDoubleBondFlags=False)
    self.assertEqual(m2.GetBondWithIdx(1).GetStereo(), Chem.BondStereo.STEREOANY)
    self.assertEqual(m2.GetBondWithIdx(0).GetBondDir(), Chem.BondDir.UNKNOWN)

    m = Chem.MolFromSmiles("CC=CC=CC=CC")
    m.GetBondWithIdx(1).SetStereoAtoms(0, 3)
    m.GetBondWithIdx(1).SetStereo(Chem.BondStereo.STEREOCIS)
    m.GetBondWithIdx(3).SetStereo(Chem.BondStereo.STEREOANY)
    m.GetBondWithIdx(5).SetStereoAtoms(4, 7)
    m.GetBondWithIdx(5).SetStereo(Chem.BondStereo.STEREOCIS)
    m2 = Chem.Mol(m)
    Chem.AddWavyBondsForStereoAny(m2)
    self.assertEqual(m2.GetBondWithIdx(3).GetStereo(), Chem.BondStereo.STEREOANY)
    self.assertEqual(m2.GetBondWithIdx(0).GetBondDir(), Chem.BondDir.NONE)

    Chem.AddWavyBondsForStereoAny(
      m2, addWhenImpossible=Chem.StereoBondThresholds.DBL_BOND_SPECIFIED_STEREO)
    self.assertEqual(m2.GetBondWithIdx(3).GetStereo(), Chem.BondStereo.STEREONONE)
    self.assertEqual(m2.GetBondWithIdx(2).GetBondDir(), Chem.BondDir.UNKNOWN)

  def testSmartsParseParams(self):
    smi = "CCC |$foo;;bar$| ourname"
    m = Chem.MolFromSmarts(smi)
    self.assertTrue(m is not None)
    ps = Chem.SmartsParserParams()
    ps.allowCXSMILES = False
    ps.parseName = False
    m = Chem.MolFromSmarts(smi, ps)
    self.assertTrue(m is None)
    ps.allowCXSMILES = True
    ps.parseName = True
    m = Chem.MolFromSmarts(smi, ps)
    self.assertTrue(m is not None)
    self.assertTrue(m.GetAtomWithIdx(0).HasProp('atomLabel'))
    self.assertEqual(m.GetAtomWithIdx(0).GetProp('atomLabel'), "foo")
    self.assertTrue(m.HasProp('_Name'))
    self.assertEqual(m.GetProp('_Name'), "ourname")
    self.assertEqual(m.GetProp("_CXSMILES_Data"), "|$foo;;bar$|")

  def testSmilesWriteParams(self):
    m = Chem.MolFromSmiles('C[C@H](F)Cl')
    ps = Chem.SmilesWriteParams()
    ps.rootedAtAtom = 1
    self.assertEqual(Chem.MolToSmiles(m, ps), "[C@@H](C)(F)Cl")
    self.assertEqual(Chem.MolToCXSmiles(m, ps), "[C@@H](C)(F)Cl")

  def testCXSmilesOptions(self):
    # just checking that the fields parameter gets used.
    # we've tested the individual values on the C++ side
    m = Chem.MolFromSmiles("OC1CCC(F)C1 |LN:1:1.3.2.6|")
    ps = Chem.SmilesWriteParams()
    ps.rootedAtAtom = 1
    self.assertEqual(
      Chem.MolToCXSmiles(m, ps, (Chem.CXSmilesFields.CX_ALL ^ Chem.CXSmilesFields.CX_LINKNODES)),
      "C1(O)CCC(F)C1")
    self.assertTrue(hasattr(Chem.CXSmilesFields, 'CX_BOND_CFG'))
    self.assertTrue(hasattr(Chem.CXSmilesFields, 'CX_ALL_BUT_COORDS'))

  def testKekulizeIfPossible(self):
    m = Chem.MolFromSmiles('c1cccn1', sanitize=False)
    m.UpdatePropertyCache(strict=False)
    Chem.KekulizeIfPossible(m)
    for atom in m.GetAtoms():
      self.assertTrue(atom.GetIsAromatic())
    for bond in m.GetBonds():
      self.assertTrue(bond.GetIsAromatic())
      self.assertEqual(bond.GetBondType(), Chem.BondType.AROMATIC)

  def testHasValenceViolation(self):
    mol = Chem.MolFromSmiles('C(C)(C)(C)(C)C', sanitize=False)
    mol.UpdatePropertyCache(strict=False)
    self.assertTrue(any(a.HasValenceViolation() for a in mol.GetAtoms()))

  def testgithub4992(self):
    if not hasattr(Chem, "Chem.MultithreadedSDMolSupplier"):
      return

    good1 = Chem.MolFromSmiles('C')
    #good1.SetProp('molname', 'good1')
    good2 = Chem.MolFromSmiles('CC')
    #good2.SetProp('molname', 'good2')
    good3 = Chem.MolFromSmiles('CCC')
    #good3.SetProp('molname', 'good3')
    bad = Chem.MolFromSmiles('CN(C)(C)C', sanitize=False)
    #bad.SetProp('molname', 'bad')

    with Chem.SDWriter("good1_good2_good3.sdf") as w:
      w.write(good1)
      w.write(good2)
      w.write(good3)
      w.write(good1)
      w.write(good2)
      w.write(good3)

    with Chem.SDWriter("good1_good2_good3_bad.sdf") as w:
      w.write(good1)
      w.write(good2)
      w.write(good3)
      w.write(bad)

    with Chem.SDWriter("good1_good2_bad_good3.sdf") as w:
      w.write(good1)
      w.write(good2)
      w.write(bad)
      w.write(good3)

    with Chem.SDWriter("bad_good1_good2_good3.sdf") as w:
      w.write(bad)
      w.write(good1)
      w.write(good2)
      w.write(good3)

    def read_mols(supplier, filename):
      i = 0
      with supplier(filename) as sdSuppl:
        count = -1
        for count, mol in enumerate(sdSuppl):
          if mol is not None:
            i += 1
      return i

    counts = []
    for i, s in enumerate((Chem.SDMolSupplier, Chem.MultithreadedSDMolSupplier)):
      for j, f in enumerate(('good1_good2_good3.sdf', 'good1_good2_bad_good3.sdf',
                             'good1_good2_bad_good3.sdf', 'bad_good1_good2_good3.sdf')):
        print(f'---------------\n{s.__name__} {f}')
        if i == 0:
          counts.append(read_mols(s, f))
        else:
          self.assertEqual(counts[j], read_mols(s, f))

  def test_blocklogs(self):
    with capture_logging(logging.INFO) as log_stream:

      # Logging is known to be problematic with static linked libs
      # so let's skip the test if we can't see the expected error
      # log before we block logging.
      Chem.MolFromSmiles('garbage_0')
      if 'garbage_0' not in log_stream.getvalue():
        self.skipTest('cannot fetch log msgs')
      else:
        log_stream.truncate(0)

      with rdBase.BlockLogs():
        Chem.MolFromSmiles('garbage_1')
      self.assertNotIn('garbage_1', log_stream.getvalue())

      Chem.MolFromSmiles('garbage_2')
      self.assertIn('garbage_2', log_stream.getvalue())
      log_stream.truncate(0)

      block = rdBase.BlockLogs()
      self.assertIsNotNone(block)

      Chem.MolFromSmiles('garbage_3')
      self.assertNotIn('garbage_3', log_stream.getvalue())

      del block

      Chem.MolFromSmiles('garbage_4')
      self.assertIn('garbage_4', log_stream.getvalue())
      log_stream.truncate(0)

  def testDisableNontetrahedralStereo(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'test_data',
                         'nontetrahedral_3d.sdf')
    origVal = Chem.GetAllowNontetrahedralChirality()
    Chem.SetAllowNontetrahedralChirality(True)
    suppl = Chem.SDMolSupplier(fileN, sanitize=False)
    for mol in suppl:
      Chem.AssignStereochemistryFrom3D(mol)
      ct = mol.GetProp("ChiralType")
      at = mol.GetAtomWithIdx(0)
      if ct == "SP":
        self.assertEqual(at.GetChiralTag(), Chem.ChiralType.CHI_SQUAREPLANAR)
      elif ct == "TB":
        self.assertEqual(at.GetChiralTag(), Chem.ChiralType.CHI_TRIGONALBIPYRAMIDAL)
      elif ct == "OH":
        self.assertEqual(at.GetChiralTag(), Chem.ChiralType.CHI_OCTAHEDRAL)
      elif ct == "TH":
        self.assertEqual(at.GetChiralTag(), Chem.ChiralType.CHI_TETRAHEDRAL)
    Chem.SetAllowNontetrahedralChirality(False)
    suppl = Chem.SDMolSupplier(fileN, sanitize=False)
    for mol in suppl:
      Chem.AssignStereochemistryFrom3D(mol)
      ct = mol.GetProp("ChiralType")
      at = mol.GetAtomWithIdx(0)
      if ct == "TH":
        self.assertEqual(at.GetChiralTag(), Chem.ChiralType.CHI_TETRAHEDRAL)
      else:
        self.assertEqual(at.GetChiralTag(), Chem.ChiralType.CHI_UNSPECIFIED)
    Chem.SetAllowNontetrahedralChirality(origVal)

  def test_legacyStereochemGlobal(self):
    origVal = Chem.GetUseLegacyStereoPerception()
    Chem.SetUseLegacyStereoPerception(True)
    m = Chem.MolFromSmiles("C[C@H]1CCC2(CC1)CC[C@H](C)C(C)C2")
    self.assertEqual(m.GetAtomWithIdx(1).GetChiralTag(), Chem.ChiralType.CHI_UNSPECIFIED)
    self.assertNotEqual(m.GetAtomWithIdx(9).GetChiralTag(), Chem.ChiralType.CHI_UNSPECIFIED)

    Chem.SetUseLegacyStereoPerception(False)
    m = Chem.MolFromSmiles("C[C@H]1CCC2(CC1)CC[C@H](C)C(C)C2")
    self.assertEqual(m.GetAtomWithIdx(1).GetChiralTag(), Chem.ChiralType.CHI_UNSPECIFIED)
    self.assertNotEqual(m.GetAtomWithIdx(9).GetChiralTag(), Chem.ChiralType.CHI_UNSPECIFIED)

    Chem.SetUseLegacyStereoPerception(origVal)

  def testValidationMrvRxn(self):
    Chem.SetUseLegacyStereoPerception(False)

    mrvBlock = ""
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Wrap', 'test_data', 'DoubleBondRxn.mrv')
    with open(fileN) as f:
      mrvBlock = f.read()
    rxn = AllChem.ReactionFromMrvBlock(mrvBlock, True, True)
    
    smi = AllChem.ReactionToSmiles(rxn)
    self.assertEqual(smi,"CC=C1COC(C)OC1>>CC=C1COC(C)OC1")
    
    rxn = AllChem.ReactionFromMrvFile(fileN, True, True)
    smi = AllChem.ReactionToSmiles(rxn)
    self.assertEqual(smi,"CC=C1COC(C)OC1>>CC=C1COC(C)OC1")

   
    Chem.SetUseLegacyStereoPerception(True)

  def testValidationMrvRxn2(self):
    Chem.SetUseLegacyStereoPerception(False)

    mrvBlock = ""
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Wrap', 'test_data', 'DoubleBondRxn2.mrv')

    with open(fileN) as f:
      mrvBlock = f.read()
    rxn = AllChem.ReactionFromMrvBlock(mrvBlock, True,True)
    
    smi = AllChem.ReactionToSmiles(rxn)
    self.assertEqual(smi,"C=C1CC(C)(C(C)(C)C)CC(=C)C1=Cc1cc(C)c(C)c(CCC)c1>>C=C1CC(C)(C(C)(C)C)CC(=C)C1=Cc1cc(C)c(Cl)c(CCC)c1")

    Chem.SetUseLegacyStereoPerception(True)

  def testValidationMrvRxn3(self):
    Chem.SetUseLegacyStereoPerception(False)

    mrvBlock = ""
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Wrap', 'test_data', 'DoubleBondAndChiralRxn.mrv')
    with open(fileN) as f:
      mrvBlock = f.read()
    rxn = AllChem.ReactionFromMrvBlock(mrvBlock, True, True)

    smi = AllChem.ReactionToSmiles(rxn)
    self.assertEqual(smi,"C=C1C[C@](C)(C(C)(C)C)CC(=C)/C1=C/c1cc(C)c(C)c(CCC)c1>>C=C1CC(C)(C(C)(C)C)CC(=C)C1=Cc1cc(C)c(Cl)c(C[C@H](C)Cl)c1")

   
    Chem.SetUseLegacyStereoPerception(True)

  def testValidationRxn(self):
    Chem.SetUseLegacyStereoPerception(False)

    rxnBlock = ""
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Wrap', 'test_data', 'DoubleBondAndChiralRxn.rxn')
    with open(fileN) as f:
      rxnBlock = f.read()

    rxn = AllChem.ReactionFromRxnBlock(rxnBlock, True,  True,  False)
    smi = AllChem.ReactionToSmiles(rxn)
    self.assertEqual(smi,"C=C1C[C@](C)(C(C)(C)C)CC(=C)/C1=C/c1cc(C)c(C)c(CCC)c1>>C=C1CC(C)(C(C)(C)C)CC(=C)C1=Cc1cc(C)c(Cl)c(C[C@H](C)Cl)c1")
  
    rxn = AllChem.ReactionFromRxnFile(fileN, True, True,  False)

    smi = AllChem.ReactionToSmiles(rxn)
    self.assertEqual(smi,"C=C1C[C@](C)(C(C)(C)C)CC(=C)/C1=C/c1cc(C)c(C)c(CCC)c1>>C=C1CC(C)(C(C)(C)C)CC(=C)C1=Cc1cc(C)c(Cl)c(C[C@H](C)Cl)c1")

   
    Chem.SetUseLegacyStereoPerception(True)


  def testValidationSmiles(self):
    Chem.SetUseLegacyStereoPerception(False)

    smi = "C=C1CC(C)(C(C)(C)C)CC(=C)/C1=C/c1cc(C)c(Cl)c(C[C@H](C)Cl)c1"
    ps = Chem.SmilesParserParams()
    ps.allowCXSMILES = True
    ps.parseName = False
    ps.sanitize = True
    ps.removeHs = False
    mol = Chem.MolFromSmiles(smi, ps)
    self.assertIsNotNone(mol)

    of = Chem.SmilesWriteParams()
    outSmi =       Chem.MolToCXSmiles(mol, of, Chem.CXSmilesFields.CX_ALL)

    self.assertEqual(outSmi,"C=C1CC(C)(C(C)(C)C)CC(=C)C1=Cc1cc(C)c(Cl)c(C[C@H](C)Cl)c1")

    smi = "C=C1C[C@](C)(C(C)(C)C)CC(=C)/C1=C/c1cc(C)c(Cl)c(C[C@H](C)Cl)c1"
    mol = Chem.MolFromSmiles(smi, ps)
    self.assertIsNotNone(mol)


    outSmi =       Chem.MolToCXSmiles(mol, of, Chem.CXSmilesFields.CX_ALL),

    self.assertEqual(outSmi[0],"C=C1C[C@@](C)(C(C)(C)C)CC(=C)/C1=C\\c1cc(C)c(Cl)c(C[C@H](C)Cl)c1")

    Chem.SetUseLegacyStereoPerception(True)

  def testValidationMol(self):
    Chem.SetUseLegacyStereoPerception(False)
    molBlock = ""
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Wrap', 'test_data',  'DoubleBond2000.mol')
    with open(fileN) as f:
      molBlock = f.read()
    mol =Chem.MolFromMolBlock(molBlock, True, True, False)
    
    smi = Chem.MolToSmiles(mol)
    self.assertEqual(smi,"C=C1CC(C)(C(C)(C)C)CC(=C)C1=Cc1cc(C)c(Cl)c(C[C@H](C)Cl)c1")
    
    rxn = Chem.MolFromMolFile(fileN, True,  True,  False)
    
    smi = Chem.MolToSmiles(rxn)
    self.assertEqual(smi,"C=C1CC(C)(C(C)(C)C)CC(=C)C1=Cc1cc(C)c(Cl)c(C[C@H](C)Cl)c1")

    # nowthe V3000 version

    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Wrap', 'test_data',  'DoubleBond3000.mol')
    with open(fileN) as f:
      molBlock = f.read()
    mol =Chem.MolFromMolBlock(molBlock, True, True, False)
    
    smi = Chem.MolToSmiles(mol)
    self.assertEqual(smi,"C=C1CC(C)(C(C)(C)C)CC(=C)C1=Cc1cc(C)c(Cl)c(C[C@H](C)Cl)c1")
    
    mol = Chem.MolFromMolFile(fileN, True, True,  False)
    smi = Chem.MolToSmiles(mol)
    self.assertEqual(smi,"C=C1CC(C)(C(C)(C)C)CC(=C)C1=Cc1cc(C)c(Cl)c(C[C@H](C)Cl)c1")

   # now the chiral versions - the double bond should be chiral

    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Wrap', 'test_data',  'DoubleBond2000Chiral.mol')
    with open(fileN) as f:
      molBlock = f.read()
    mol =Chem.MolFromMolBlock(molBlock, True, True, False)
    
    smi = Chem.MolToSmiles(mol)
    self.assertEqual(smi,"C=C1C[C@@](C)(C(C)(C)C)CC(=C)/C1=C/c1cc(C)c(Cl)c(C[C@H](C)Cl)c1")
    
    rxn = Chem.MolFromMolFile(fileN, True,  True,  False)
    
    smi = Chem.MolToSmiles(rxn)
    self.assertEqual(smi,"C=C1C[C@@](C)(C(C)(C)C)CC(=C)/C1=C/c1cc(C)c(Cl)c(C[C@H](C)Cl)c1")

    # now the V3000 version

    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Wrap', 'test_data',  'DoubleBond3000Chiral.mol')
    with open(fileN) as f:
      molBlock = f.read()
    mol =Chem.MolFromMolBlock(molBlock, True, True, False)
    
    smi = Chem.MolToSmiles(mol)
    self.assertEqual(smi,"C=C1C[C@@](C)(C(C)(C)C)CC(=C)/C1=C/c1cc(C)c(Cl)c(C[C@H](C)Cl)c1")
    
    mol = Chem.MolFromMolFile(fileN, True, True,  False)
    smi = Chem.MolToSmiles(mol)
    self.assertEqual(smi,"C=C1C[C@@](C)(C(C)(C)C)CC(=C)/C1=C/c1cc(C)c(Cl)c(C[C@H](C)Cl)c1")

    Chem.SetUseLegacyStereoPerception(True)

  def testValidationMrv(self):
    Chem.SetUseLegacyStereoPerception(False)
    mrvBlock = ""
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Wrap', 'test_data', 'DoubleBond.mrv')
    with open(fileN) as f:
      mrvBlock = f.read()
    mol =Chem.MolFromMrvBlock(mrvBlock, True, True)
    
    smi = Chem.MolToSmiles(mol)
    self.assertEqual(smi,"C=C1CC(C)(C(C)(C)C)CC(=C)C1=Cc1cc(C)c(Cl)c(C[C@H](C)Cl)c1")
    
    mol = Chem.MolFromMrvFile(fileN, True, True)
    smi = Chem.MolToSmiles(mol)
    self.assertEqual(smi,"C=C1CC(C)(C(C)(C)C)CC(=C)C1=Cc1cc(C)c(Cl)c(C[C@H](C)Cl)c1")


    Chem.SetUseLegacyStereoPerception(True)





  def test_picklingWithAddedAttribs(self):
    m = Chem.MolFromSmiles("C")
    m.foo = 1
    m.SetIntProp("bar", 2)
    pkl = pickle.dumps(m)
    nm = pickle.loads(pkl)
    self.assertEqual(nm.GetIntProp("bar"), 2)
    self.assertEqual(nm.foo, 1)

  def testGithubIssue6306(self):
    # test of unpickling
    props = Chem.GetDefaultPickleProperties()
    try:
      Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
      mols = [Chem.MolFromSmiles(s) for s in ["C", "CC"]]
      scaffolds = [MurckoScaffold.GetScaffoldForMol(m) for m in mols]
      # this shouldn't throw an exception
      unpickler = [pickle.loads(pickle.dumps(m)) for m in mols]
    finally:
      Chem.SetDefaultPickleProperties(props)

  @unittest.skipIf(not hasattr(Chem, 'MaeWriter'), "not build with MAEParser support")
  def testMaeWriter(self):
    mol = Chem.MolFromSmiles("C1CCCCC1")
    self.assertTrue(mol)
    title = "random test mol"
    mol.SetProp('_Name', title)

    ofile = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Wrap', 'test_data',
                         'outMaeWriter.mae')
    writer1 = Chem.MaeWriter(ofile)

    osio = StringIO()
    writer2 = Chem.MaeWriter(osio)

    for writer in (writer1, writer2):
      writer.write(mol)
      writer.close()
      del writer

    with open(ofile) as f:
      maefile = f.read()

    self.assertEqual(maefile, osio.getvalue())

    self.assertIn('s_m_m2io_version', maefile)

    self.assertIn('f_m_ct', maefile)

    self.assertIn('s_m_title', maefile)
    self.assertIn(title, maefile)

    self.assertIn(f' m_atom[{mol.GetNumAtoms()}] {{', maefile)

    self.assertEqual(maefile.count("A0A0A0"), 6)  # 6 grey-colored heavy atoms

    self.assertTrue(f' m_bond[{mol.GetNumBonds()}] {{', maefile)

  @unittest.skipIf(not hasattr(Chem, 'MaeWriter'), "not build with MAEParser support")
  def testMaeWriterProps(self):
    mol = Chem.MolFromSmiles("C1CCCCC1")
    self.assertTrue(mol)

    title = "random test mol"
    mol.SetProp('_Name', title)

    boolProp = False
    intProp = 123454321
    realProp = 2.718282  # Mae files have a predefined precision of 6 digits!
    strProp = r"This is a dummy prop, yay!"

    ignored_prop = 'ignored_prop'
    str_dummy_prop = 'str_dummy_prop'
    mol_prop = 'mol_prop'
    atom_prop = 'atom_prop'
    bond_prop = 'bond_prop'
    exported_props = [str_dummy_prop, mol_prop, atom_prop, bond_prop]

    mol.SetIntProp(mol_prop, intProp)
    mol.SetProp(str_dummy_prop, strProp)
    mol.SetProp(ignored_prop, ignored_prop)

    atomIdx = 2
    at = mol.GetAtomWithIdx(atomIdx)
    at.SetDoubleProp(atom_prop, realProp)
    at.SetProp(str_dummy_prop, strProp)
    at.SetProp(ignored_prop, ignored_prop)

    bondIdx = 4
    b = mol.GetBondWithIdx(bondIdx)
    b.SetBoolProp(bond_prop, boolProp)
    b.SetProp(str_dummy_prop, strProp)
    b.SetProp(ignored_prop, ignored_prop)

    heavyAtomColor = "767676"

    osio = StringIO()
    with Chem.MaeWriter(osio) as w:
      w.SetProps(exported_props)
      w.write(mol, heavyAtomColor=heavyAtomColor)

    maestr = osio.getvalue()

    ctBlockStart = maestr.find('f_m_ct')
    atomBlockStart = maestr.find(' m_atom[')
    bondBlockStart = maestr.find(' m_bond[')

    self.assertNotEqual(ctBlockStart, -1)
    self.assertNotEqual(atomBlockStart, -1)
    self.assertNotEqual(bondBlockStart, -1)

    self.assertGreater(bondBlockStart, atomBlockStart)
    self.assertGreater(atomBlockStart, ctBlockStart)

    # structure properties
    self.assertIn(mol_prop, maestr[ctBlockStart:atomBlockStart])
    self.assertIn(str(intProp), maestr[ctBlockStart:atomBlockStart])

    self.assertIn(str_dummy_prop, maestr[ctBlockStart:atomBlockStart])
    self.assertIn(strProp, maestr[ctBlockStart:atomBlockStart])

    self.assertNotIn(ignored_prop, maestr[ctBlockStart:atomBlockStart])

    # atom properties
    self.assertIn(atom_prop, maestr[atomBlockStart:bondBlockStart])
    self.assertIn(str_dummy_prop, maestr[atomBlockStart:bondBlockStart])

    self.assertNotIn(ignored_prop, maestr[atomBlockStart:bondBlockStart])

    for line in maestr[atomBlockStart:bondBlockStart].split('\n'):
      if line.strip().startswith(str(atomIdx + 1)):
        break
    self.assertIn(str(realProp), line)
    self.assertIn(strProp, line)
    self.assertIn(heavyAtomColor, line)

    # bond properties
    self.assertIn(bond_prop, maestr[bondBlockStart:])
    self.assertIn(str_dummy_prop, maestr[bondBlockStart:])

    self.assertNotIn(ignored_prop, maestr[bondBlockStart:])

    for line in maestr[bondBlockStart:].split('\n'):
      if line.strip().startswith(str(bondIdx + 1)):
        break
    self.assertIn(str(int(boolProp)), line)
    self.assertIn(strProp, line)

  @unittest.skipIf(not hasattr(Chem, 'MaeWriter'), "not build with MAEParser support")
  def testMaeWriterRoundtrip(self):
    smiles = "C1CCCCC1"
    mol = Chem.MolFromSmiles(smiles)
    self.assertTrue(mol)

    osio = StringIO()
    with Chem.MaeWriter(osio) as w:
      w.write(mol)

    isio = BytesIO(osio.getvalue().encode())
    with Chem.MaeMolSupplier(isio) as r:
      roundtrip_mol = next(r)
    self.assertTrue(roundtrip_mol)

    self.assertEqual(Chem.MolToSmiles(roundtrip_mol), smiles)

  @unittest.skipIf(not hasattr(Chem, 'MaeWriter'), "not build with MAEParser support")
  def testMaeWriterGetText(self):
    smiles = "C1CCCCC1"
    mol = Chem.MolFromSmiles(smiles)
    self.assertTrue(mol)

    dummy_prop = 'dummy_prop'
    another_dummy_prop = 'another_dummy_prop'
    mol.SetProp(dummy_prop, dummy_prop)
    mol.SetProp(another_dummy_prop, another_dummy_prop)

    heavyAtomColor = "767676"

    osio = StringIO()
    with Chem.MaeWriter(osio) as w:
      w.SetProps([dummy_prop])
      w.write(mol, heavyAtomColor=heavyAtomColor)

    iomae = osio.getvalue()

    ctBlockStart = iomae.find('f_m_ct')
    self.assertNotEqual(ctBlockStart, -1)

    self.assertIn(dummy_prop, iomae)
    self.assertNotIn(another_dummy_prop, iomae)

    mae = Chem.MaeWriter.GetText(mol, heavyAtomColor, -1, [dummy_prop])

    self.assertEqual(mae, iomae[ctBlockStart:])

  def test3dChiralMolFile(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'Cubane.sdf')
    with open(fileN, 'r') as inF:
      inD = inF.read()

    fileWedges = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'CubaneWedges.cxsmi')
    with open(fileWedges, 'r') as inF:
      inWedges = inF.read()

    fileNoWedges = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'FileParsers', 'test_data',
                         'CubaneNoWedges.cxsmi')
    with open(fileNoWedges, 'r') as inF:
      inNoWedges = inF.read()

    m1 = Chem.MolFromMolBlock(inD, sanitize=False, removeHs=False, strictParsing=True)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 16)
    smi = Chem.MolToCXSmiles(m1)
    self.assertTrue(smi == inWedges)

    m1 = Chem.MolFromMolFile(fileN, sanitize=False, removeHs=False, strictParsing=True)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 16)
    smi = Chem.MolToCXSmiles(m1)
    self.assertTrue(smi == inWedges)

    m1 = Chem.MolFromMolBlock(inD, sanitize=False, removeHs=False, strictParsing=True)
    Chem.RemoveNonExplicit3DChirality(m1)

    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 16)
    smi = Chem.MolToCXSmiles(m1)
    self.assertTrue(smi == inNoWedges)

    m1 = Chem.MolFromMolFile(fileN, sanitize=False, removeHs=False, strictParsing=True)
    Chem.RemoveNonExplicit3DChirality(m1)

    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 16)
    smi = Chem.MolToCXSmiles(m1)
    self.assertTrue(smi == inNoWedges)

  def test3dChiralMrvFile(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'MarvinParse', 'test_data',
                         'Cubane.mrv')
    with open(fileN, 'r') as inF:
      inD = inF.read()

    fileWedges = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'MarvinParse', 'test_data',
                         'CubaneWedges.cxsmi')
    with open(fileWedges, 'r') as inF:
      inWedges = inF.read()

    fileNoWedges = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'MarvinParse', 'test_data',
                         'CubaneNoWedges.cxsmi')
    with open(fileNoWedges, 'r') as inF:
      inNoWedges = inF.read()


    m1 = Chem.MolFromMrvBlock(inD, sanitize=False, removeHs=False)

    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 16)
    smi = Chem.MolToCXSmiles(m1)
    sys.stdout.flush()
    self.assertTrue(smi == inWedges)

    m1 = Chem.MolFromMrvFile(fileN, sanitize=False, removeHs=False)

    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 16)
    smi = Chem.MolToCXSmiles(m1)
    self.assertTrue(smi == inWedges)

    m1 = Chem.MolFromMrvBlock(inD, sanitize=False, removeHs=False)
    Chem.RemoveNonExplicit3DChirality(m1)

    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 16)
    smi = Chem.MolToCXSmiles(m1)
    self.assertTrue(smi == inNoWedges)

    m1 = Chem.MolFromMrvFile(fileN, sanitize=False, removeHs=False)
    Chem.RemoveNonExplicit3DChirality(m1)

    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 16)
    smi = Chem.MolToCXSmiles(m1)
    self.assertTrue(smi == inNoWedges)

  def test3dChiralCxsmiles(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'SmilesParse', 'test_data',
                         'Cubane.cxsmi')
    with open(fileN, 'r') as inF:
      inD = inF.read()

    fileWedges = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'SmilesParse', 'test_data',
                         'Cubane.cxsmi.expected3D.cxsmi')
    with open(fileWedges, 'r') as inF:
      inWedges = inF.read()

    fileNoWedges = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'SmilesParse', 'test_data',
                         'Cubane.cxsmi.expected3D2.cxsmi')
    with open(fileNoWedges, 'r') as inF:
      inNoWedges = inF.read()

    ps = Chem.SmilesParserParams()
    ps.allowCXSMILES = True
    ps.parseName = False
    ps.sanitize = False
    ps.removeHs = False

    m1 = Chem.MolFromSmiles(inD, ps)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 16)
    smi = Chem.MolToCXSmiles(m1)
    sys.stdout.flush()
    self.assertTrue(smi == inWedges)


    m1 = Chem.MolFromSmiles(inD, ps)
    Chem.RemoveNonExplicit3DChirality(m1)

    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms() == 16)
    smi = Chem.MolToCXSmiles(m1)
    print('smi: ', smi)
    print('inWedges: ', inNoWedges)

    self.assertTrue(smi == inNoWedges)

  def testReapplyMolBlockWedging(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'MarvinParse', 'test_data',
                         'JDQ443_atrop1.mrv')
  
    fileReapplied = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'MarvinParse', 'test_data',
                         'JDQ443_atrop1.mrv.expected.sdf')
    with open(fileReapplied, 'r') as inF:
      isReapplied = inF.read()

    fileNotReapplied = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'MarvinParse', 'test_data',
                         'JDQ443_atrop1.mrv.expected2.sdf')
    with open(fileNotReapplied, 'r') as inF:
      isNotReapplied = inF.read()

    m = Chem.MolFromMrvFile(fileN, False, False)
    self.assertTrue(m is  not None)
    self.assertTrue(m.GetNumAtoms() == 38)
    mBlock = Chem.MolToMolBlock(m, False, -1, True, True)

    sys.stdout.flush()
    self.assertTrue(mBlock == isNotReapplied)
    Chem.ReapplyMolBlockWedging(m)
    
    mBlock = Chem.MolToMolBlock(m, False, -1, True, True)
    sys.stdout.flush()

    self.assertTrue(mBlock == isReapplied)

  def testReapplyMolBlockWedgingAllBondTypes(self):
    m = Chem.MolFromMolBlock('''
  Mrv2311 04232413302D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 S -11.583 11.3533 0 0
M  V30 2 C -12.9167 10.5833 0 0
M  V30 3 O -11.583 12.8933 0 0
M  V30 4 C -10.2493 10.5833 0 0
M  V30 5 C -10.2493 9.0433 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 4 5
M  V30 2 1 2 1
M  V30 3 1 1 4
M  V30 4 2 1 3 CFG=1
M  V30 END BOND
M  V30 END CTAB
M  END
''')
    self.assertEqual(m.GetBondWithIdx(3).GetBondType(), Chem.BondType.DOUBLE)
    Chem.ReapplyMolBlockWedging(m)
    self.assertEqual(m.GetBondWithIdx(3).GetBondDir(), Chem.BondDir.BEGINWEDGE)
    Chem.ReapplyMolBlockWedging(m, False)
    self.assertEqual(m.GetBondWithIdx(3).GetBondDir(), Chem.BondDir.NONE)

  def testAtropisomerWedging(self):
    m =Chem.MolFromSmiles('CC1=C(N2C=CC=C2[C@H](C)Cl)C(C)CCC1 |(2.679,0.4142,;1.3509,1.181,;0.0229,0.4141,;0.0229,-1.1195,;1.2645,-2.0302,;0.7901,-3.4813,;-0.7446,-3.4813,;-1.219,-2.0302,;-2.679,-1.5609,;-3.0039,-0.0556,;-3.8202,-2.595,;-1.3054,1.1809,;-2.6335,0.4141,;-1.3054,2.7145,;0.0229,3.4813,;1.3509,2.7146,),wD:2.11,wU:8.10,&1:8|')
    self.assertTrue(m is  not None)
    self.assertTrue(m.GetNumAtoms() == 16)
    
    sys.stdout.flush()
    flags = Chem.CXSmilesFields.CX_COORDS | \
                        Chem.CXSmilesFields.CX_MOLFILE_VALUES | \
                        Chem.CXSmilesFields.CX_ATOM_PROPS | \
                        Chem.CXSmilesFields.CX_BOND_CFG | \
                        Chem.CXSmilesFields.CX_ENHANCEDSTEREO
    
    ps = Chem.SmilesWriteParams()
    ps.canonical = True
    smi = Chem.MolToCXSmiles(m, ps, flags, Chem.RestoreBondDirOption.RestoreBondDirOptionTrue)
    self.assertTrue(smi == 'CC1=C(n2cccc2[C@H](C)Cl)C(C)CCC1 |(2.679,0.4142,;1.3509,1.181,;0.0229,0.4141,;0.0229,-1.1195,;1.2645,-2.0302,;0.7901,-3.4813,;-0.7446,-3.4813,;-1.219,-2.0302,;-2.679,-1.5609,;-3.0039,-0.0556,;-3.8202,-2.595,;-1.3054,1.1809,;-2.6335,0.4141,;-1.3054,2.7145,;0.0229,3.4813,;1.3509,2.7146,),wD:2.11,wU:8.10,&1:8|')

    flags = Chem.CXSmilesFields.CX_COORDS | \
                        Chem.CXSmilesFields.CX_MOLFILE_VALUES | \
                        Chem.CXSmilesFields.CX_ATOM_PROPS | \
                        Chem.CXSmilesFields.CX_BOND_ATROPISOMER | \
                        Chem.CXSmilesFields.CX_ENHANCEDSTEREO
    
    
    smi = Chem.MolToCXSmiles(m, ps, flags, Chem.RestoreBondDirOption.RestoreBondDirOptionTrue)

    self.assertTrue(smi == 'CC1=C(n2cccc2[C@H](C)Cl)C(C)CCC1 |(2.679,0.4142,;1.3509,1.181,;0.0229,0.4141,;0.0229,-1.1195,;1.2645,-2.0302,;0.7901,-3.4813,;-0.7446,-3.4813,;-1.219,-2.0302,;-2.679,-1.5609,;-3.0039,-0.0556,;-3.8202,-2.595,;-1.3054,1.1809,;-2.6335,0.4141,;-1.3054,2.7145,;0.0229,3.4813,;1.3509,2.7146,),wD:2.11,&1:8|')

    flags = Chem.CXSmilesFields.CX_COORDS | \
                        Chem.CXSmilesFields.CX_MOLFILE_VALUES | \
                        Chem.CXSmilesFields.CX_ATOM_PROPS | \
                        Chem.CXSmilesFields.CX_ENHANCEDSTEREO
    
    
    smi = Chem.MolToCXSmiles(m, ps, flags, Chem.RestoreBondDirOption.RestoreBondDirOptionTrue)
    self.assertTrue(smi == 'CC1=C(n2cccc2[C@H](C)Cl)C(C)CCC1 |(2.679,0.4142,;1.3509,1.181,;0.0229,0.4141,;0.0229,-1.1195,;1.2645,-2.0302,;0.7901,-3.4813,;-0.7446,-3.4813,;-1.219,-2.0302,;-2.679,-1.5609,;-3.0039,-0.0556,;-3.8202,-2.595,;-1.3054,1.1809,;-2.6335,0.4141,;-1.3054,2.7145,;0.0229,3.4813,;1.3509,2.7146,),&1:8|')

  def test_picklingWithAddedAttribs(self):
    m = Chem.MolFromSmiles("C")
    m.foo = 1
    m.SetIntProp("bar", 2)
    pkl = pickle.dumps(m)
    nm = pickle.loads(pkl)
    self.assertEqual(nm.GetIntProp("bar"), 2)
    self.assertEqual(nm.foo, 1)

  def testGithubIssue6306(self):
    # test of unpickling
    props = Chem.GetDefaultPickleProperties()
    try:
      Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
      mols = [Chem.MolFromSmiles(s) for s in ["C", "CC"]]
      scaffolds = [MurckoScaffold.GetScaffoldForMol(m) for m in mols]
      # this shouldn't throw an exception
      unpickler = [pickle.loads(pickle.dumps(m)) for m in mols]
    finally:
      Chem.SetDefaultPickleProperties(props)

  @unittest.skipIf(not hasattr(Chem, 'MaeWriter'), "not built with MAEParser support")
  def testMaeWriter(self):
    mol = Chem.MolFromSmiles("C1CCCCC1")
    self.assertTrue(mol)
    title = "random test mol"
    mol.SetProp('_Name', title)

    ofile = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Wrap', 'test_data',
                         'outMaeWriter.mae')
    writer1 = Chem.MaeWriter(ofile)

    osio = StringIO()
    writer2 = Chem.MaeWriter(osio)

    for writer in (writer1, writer2):
      writer.write(mol)
      writer.close()
      del writer

    with open(ofile) as f:
      maefile = f.read()

    self.assertEqual(maefile, osio.getvalue())

    self.assertIn('s_m_m2io_version', maefile)

    self.assertIn('f_m_ct', maefile)

    self.assertIn('s_m_title', maefile)
    self.assertIn(title, maefile)

    self.assertIn(f' m_atom[{mol.GetNumAtoms()}] {{', maefile)

    self.assertTrue(f' m_bond[{mol.GetNumBonds()}] {{', maefile)

  @unittest.skipIf(not hasattr(Chem, 'MaeWriter'), "not built with MAEParser support")
  def testMaeWriterProps(self):
    mol = Chem.MolFromSmiles("C1CCCCC1")
    self.assertTrue(mol)

    title = "random test mol"
    mol.SetProp('_Name', title)

    boolProp = False
    intProp = 123454321
    realProp = 2.718282  # Mae files have a predefined precision of 6 digits!
    strProp = r"This is a dummy prop, yay!"

    ignored_prop = 'ignored_prop'
    str_dummy_prop = 'str_dummy_prop'
    mol_prop = 'mol_prop'
    atom_prop = 'atom_prop'
    bond_prop = 'bond_prop'
    exported_props = [str_dummy_prop, mol_prop, atom_prop, bond_prop]

    mol.SetIntProp(mol_prop, intProp)
    mol.SetProp(str_dummy_prop, strProp)
    mol.SetProp(ignored_prop, ignored_prop)

    atomIdx = 2
    at = mol.GetAtomWithIdx(atomIdx)
    at.SetDoubleProp(atom_prop, realProp)
    at.SetProp(str_dummy_prop, strProp)
    at.SetProp(ignored_prop, ignored_prop)

    bondIdx = 4
    b = mol.GetBondWithIdx(bondIdx)
    b.SetBoolProp(bond_prop, boolProp)
    b.SetProp(str_dummy_prop, strProp)
    b.SetProp(ignored_prop, ignored_prop)

    osio = StringIO()
    with Chem.MaeWriter(osio) as w:
      w.SetProps(exported_props)
      w.write(mol)

    maestr = osio.getvalue()

    ctBlockStart = maestr.find('f_m_ct')
    atomBlockStart = maestr.find(' m_atom[')
    bondBlockStart = maestr.find(' m_bond[')

    self.assertNotEqual(ctBlockStart, -1)
    self.assertNotEqual(atomBlockStart, -1)
    self.assertNotEqual(bondBlockStart, -1)

    self.assertGreater(bondBlockStart, atomBlockStart)
    self.assertGreater(atomBlockStart, ctBlockStart)

    # structure properties
    self.assertIn(mol_prop, maestr[ctBlockStart:atomBlockStart])
    self.assertIn(str(intProp), maestr[ctBlockStart:atomBlockStart])

    self.assertIn(str_dummy_prop, maestr[ctBlockStart:atomBlockStart])
    self.assertIn(strProp, maestr[ctBlockStart:atomBlockStart])

    self.assertNotIn(ignored_prop, maestr[ctBlockStart:atomBlockStart])

    # atom properties
    self.assertIn(atom_prop, maestr[atomBlockStart:bondBlockStart])
    self.assertIn(str_dummy_prop, maestr[atomBlockStart:bondBlockStart])

    self.assertNotIn(ignored_prop, maestr[atomBlockStart:bondBlockStart])

    for line in maestr[atomBlockStart:bondBlockStart].split('\n'):
      if line.strip().startswith(str(atomIdx + 1)):
        break
    self.assertIn(str(realProp), line)
    self.assertIn(strProp, line)

    # bond properties
    self.assertIn(bond_prop, maestr[bondBlockStart:])
    self.assertIn(str_dummy_prop, maestr[bondBlockStart:])

    self.assertNotIn(ignored_prop, maestr[bondBlockStart:])

    for line in maestr[bondBlockStart:].split('\n'):
      if line.strip().startswith(str(bondIdx + 1)):
        break
    self.assertIn(str(int(boolProp)), line)
    self.assertIn(strProp, line)

  @unittest.skipIf(not hasattr(Chem, 'MaeWriter'), "not built with MAEParser support")
  def testMaeWriterRoundtrip(self):
    smiles = "C1CCCCC1"
    mol = Chem.MolFromSmiles(smiles)
    self.assertTrue(mol)

    osio = StringIO()
    with Chem.MaeWriter(osio) as w:
      w.write(mol)

    isio = BytesIO(osio.getvalue().encode())
    with Chem.MaeMolSupplier(isio) as r:
      roundtrip_mol = next(r)
    self.assertTrue(roundtrip_mol)

    self.assertEqual(Chem.MolToSmiles(roundtrip_mol), smiles)

  @unittest.skipIf(not hasattr(Chem, 'MaeWriter'), "not built with MAEParser support")
  def testMaeWriterGetText(self):
    smiles = "C1CCCCC1"
    mol = Chem.MolFromSmiles(smiles)
    self.assertTrue(mol)

    dummy_prop = 'dummy_prop'
    another_dummy_prop = 'another_dummy_prop'
    mol.SetProp(dummy_prop, dummy_prop)
    mol.SetProp(another_dummy_prop, another_dummy_prop)

    osio = StringIO()
    with Chem.MaeWriter(osio) as w:
      w.SetProps([dummy_prop])
      w.write(mol)

    iomae = osio.getvalue()

    ctBlockStart = iomae.find('f_m_ct')
    self.assertNotEqual(ctBlockStart, -1)

    self.assertIn(dummy_prop, iomae)
    self.assertNotIn(another_dummy_prop, iomae)

    mae = Chem.MaeWriter.GetText(mol, -1, [dummy_prop])

    self.assertEqual(mae, iomae[ctBlockStart:])

  def test_HapticBondsToDative(self):
    fefile = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'MolStandardize', 'test_data',
                          'ferrocene.mol')
    femol = Chem.MolFromMolFile(fefile)
    newfemol = Chem.rdmolops.HapticBondsToDative(femol)
    self.assertEqual(Chem.MolToSmiles(newfemol),
                     'c12->[Fe+2]3456789(<-c1c->3[cH-]->4c->52)<-c1c->6c->7[cH-]->8c->91')

  def test_DativeBondsToHaptic(self):
    fefile = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'MolStandardize', 'test_data',
                          'ferrocene.mol')
    femol = Chem.MolFromMolFile(fefile)
    newfemol = Chem.rdmolops.HapticBondsToDative(femol)
    backfemol = Chem.rdmolops.DativeBondsToHaptic(newfemol)
    self.assertEqual(Chem.MolToSmiles(femol), Chem.MolToSmiles(backfemol))

  def testTranslateChiralFlag(self):
    mol = Chem.MolFromSmiles("C[C@@](N)(F)C[C@](C)(O)F |a:1|")
    flagMol = Chem.Mol(mol)
    flagMol.SetIntProp("_MolFileChiralFlag", 1)
    Chem.TranslateChiralFlagToStereoGroups(flagMol)
    sgs = flagMol.GetStereoGroups()
    self.assertEqual(len(sgs), 1)
    self.assertEqual(len(sgs[0].GetAtoms()), 2)
    self.assertEqual(sgs[0].GetGroupType(), Chem.StereoGroupType.STEREO_ABSOLUTE)

    flagMol = Chem.Mol(mol)
    flagMol.SetIntProp("_MolFileChiralFlag", 0)
    Chem.TranslateChiralFlagToStereoGroups(flagMol)
    sgs = flagMol.GetStereoGroups()
    self.assertEqual(len(sgs), 2)
    self.assertEqual(sgs[0].GetGroupType(), Chem.StereoGroupType.STEREO_ABSOLUTE)
    self.assertEqual(len(sgs[0].GetAtoms()), 1)

    self.assertEqual(sgs[1].GetGroupType(), Chem.StereoGroupType.STEREO_AND)
    self.assertEqual(len(sgs[1].GetAtoms()), 1)

    flagMol = Chem.Mol(mol)
    flagMol.SetIntProp("_MolFileChiralFlag", 0)
    Chem.TranslateChiralFlagToStereoGroups(flagMol, Chem.StereoGroupType.STEREO_OR)
    sgs = flagMol.GetStereoGroups()
    self.assertEqual(len(sgs), 2)
    self.assertEqual(sgs[0].GetGroupType(), Chem.StereoGroupType.STEREO_ABSOLUTE)
    self.assertEqual(len(sgs[0].GetAtoms()), 1)

    self.assertEqual(sgs[1].GetGroupType(), Chem.StereoGroupType.STEREO_OR)
    self.assertEqual(len(sgs[1].GetAtoms()), 1)

  def testHasQueryHs(self):
    for sma, hasQHs in [
      ("[#1]", (True, False)), ("[#1,N]", (True, True)), ("[$(C-[H])]", (True, False)),
      ("[$([C,#1])]", (True, True)),
      ("[$(c([C;!R;!$(C-[N,O,S]);!$(C-[H])](=O))1naaaa1),$(c([C;!R;!$(C-[N,O,S]);!$(C-[H])](=O))1naa[n,s,o]1)]",
       (True, False))
    ]:
      pat = Chem.MolFromSmarts(sma)
      self.assertEqual(Chem.HasQueryHs(pat), hasQHs)

  def testMolHasQuery(self):
    m1 = Chem.MolFromSmiles("CCO")
    self.assertFalse(m1.HasQuery())

    m2 = Chem.MolFromSmarts("[#6][#6][#8]")
    self.assertTrue(m2.HasQuery())

    m3 = Chem.MolFromSmiles("CC~O")
    self.assertTrue(m3.HasQuery())

  def testMrvHandling(self):
    fn1 = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'MarvinParse', 'test_data',
                       'aspirin.mrv')
    mol = Chem.MolFromMrvFile(fn1)
    self.assertIsNotNone(mol)
    self.assertEqual(mol.GetNumAtoms(), 13)
    mrv = Chem.MolToMrvBlock(mol)
    self.assertTrue('<molecule molID="m1">' in mrv)
    self.assertFalse('<reaction>' in mrv)

    fName = tempfile.NamedTemporaryFile(suffix='.mrv').name
    self.assertFalse(os.path.exists(fName))
    Chem.MolToMrvFile(mol, fName)
    self.assertTrue(os.path.exists(fName))
    os.unlink(fName)

    with open(fn1, 'r') as inf:
      ind = inf.read()
    mol = Chem.MolFromMrvBlock(ind)
    self.assertIsNotNone(mol)
    self.assertEqual(mol.GetNumAtoms(), 13)

  def testAtomMapsInCanonicalization(self):
    mol = Chem.MolFromSmiles("[F:1]C([F:2])O")
    ranks = Chem.CanonicalRankAtoms(mol,breakTies=False,includeAtomMaps=True)
    self.assertNotEqual(ranks[0],ranks[2])
    ranks = Chem.CanonicalRankAtoms(mol,breakTies=False,includeAtomMaps=False)
    self.assertEqual(ranks[0],ranks[2])

  def testExpandAndCollapseAttachmentPoints(self):
    mol = Chem.MolFromSmarts("*CO")
    self.assertEqual(mol.GetNumAtoms(), 3)
    mol.GetAtomWithIdx(2).SetIntProp("molAttchpt", 1)
    Chem.ExpandAttachmentPoints(mol)
    self.assertEqual(mol.GetNumAtoms(), 4)
    self.assertTrue(mol.GetAtomWithIdx(3).HasQuery())

    Chem.CollapseAttachmentPoints(mol)
    self.assertEqual(mol.GetNumAtoms(), 3)

    Chem.ExpandAttachmentPoints(mol,addAsQueries=False)
    self.assertEqual(mol.GetNumAtoms(), 4)
    self.assertFalse(mol.GetAtomWithIdx(3).HasQuery())
    Chem.CollapseAttachmentPoints(mol,markedOnly=False)
    self.assertEqual(mol.GetNumAtoms(), 2)

  def testAddStereoAnnotations(self):
    mol = Chem.MolFromSmiles(
      "C[C@@H]1N[C@H](C)[C@@H]([C@H](C)[C@@H]1C)C1[C@@H](C)O[C@@H](C)[C@@H](C)[C@H]1C/C=C/C |a:5,o1:1,8,o2:14,16,&1:11,18,&2:3,6,r|"
    )
    self.assertIsNotNone(mol)
    Chem.rdCIPLabeler.AssignCIPLabels(mol)
    Chem.AddStereoAnnotations(mol)
    self.assertEqual(mol.GetAtomWithIdx(5).GetProp("atomNote"), "abs (S)")
    self.assertEqual(mol.GetAtomWithIdx(3).GetProp("atomNote"), "and2")

  def testIsRingFused(self):
    molOrig = Chem.MolFromSmiles("C1C(C2CC3CCCCC3C12)C1CCCCC1")
    mol = Chem.RWMol(molOrig)
    ri = mol.GetRingInfo()
    self.assertEqual(ri.NumRings(), 4)
    fusedRings = [ri.IsRingFused(i) for i in range(ri.NumRings())]
    self.assertEqual(fusedRings.count(True), 3)
    self.assertEqual(fusedRings.count(False), 1)
    atoms = mol.GetSubstructMatch(Chem.MolFromSmarts("[$(C1CCC1)]-@[$(C1CCCCC1)]"))
    mol.RemoveBond(*atoms)
    Chem.SanitizeMol(mol)
    self.assertEqual(Chem.MolToSmiles(mol), "C1CCC(CC2CCC2C2CCCCC2)CC1")
    self.assertEqual(ri.NumRings(), 3)
    fusedRings = [ri.IsRingFused(i) for i in range(ri.NumRings())]
    self.assertEqual(fusedRings.count(True), 0)
    self.assertEqual(fusedRings.count(False), 3)
    mol = Chem.RWMol(molOrig)
    ri = mol.GetRingInfo()
    self.assertEqual(ri.NumRings(), 4)
    fusedRings = [ri.IsRingFused(i) for i in range(ri.NumRings())]
    self.assertEqual(fusedRings.count(True), 3)
    self.assertEqual(fusedRings.count(False), 1)
    fusedBonds = [ri.NumFusedBonds(i) for i in range(ri.NumRings())]
    self.assertEqual(fusedBonds.count(0), 1)
    self.assertEqual(fusedBonds.count(1), 2)
    self.assertEqual(fusedBonds.count(2), 1)
    atoms = mol.GetSubstructMatch(Chem.MolFromSmarts("[$(C1CCCCC1-!@[CX4;R1;r4])].[$(C1C(-!@[CX4;R1;r6])CC1)]"))
    mol.AddBond(*atoms, Chem.BondType.SINGLE)
    Chem.SanitizeMol(mol)
    self.assertEqual(Chem.MolToSmiles(mol), "C1CCC2C(C1)CC1C2C2C3CCCCC3C12")
    self.assertEqual(ri.NumRings(), 5)
    fusedRings = [ri.IsRingFused(i) for i in range(ri.NumRings())]
    self.assertEqual(fusedRings.count(True), 5)
    self.assertEqual(fusedRings.count(False), 0)
    fusedBonds = [ri.NumFusedBonds(i) for i in range(ri.NumRings())]
    self.assertEqual(fusedBonds.count(0), 0)
    self.assertEqual(fusedBonds.count(1), 2)
    self.assertEqual(fusedBonds.count(2), 3)


if __name__ == '__main__':
  if "RDTESTCASE" in os.environ:
    suite = unittest.TestSuite()
    testcases = os.environ["RDTESTCASE"]
    for name in testcases.split(':'):
      suite.addTest(TestCase(name))

    runner = unittest.TextTestRunner()
    runner.run(suite)
  else:
    unittest.main()
