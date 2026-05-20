# $Id$
#
#  Copyright (C) 2006  Greg Landrum
#
import os
import pickle
import sys
import unittest

from rdkit import Chem, DataStructs, RDConfig
from rdkit.Chem import MolCatalog


class TestCase(unittest.TestCase):

  def test1(self):
    cat = MolCatalog.CreateMolCatalog()
    es = []
    for smi in ('C1CCC1OC', 'C1CCC1', 'C'):
      m = Chem.MolFromSmiles(smi)
      entry = MolCatalog.MolCatalogEntry()
      entry.SetMol(m)
      self.assertTrue(entry.GetMol())
      eSmi = Chem.MolToSmiles(entry.GetMol())
      self.assertTrue(eSmi == Chem.MolToSmiles(m))
      entry.SetDescription(smi)
      self.assertTrue(entry.GetDescription() == smi)
      es.append(entry)

    v = cat.AddEntry(es[0])
    self.assertTrue(v == 0)
    self.assertTrue(cat.GetNumEntries() == 1)

    v = cat.AddEntry(es[1])
    self.assertTrue(v == 1)
    self.assertTrue(cat.GetNumEntries() == 2)

    v = cat.AddEntry(es[2])
    self.assertTrue(v == 2)
    self.assertTrue(cat.GetNumEntries() == 3)

    cat.AddEdge(0, 1)
    cat.AddEdge(0, 2)
    cat.AddEdge(1, 2)

    d = pickle.dumps(cat)
    es = None
    entry = None
    cat = None

    cat = pickle.loads(d)
    self.assertTrue(cat.GetNumEntries() == 3)
    cat = None

  def test2(self):
    cat = MolCatalog.CreateMolCatalog()
    for smi in ('C1CCC1OC', 'C1CCC1', 'C'):
      m = Chem.MolFromSmiles(smi)
      entry = MolCatalog.MolCatalogEntry()
      entry.SetMol(m)
      entry.SetDescription(smi)
      cat.AddEntry(entry)
    cat.AddEdge(0, 1)
    cat.AddEdge(0, 2)

    pkl = cat.Serialize()
    self.assertIsInstance(pkl, bytes)
    cat2 = MolCatalog.MolCatalog(pkl)
    self.assertEqual(cat2.GetNumEntries(), 3)

  def test3(self):
    m = Chem.MolFromSmiles('C1CCC1')
    entry = MolCatalog.MolCatalogEntry()
    entry.SetMol(m)
    entry.SetDescription('cyclobutane')

    pkl = entry.Serialize()
    self.assertIsInstance(pkl, bytes)
    entry2 = MolCatalog.MolCatalogEntry(pkl)
    self.assertEqual(entry2.GetDescription(), 'cyclobutane')

    d = pickle.dumps(entry)
    entry = None
    entry2 = pickle.loads(d)
    self.assertEqual(entry2.GetDescription(), 'cyclobutane')


if __name__ == '__main__':
  unittest.main()
