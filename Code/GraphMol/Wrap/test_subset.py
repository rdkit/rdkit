#
#  Copyright (C) 2025 Hussein Faara, Brian Kelley and other RDKit contributors
#         All Rights Reserved
#
import doctest
import gc
import gzip
import logging
import os
import sys
import tempfile
import unittest
from io import StringIO

from rdkit import Chem

def check(info, selection):
  for i,v in enumerate(selection):
    assert (v and i in info) or (i not in info)

class TestCase(unittest.TestCase):

  def test_subset(self):
    smi = "c1ccccc1CCN"
    m = Chem.MolFromSmiles(smi)
    N = list(range(6))
    sub = Chem.CopyMolSubset(m, N)
    self.assertEqual(Chem.MolToSmiles(sub), "c1ccccc1")
    
    info = Chem.SubsetInfo()
    sub = Chem.CopyMolSubset(m, N, info)
    check(info.atomMapping, [True]*6 + [False]*3)
    check(info.bondMapping, 
                     [True, True, True, True, True, False, False, False, True])
    atoms = [(0,0),(1,1),(2,2),(3,3),(4,4),(5,5)]
    bonds = [(0,0),(1,1),(2,2),(3,3),(4,4),(8,5)]

    for src,dst in atoms:
      self.assertEqual(info.atomMapping[src], dst)
    for src,dst in bonds:
      self.assertEqual(info.bondMapping[src], dst)

    opts = Chem.SubsetOptions()
    opts.method = Chem.SubsetMethod.BONDS;
    sub = Chem.CopyMolSubset(m, N, opts)
    self.assertEqual(Chem.MolToSmiles(sub), "ccccccC")

    sub = Chem.CopyMolSubset(m, N, info, opts)
    check(info.atomMapping, [True]*7 + [False]*2)
    check(info.bondMapping, [True]*6 + [False]*3)
    for i in N:
      self.assertEqual(info.atomMapping[i], i)
      self.assertEqual(info.bondMapping[i], i)                     
    
    N = list(range(6,11))
    sub = Chem.CopyMolSubset(m, N)
    self.assertEqual(Chem.MolToSmiles(sub), "CCN")
    info = Chem.SubsetInfo()
    sub = Chem.CopyMolSubset(m, N, info)
    self.assertEqual(Chem.MolToSmiles(sub), "CCN")
    check(info.atomMapping, [False]*6 + [True]*3)
    check(info.bondMapping, [False]*6+ [True]*2 + [False])

    sub = Chem.CopyMolSubset(m, N, info, opts)
    self.assertEqual(Chem.MolToSmiles(sub), "CCN.cc")
    check(info.atomMapping, [True, False, False, False, False, True, True, True, True])
    check(info.bondMapping, [False, False, False, False, False, False, True, True, True])

    m.SetProp("foo", "bar", computed=True)
    sub = Chem.CopyMolSubset(m, N)
    self.assertTrue(sub.HasProp("foo"))

    opts.clearComputedProps = True
    sub = Chem.CopyMolSubset(m, N, opts)
    self.assertFalse(sub.HasProp("foo"))

  def test_explicit_atoms_bonds_inconsistent_raises(self):
    m = Chem.MolFromSmiles("CCC")  # atoms 0-1-2, bonds: 0 (0-1), 1 (1-2)
    # select atoms 0 and 2, but bond 1 references atom 1 (missing)
    with self.assertRaises(ValueError):
      Chem.CopyMolSubset(m, [0, 2], [1])

  def test_explicit_atoms_bonds_consistent_ok(self):
    m = Chem.MolFromSmiles("CCC")
    sub = Chem.CopyMolSubset(m, [1, 2], [1])
    self.assertEqual(Chem.MolToSmiles(sub), "CC")

  def test_explicit_atoms_out_of_range_raises(self):
    m = Chem.MolFromSmiles("CCC")
    with self.assertRaises(ValueError):
      Chem.CopyMolSubset(m, [0, 3], [0])  # atom index 3 is OOR

  def test_explicit_bonds_out_of_range_raises(self):
    m = Chem.MolFromSmiles("CCC")
    with self.assertRaises(ValueError):
      Chem.CopyMolSubset(m, [0, 1], [2])  # bond index 2 is OOR
