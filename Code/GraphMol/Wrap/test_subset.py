#
#  Copyright (C) 2003-2021  Greg Landrum and other RDKit contributors
#         All Rights Reserved
#
""" This is a rough coverage test of the python wrapper

it's intended to be shallow, but broad

"""

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


class TestCase(unittest.TestCase):

  def test_subset(self):
    smi = "c1ccccc1CCN"
    m = Chem.MolFromSmiles(smi)
    N = list(range(6))
    sub = Chem.CopyMolSubset(m, N)
    self.assertEqual(Chem.MolToSmiles(sub), "c1ccccc1")
    
    info = Chem.SubsetInfo()
    sub = Chem.CopyMolSubset(m, N, info)
    self.assertEqual(list(info.selectedAtoms), [True]*6 + [False]*3)
    self.assertEqual(list(info.selectedBonds),
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
    
    info = Chem.SubsetInfo()
    sub = Chem.CopyMolSubset(m, N, info, opts)
    self.assertEqual(list(info.selectedAtoms), [True]*7 + [False]*2)
    self.assertEqual(list(info.selectedBonds), [True]*6 + [False]*3)
    for i in N:
      self.assertEqual(info.atomMapping[i], i)
      self.assertEqual(info.bondMapping[i], i)                     
    
    # try the trailing edge
    N = list(range(6,11))
    sub = Chem.CopyMolSubset(m, N)
    self.assertEqual(Chem.MolToSmiles(sub), "CCN")
    sub = Chem.CopyMolSubset(m, N, info)
    self.assertEqual(Chem.MolToSmiles(sub), "CCN")
    self.assertEqual(list(info.selectedAtoms), [False]*6 + [True]*3)
    self.assertEqual(list(info.selectedBonds), [False]*6+ [True]*2 + [False])

    sub = Chem.CopyMolSubset(m, N, info, opts)
    self.assertEqual(Chem.MolToSmiles(sub), "CCN.cc")
    self.assertEqual(list(info.selectedAtoms), [True, False, False, False, False, True, True, True, True])
    self.assertEqual(list(info.selectedBonds), [False, False, False, False, False, False, True, True, True])

    m.SetProp("foo", "bar", computed=True)
    sub = Chem.CopyMolSubset(m, N)
    self.assertTrue(sub.HasProp("foo"))

    opts.clearComputedProps = True
    sub = Chem.CopyMolSubset(m, N, opts)
    self.assertFalse(sub.HasProp("foo"))
    
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
