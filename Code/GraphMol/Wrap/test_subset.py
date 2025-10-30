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

def check(info, selection):
  for i,v in enumerate(selection):
    if v: assert i in info
    else: assert i not in info

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

    # reuse the existing subsetinfo, it should be properly cleared
    sub = Chem.CopyMolSubset(m, N, info, opts)
    check(info.atomMapping, [True]*7 + [False]*2)
    check(info.bondMapping, [True]*6 + [False]*3)
    for i in N:
      self.assertEqual(info.atomMapping[i], i)
      self.assertEqual(info.bondMapping[i], i)                     
    
    # try the trailing edge
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
