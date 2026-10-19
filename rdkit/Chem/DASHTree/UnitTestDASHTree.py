#
#  Copyright (C) 2026 Marc Lehner and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""Tests for rdkit.Chem.DASHTree, the pure-python side: opening a tree, the
synthetic fixture and the pattern writer. The extension itself is tested in
Code/GraphMol/DASHTree/Wrap/testDASHTree.py.

_writeSyntheticTree() writes test_data/dash_test_tree.dash, a tree that knows
ethanol and nothing else, checked in so the C++ and the extension tests can
share it; the first test asserts the checked-in bytes are what it writes.
"""

import os
import tempfile
import unittest

import numpy as np

from rdkit import Chem, RDConfig
from rdkit.Chem import DASHTree, rdDASHTree
from rdkit.Chem.DASHTree import container

# ethanol with explicit hydrogens is C0 H3 H4 H5, C1 H6 H7, O2 H8. The synthetic
# tree carries exactly its four atom classes, so these are their branch indices;
# in the published trees the same classes are branches 34, 33, 83 and 37.
M, P, O, H = 0, 1, 2, 3
FEATURES = [(6, 4, 0, False, 3), (6, 4, 0, False, 2), (8, 2, 0, False, 1), (1, 1, 0, False, 0)]
STOP = 0x01
NAN = float("nan")

# The synthetic tree, one branch per class. A node is (atomFeature, conAtom,
# conType, firstChild, numChildren, attention, flags, result, std), with
# firstChild relative to the branch root. The values are arbitrary but distinct,
# so a wrong node shows up as a wrong number.
BRANCHES = {
  # methyl C: C0 -> C1 -> O2 -> H8, with no value on the deepest node so the
  # property walk has to fall back to the O node
  M: [(M, -1, -1, 1, 1, 0.0, 0, -0.30, 0.05), (P, 0, 1, 2, 1, 1.0, 0, -0.40, 0.04),
      (O, 1, 1, 3, 1, 1.0, 0, -0.45, 0.03), (H, 2, 1, 0, 0, 1.0, 0, NAN, NAN)],
  # CH2: C1 -> C0 -> O2; the root's second child is never taken
  P: [(P, -1, -1, 1, 2, 0.0, 0, -0.10, 0.10), (M, 0, 1, 3, 1, 1.0, 0, -0.15, 0.10),
      (H, 0, 1, 0, 0, 1.0, 0, -0.99, 0.10), (O, 0, 1, 0, 0, 1.0, 0, -0.18, 0.02)],
  # O: the root's first child is the hydrogen, so a descent that scanned the
  # candidates rather than the children would take C1 instead. Its attention
  # crosses the default threshold, so a default descent stops there and only a
  # caller raising the threshold reaches the node below it.
  O: [(O, -1, -1, 1, 2, 0.0, 0, -0.60, 0.06), (H, 0, 1, 3, 1, 11.0, STOP, -0.65, 0.05),
      (P, 0, 1, 0, 0, 1.0, 0, -0.99, 0.05), (P, 0, 1, 0, 0, 1.0, STOP, -0.70, 0.04)],
  # H: matched through the heavy neighbour, which attaches to nothing
  H: [(H, -1, -1, 1, 3, 0.0, 0, 0.05, 0.02), (M, -1, -1, 4, 1, 1.0, 0, 0.10, 0.02),
      (P, -1, -1, 0, 0, 1.0, 0, 0.08, 0.02), (O, -1, -1, 5, 1, 1.0, 0, 0.30, 0.02),
      (P, 0, 1, 0, 0, 1.0, 0, 0.12, 0.02), (P, 0, 1, 0, 0, 1.0, 0, 0.31, 0.02)],
}


def _writeSyntheticTree(path):
  """Writes BRANCHES as a '.dash' container: result and std float64, and size
  int32, the node id, missing where result is."""
  roots, nodes = [], []
  for b in range(len(FEATURES)):
    roots.append(len(nodes))
    for feat, conAtom, conType, first, nChildren, attn, flags, result, std in BRANCHES[b]:
      nodes.append((container.packKey(feat, conAtom, conType), nChildren, flags,
                    roots[-1] + first if nChildren else 0, attn, result, std))
  records = np.zeros(len(nodes), dtype=container.RECORD_DTYPE)
  for field, column in zip(("key", "nChildren", "flags", "firstChild"), zip(*nodes)):
    records[field] = column
  _, _, _, _, attn, result, std = (np.array(column) for column in zip(*nodes))
  size = np.where(np.isnan(result), container.MISSING_INT32, np.arange(len(nodes)))
  features = np.array([(*f, (0, 0, 0)) for f in FEATURES], dtype=container.FEATURE_DTYPE)
  props = [("result", container.DTYPE_F64, result), ("std", container.DTYPE_F64, std),
           ("size", container.DTYPE_I32, size)]
  container.writeContainer(path, roots, features, records, attn, props, 10.0)


def _fixturePath():
  return os.path.join(RDConfig.RDBaseDir, "Code", "GraphMol", "DASHTree", "test_data",
                      "dash_test_tree.dash")


class TestSyntheticTree(unittest.TestCase):

  def setUp(self):
    self.tree = rdDASHTree.DASHTree(_fixturePath(), ["result", "std"])
    self.mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))

  def testFixtureIsReproducible(self):
    with tempfile.TemporaryDirectory() as d:
      path = os.path.join(d, "tree.dash")
      _writeSyntheticTree(path)
      with open(path, "rb") as f, open(_fixturePath(), "rb") as g:
        self.assertEqual(f.read(), g.read())

  def testPatterns(self):
    self.assertEqual(DASHTree.GetAtomMatchSmarts(self.tree, self.mol, 0),
                     "[#6X4H3+0:1]-[#6X4H2+0]-[#8X2H1+0]")
    path = self.tree.GetAtomNodePath(self.mol, 0)
    self.assertEqual(DASHTree.NodePathToSmarts(self.tree, path, foldHydrogens=False),
                     "[#6X4H3+0:1]-[#6X4H2+0]-[#8X2H1+0]-[#1X1H0+0]")
    self.assertEqual(DASHTree.GetAtomMatchSmarts(self.tree, self.mol, 8),
                     "[#1X1H0+0:1]-[#8X2H1+0]-[#6X4H2+0]")
    self.assertTrue(
      Chem.MolFromSmiles("CCO").HasSubstructMatch(
        Chem.MolFromSmarts(DASHTree.GetAtomMatchSmarts(self.tree, self.mol, 0))))


if __name__ == '__main__':
  unittest.main()
