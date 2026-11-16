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
synthetic fixture, the pattern writer and the pruner. The extension itself is
tested in Code/GraphMol/DASHTree/Wrap/testDASHTree.py.

_writeSyntheticTree() writes test_data/dash_test_tree.dash, a tree that knows
ethanol and nothing else, checked in so the C++ and the extension tests can
share it; the first test asserts the checked-in bytes are what it writes.

The tests that need real chemistry run on the pruned tree in Data/DASHTree.
"""

import math
import os
import tempfile
import unittest

import numpy as np

from rdkit import Chem, RDConfig
from rdkit.Chem import DASHTree, rdDASHTree
from rdkit.Chem.DASHTree import container, convert, prune

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

  def testPrunedTreeFallsBack(self):
    # at 0.1 the O node below the methyl's CH2 node is cut, so atom 0's descent
    # stops at the CH2 node and reads its value
    with tempfile.TemporaryDirectory() as d:
      out = os.path.join(d, "pruned.dash")
      prune.main([_fixturePath(), out, "--tolerance", "0.1"])
      pruned = rdDASHTree.DASHTree(out, ["result"])
      self.assertEqual(list(pruned.GetAtomNodePath(self.mol, 0)), [M, 0, 1])
      self.assertEqual(pruned.GetAtomProperty(self.mol, 0, "result"), -0.40)
      # the CH2's methyl child is cut but stays as a placeholder ahead of the
      # kept H, so atom 1 stops there and reads the root
      self.assertEqual(pruned.GetAtomProperty(self.mol, 1, "result"), -0.10)
      del pruned  # unmaps the file before the directory goes

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


class TestBundledTree(unittest.TestCase):
  """The pruned tree in Data/DASHTree, which GetDASHTree() falls back to."""

  @classmethod
  def setUpClass(cls):
    cls.saved = os.environ.pop("RDKIT_DASH_TREE", None)
    cls.tree = DASHTree.GetDASHTree()

  @classmethod
  def tearDownClass(cls):
    if cls.saved is not None:
      os.environ["RDKIT_DASH_TREE"] = cls.saved

  def testIsThePublishedTreePruned(self):
    self.assertEqual(self.tree.GetFilename(), DASHTree.bundledTreePath)
    self.assertEqual(self.tree.GetNumBranches(), 122)
    expected = [(z, d, q, bool(c), h) for z, d, q, c, h in convert.STANDARD_FEATURES]
    self.assertEqual([self.tree.GetAtomFeature(i) for i in range(122)], expected)
    self.assertLess(self.tree.GetNumNodes(), 2500000)  # the full tree has 9.4 million
    self.assertEqual(sorted(self.tree.GetAvailablePropertyNames()), ["result", "std"])

  def testEnvironmentOverrides(self):
    os.environ["RDKIT_DASH_TREE"] = _fixturePath()
    try:
      self.assertEqual(DASHTree.GetDASHTree().GetNumBranches(), 4)
    finally:
      del os.environ["RDKIT_DASH_TREE"]

  def testChargesStayCloseToThePythonPackage(self):
    # test_data/dash_ref_values.txt holds what the DASH-tree python package
    # computes with the full tree; pruning at 0.02 e was measured to move the
    # std-weighted charges by 0.0065 e RMS and 0.044 e at most. The package
    # takes the first of several atoms that fit a node, hence the legacy rule.
    testData = os.path.join(RDConfig.RDBaseDir, "Code", "GraphMol", "DASHTree", "test_data")
    mols = [
      m for m in Chem.SDMolSupplier(os.path.join(testData, "dash_ref_mols.sdf"), removeHs=False)
      if m is not None
    ]
    self.assertEqual(len(mols), 100)
    options = rdDASHTree.NormalizationOptions()
    options.params.legacyTieBreaking = True
    charges = [self.tree.GetPartialCharges(m, options) for m in mols]
    diffs = []
    with open(os.path.join(testData, "dash_ref_values.txt")) as f:
      for line in f:
        parts = line.split()
        if not parts or parts[0] == "#":
          continue
        molIdx, atomIdx, reference = int(parts[0]), int(parts[1]), float(parts[7])
        diffs.append(charges[molIdx][atomIdx] - reference)
    self.assertEqual(len(diffs), sum(m.GetNumAtoms() for m in mols))
    self.assertLess(max(abs(d) for d in diffs), 0.06)
    self.assertLess(math.sqrt(sum(d * d for d in diffs) / len(diffs)), 0.01)
    for mol, q in zip(mols, charges):
      self.assertAlmostEqual(sum(q), sum(a.GetFormalCharge() for a in mol.GetAtoms()), places=10)

  def testPatternRefindsItsOwnAtom(self):
    mol = Chem.AddHs(Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1"))
    for atom in mol.GetAtoms():
      smarts = DASHTree.GetAtomMatchSmarts(self.tree, mol, atom.GetIdx())
      query = Chem.MolFromSmarts(smarts)
      mapped = next(a.GetIdx() for a in query.GetAtoms() if a.GetAtomMapNum() == 1)
      landings = {m[mapped] for m in mol.GetSubstructMatches(query, uniquify=False)}
      self.assertIn(atom.GetIdx(), landings, smarts)


if __name__ == '__main__':
  unittest.main()
