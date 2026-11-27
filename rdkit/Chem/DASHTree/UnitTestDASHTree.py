#
#  Copyright (C) 2026 Marc Lehner and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""Tests for rdkit.Chem.DASHTree, the pure-python side: GetDASHTree(), the
synthetic fixture, the pattern writer, the pruner and the populator.
Code/GraphMol/DASHTree/Wrap/testDASHTree.py tests the extension itself.

_writeSyntheticTree() writes test_data/dash_test_tree.dash, a tree that knows
ethanol and nothing else. The file is in the repository, so the C++ and the
extension tests can share it. testFixtureIsReproducible asserts that
_writeSyntheticTree() writes the same bytes as the copy in the repository.

The tests that need real chemistry run on initial_test, the pruned tree the
build downloads into Data/DASHTree.
"""

import contextlib
import io
import math
import os
import tempfile
import unittest
from unittest import mock

import numpy as np

from rdkit import Chem, RDConfig
from rdkit.Chem import DASHTree, rdDASHTree
from rdkit.Chem.DASHTree import container, convert, populate, prune

# ethanol with explicit hydrogens is C0 H3 H4 H5, C1 H6 H7, O2 H8. The synthetic
# tree carries exactly its four atom classes, so these are their branch indices.
# In the published trees the same classes are branches 34, 33, 83 and 37.
M, P, O, H = 0, 1, 2, 3
FEATURES = [(6, 4, 0, False, 3), (6, 4, 0, False, 2), (8, 2, 0, False, 1), (1, 1, 0, False, 0)]
STOP = container.NODE_FLAG_STOP
NAN = float("nan")

# The synthetic tree has one branch per class. A node is (atomFeature, conAtom,
# conType, firstChild, numChildren, attention, flags, result, std), with
# firstChild relative to the branch root. The results are arbitrary but distinct
# on the nodes a descent can take, so a wrong node gives a wrong number.
BRANCHES = {
  # methyl C: C0 -> C1 -> O2 -> H8, with no value on the deepest node, so the
  # property walk reads the O node instead
  M: [(M, -1, -1, 1, 1, 0.0, 0, -0.30, 0.05), (P, 0, 1, 2, 1, 1.0, 0, -0.40, 0.04),
      (O, 1, 1, 3, 1, 1.0, 0, -0.45, 0.03), (H, 2, 1, 0, 0, 1.0, 0, NAN, NAN)],
  # CH2: C1 -> C0 -> O2, and no descent takes the root's second child
  P: [(P, -1, -1, 1, 2, 0.0, 0, -0.10, 0.10), (M, 0, 1, 3, 1, 1.0, 0, -0.15, 0.10),
      (H, 0, 1, 0, 0, 1.0, 0, -0.99, 0.10), (O, 0, 1, 0, 0, 1.0, 0, -0.18, 0.02)],
  # O: the root's first child is the hydrogen, so a descent that scanned the
  # candidates rather than the children would take C1 instead. Its attention
  # crosses the default threshold, so a default descent stops there and only a
  # caller that raises the threshold reaches the node below it.
  O: [(O, -1, -1, 1, 2, 0.0, 0, -0.60, 0.06), (H, 0, 1, 3, 1, 11.0, STOP, -0.65, 0.05),
      (P, 0, 1, 0, 0, 1.0, 0, -0.99, 0.05), (P, 0, 1, 0, 0, 1.0, STOP, -0.70, 0.04)],
  # H: matched through the heavy neighbour, which attaches to nothing
  H: [(H, -1, -1, 1, 3, 0.0, 0, 0.05, 0.02), (M, -1, -1, 4, 1, 1.0, 0, 0.10, 0.02),
      (P, -1, -1, 0, 0, 1.0, 0, 0.08, 0.02), (O, -1, -1, 5, 1, 1.0, 0, 0.30, 0.02),
      (P, 0, 1, 0, 0, 1.0, 0, 0.12, 0.02), (P, 0, 1, 0, 0, 1.0, 0, 0.31, 0.02)],
}


def _writeSyntheticTree(path):
  """Writes BRANCHES as a '.dash' container: result and std as float64, and size
  as int32. size is the node id, or missing where result is NaN."""
  roots, nodes = [], []
  for b in range(len(FEATURES)):
    roots.append(len(nodes))
    for feat, conAtom, conType, first, nChildren, attn, flags, result, std in BRANCHES[b]:
      nodes.append((container.packKey(feat, conAtom, conType), nChildren, flags,
                    roots[-1] + first if nChildren else 0, attn, result, std))
  records = np.array([n[:4] for n in nodes], dtype=container.RECORD_DTYPE)
  attn, result, std = np.array([n[4:] for n in nodes]).T
  size = np.where(np.isnan(result), container.MISSING_INT32, np.arange(len(nodes)))
  features = np.array([(*f, (0, 0, 0)) for f in FEATURES], dtype=container.FEATURE_DTYPE)
  props = [("result", container.DTYPE_F64, result), ("std", container.DTYPE_F64, std),
           ("size", container.DTYPE_I32, size)]
  container.writeContainer(path, roots, features, records, attn, props,
                           container.DEFAULT_ATTENTION_THRESHOLD)


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
    # at 0.1, pruning cuts the O node below the methyl's CH2 node, so atom 0's
    # descent stops at the CH2 node and reads its value
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as d:
      out = os.path.join(d, "pruned.dash")
      prune.main([_fixturePath(), out, "--tolerance", "0.1"])
      pruned = rdDASHTree.DASHTree(out, ["result"])
      self.assertEqual(list(pruned.GetAtomNodePath(self.mol, 0)), [M, 0, 1])
      self.assertEqual(pruned.GetAtomProperty(self.mol, 0, "result"), -0.40)
      # pruning cuts the CH2's methyl child, which stays as a placeholder ahead
      # of the kept H, so atom 1 stops there and reads the root
      self.assertEqual(pruned.GetAtomProperty(self.mol, 1, "result"), -0.10)
      # the hydrogen rule keeps the H root's children, so atom 3 reads its methyl node
      self.assertEqual(pruned.GetAtomProperty(self.mol, 3, "result"), 0.10)
      del pruned  # unmaps the file before the directory goes

      # at 0.01, only the methyl branch's deepest node goes. It has no value,
      # so the descent reads the O node above it, as before
      prune.main([_fixturePath(), out, "--tolerance", "0.01"])
      pruned = rdDASHTree.DASHTree(out, ["result"])
      self.assertEqual(pruned.GetNumNodes(), self.tree.GetNumNodes() - 1)
      self.assertEqual(list(pruned.GetAtomNodePath(self.mol, 0)), [M, 0, 1, 2])
      self.assertEqual(pruned.GetAtomProperty(self.mol, 0, "result"), -0.45)
      del pruned

      # the same holds for an int32 column, where the leaf has no size
      prune.main([_fixturePath(), out, "--tolerance", "0.5", "--value-property", "size"])
      pruned = rdDASHTree.DASHTree(out, ["result"])
      self.assertEqual(pruned.GetNumNodes(), self.tree.GetNumNodes() - 1)
      del pruned

  def testPopulate(self):
    # two conformers: the values are the atom index, then the index + 2. Atom 8,
    # the OH hydrogen, has none, so the nodes only it reaches stay empty
    values = np.array([np.arange(9.0), np.arange(9.0) + 2])
    values[:, 8] = NAN
    columns, failed = populate.collectStatistics(_fixturePath(), [(self.mol, {"q": values})], ["q"])
    self.assertEqual(failed, 0)
    self.assertEqual([(name, tag) for name, tag, _ in columns], [("q", container.DTYPE_F32),
                                                                 ("q_std", container.DTYPE_F16),
                                                                 ("q_count", container.DTYPE_I32)])
    self.assertEqual(columns[2][2][15], container.MISSING_INT32)  # only atom 8 reaches node 15
    # the sums are the same when they go in after each record
    with mock.patch.object(populate, "_BATCH", 0):
      batched, _ = populate.collectStatistics(_fixturePath(), [(self.mol, {
        "q": values[:1]
      }), (self.mol, {
        "q": values[1]
      })], ["q"])
    for (_, _, a), (_, _, b) in zip(columns, batched):
      np.testing.assert_array_equal(a, b)

    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as d:
      out = os.path.join(d, "populated.dash")
      populate.writeTree(_fixturePath(), out, columns)
      tree = rdDASHTree.DASHTree(out)
      self.assertEqual(sorted(tree.GetAvailablePropertyNames()),
                       ["q", "q_count", "q_std", "result", "size", "std"])

      def read(i, p):
        return tree.GetAtomProperty(self.mol, i, p)

      # atom 0 is the only atom on its deepest node, with the values 0 and 2
      self.assertEqual((read(0, "q"), read(0, "q_count")), (1.0, 2))
      self.assertAlmostEqual(read(0, "q_std"), math.sqrt(2), places=3)
      # the three methyl Hs share their deepest node: 3, 4, 5, 5, 6, 7
      self.assertEqual((read(3, "q"), read(3, "q_count")), (5.0, 6))
      # atom 8 reads the H root, which the five other Hs reach
      self.assertEqual((read(8, "q"), read(8, "q_count")), (6.0, 10))
      del tree, read

      # pruning on q drops the subtrees with the parent's mean and the subtrees
      # without a q value. The H root keeps its children, so 7 nodes stay
      pruned = os.path.join(d, "pruned.dash")
      prune.main([out, pruned, "--tolerance", "0.001", "--value-property", "q"])
      tree = rdDASHTree.DASHTree(pruned)
      self.assertEqual(tree.GetNumNodes(), 7)
      self.assertEqual(list(tree.GetAtomNodePath(self.mol, 8)), [H, 0, 3])
      self.assertEqual(tree.GetAtomProperty(self.mol, 8, "q"), 6.0)
      del tree

      # a new column replaces the column of the same name
      populate.writeTree(_fixturePath(), out, [("result", ) + columns[0][1:]])
      tree = rdDASHTree.DASHTree(out)
      self.assertEqual(sorted(tree.GetAvailablePropertyNames()), ["result", "size", "std"])
      self.assertEqual(tree.GetAtomProperty(self.mol, 3, "result"), 5.0)  # the fixture has 0.12
      del tree
      with self.assertRaises(ValueError):
        populate.writeTree(out, out, columns)

    with self.assertRaises(ValueError):  # two values for nine atoms
      populate.collectStatistics(_fixturePath(), [(self.mol, {"q": [1.0, 2.0]})], ["q"])
    with self.assertRaises(ValueError):  # values for a property that is not collected
      populate.collectStatistics(_fixturePath(), [(self.mol, {"r": values})], ["q"])
    with self.assertRaises(ValueError):
      populate.collectStatistics(_fixturePath(), [(self.mol, {"q": np.full(9, np.inf)})], ["q"])
    with self.assertRaises(ValueError):  # implicit hydrogens
      populate.collectStatistics(_fixturePath(), [(Chem.MolFromSmiles("CCO"), {})], ["q"])
    with self.assertRaises(ValueError):  # 'q_std' is also a column of 'q'
      populate.collectStatistics(_fixturePath(), [], ["q", "q_std"])
    # the tree has no class for the N of methylamine, so the molecule adds
    # nothing. A molecule without atoms adds nothing, and does not fail
    amine = Chem.AddHs(Chem.MolFromSmiles("CN"))
    records = [(amine, {"q": np.ones(amine.GetNumAtoms())}), (Chem.Mol(), {})]
    columns, failed = populate.collectStatistics(_fixturePath(), records, ["q"])
    self.assertEqual(failed, 1)
    self.assertTrue((columns[2][2] == container.MISSING_INT32).all())

  def testPopulateFromSdf(self):
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as d:
      sdf = os.path.join(d, "mols.sdf")
      with Chem.SDWriter(sdf) as writer:
        for shift in (0, 2):  # the two conformers of testPopulate, one record each
          mol = Chem.Mol(self.mol)
          for atom in mol.GetAtoms():
            if atom.GetIdx() != 8:
              atom.SetDoubleProp("q", atom.GetIdx() + shift)
          Chem.CreateAtomDoublePropertyList(mol, "q")
          writer.write(mol)
        writer.write(Chem.AddHs(Chem.MolFromSmiles("CN")))  # the tree has no class for the N
      with open(sdf, "a") as f:  # a record that does not parse
        f.write("broken\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n$$$$\n")
      out = os.path.join(d, "populated.dash")
      err = io.StringIO()
      with contextlib.redirect_stderr(err):
        populate.main([_fixturePath(), out, sdf, "--props", "q", "--drop-existing"])
      self.assertIn("skipped 1 records that do not parse and 1 molecules", err.getvalue())
      tree = rdDASHTree.DASHTree(out)
      self.assertEqual(sorted(tree.GetAvailablePropertyNames()), ["q", "q_count", "q_std"])
      self.assertEqual([tree.GetAtomProperty(self.mol, i, "q") for i in (0, 3, 8)], [1.0, 5.0, 6.0])
      del tree

      other = os.path.join(d, "other.dash")
      for args in ([out, out, sdf, "--props",
                    "q"], [out, other, sdf, "--props",
                           "q,q_std"], [out, other, sdf, "--props",
                                        ","], [out, other, other, "--props", "q"]):
        with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()):
          populate.main(args)
      self.assertFalse(os.path.exists(other))

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
    with self.assertRaises(ValueError):  # the path of a hydrogen that no node describes
      DASHTree.NodePathToSmarts(self.tree, [H])


class TestBundledTree(unittest.TestCase):
  """initial_test in Data/DASHTree, which GetDASHTree() opens by default."""

  @classmethod
  def setUpClass(cls):
    cls.path = os.path.join(RDConfig.RDDataDir, "DASHTree", "initial_test.dash")
    assert os.path.isfile(cls.path), "the build downloads initial_test into Data/DASHTree"
    cls.saved = os.environ.pop("RDKIT_DASH_TREE", None)
    cls.tree = DASHTree.GetDASHTree()

  @classmethod
  def tearDownClass(cls):
    if cls.saved is not None:
      os.environ["RDKIT_DASH_TREE"] = cls.saved

  def testIsThePublishedTreePruned(self):
    self.assertEqual(self.tree.GetFilename(), self.path)
    self.assertEqual(DASHTree.GetDASHTree("initial_test").GetFilename(), self.path)
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
    # computes with the full tree. Pruning at 0.02 e was measured to move the
    # std-weighted charges by 0.0065 e RMS and 0.044 e at most. The package
    # takes the first of several atoms that fit a node, so the test sets
    # legacyTieBreaking.
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

  def testPopulateFollowsTheTieBreaking(self):
    # the two CH2 of 5-aminoisoindoline fit the same node, so the two rules
    # descend to different nodes and count the values there
    mol = Chem.AddHs(Chem.MolFromSmiles("Nc1ccc2CNCc2c1"))
    legacy = rdDASHTree.DASHParams()
    legacy.legacyTieBreaking = True
    records = [(mol, {"q": np.ones(mol.GetNumAtoms())})]
    counts = [
      populate.collectStatistics(self.path, records, ["q"], params)[0][2][2]
      for params in (None, legacy)
    ]
    self.assertFalse(np.array_equal(*counts))

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
