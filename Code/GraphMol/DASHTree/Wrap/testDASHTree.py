#
#  Copyright (C) 2026 Marc Lehner and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""Tests for rdkit.Chem.rdDASHTree.

Most run against test_data/dash_test_tree.dash, a synthetic tree written by
rdkit/Chem/DASHTree/UnitTestDASHTree.py that knows ethanol and nothing else:
its values are made up, so what it exercises is the machinery -- the descent,
the property fallback, the stop flag, the normalisations, the node API. The
tests that need real chemistry use the pruned tree in Data/DASHTree.
TestReferenceAgreement compares against the full published tree, which is not
in the repository, and runs only if DASH_TREE_FILE names it.
"""

import math
import os
import unittest

from rdkit import Chem, RDConfig
from rdkit.Chem import rdDASHTree

TEST_DATA_DIR = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'DASHTree', 'test_data')
FIXTURE = os.path.join(TEST_DATA_DIR, 'dash_test_tree.dash')

# ethanol with explicit hydrogens is C0 H3 H4 H5, C1 H6 H7, O2 H8; the synthetic
# tree carries exactly its four atom classes, in this order
M, P, O, H = 0, 1, 2, 3
FEATURES = [(6, 4, 0, False, 3), (6, 4, 0, False, 2), (8, 2, 0, False, 1), (1, 1, 0, False, 0)]

SHIPPED_TREE = os.path.join(RDConfig.RDDataDir, 'DASHTree', 'default_pruned.dash')
PUBLISHED_TREE = os.environ.get('DASH_TREE_FILE')


def molWithHs(smiles):
  mol = Chem.MolFromSmiles(smiles)
  assert mol is not None, smiles
  return Chem.AddHs(mol)


needsPublishedTree = unittest.skipUnless(
  PUBLISHED_TREE, 'set DASH_TREE_FILE to the full published tree, converted by '
  'python -m rdkit.Chem.DASHTree.convert')


class TestSyntheticTree(unittest.TestCase):
  """The machinery, on the fixture: every value it holds is known."""

  def setUp(self):
    self.tree = rdDASHTree.DASHTree(FIXTURE, ["result", "std"])
    self.mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))

  def testFeatureIndices(self):
    # the classes come from the table in the file, not from the code
    got = [self.tree.GetAtomFeatureIndex(a) for a in self.mol.GetAtoms()]
    self.assertEqual(got, [M, P, O, H, H, H, H, H, H])
    self.assertEqual([self.tree.GetAtomFeature(i) for i in range(4)], FEATURES)
    benzene = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    self.assertEqual(self.tree.GetAtomFeatureIndex(benzene.GetAtomWithIdx(0)), -1)
    with self.assertRaises(IndexError):
      self.tree.GetAtomFeature(4)

  def testDescent(self):
    self.assertEqual(self.tree.GetNumBranches(), 4)
    self.assertEqual(self.tree.GetNumNodes(), 18)
    expected = {
      0: [M, 0, 1, 2, 3],
      1: [P, 0, 1, 3],
      2: [O, 0, 1],
      3: [H, 0, 1, 4],
      6: [H, 0, 2],
      8: [H, 0, 3, 5]
    }
    for atom, path in expected.items():
      self.assertEqual(self.tree.GetAtomNodePath(self.mol, atom), path)
    self.assertEqual(self.tree.GetMatchedSubstructure(self.mol, 0), (0, 1, 2, 8))
    self.assertEqual(self.tree.GetMatchedSubstructure(self.mol, 3), (0, 1))
    # the O descent stops at the flagged node; told to read the attentions
    # instead of the flag, it carries on to the node below
    params = rdDASHTree.DASHParams()
    params.attentionThreshold = 100.0
    self.assertEqual(self.tree.GetAtomNodePath(self.mol, 2, params), [O, 0, 1, 3])
    self.assertAlmostEqual(self.tree.GetAtomProperty(self.mol, 2, "result"), -0.65)
    self.assertAlmostEqual(self.tree.GetAtomProperty(self.mol, 2, "result", params), -0.70)
    params.attentionIncrementThreshold = 2.0  # every step on branch M adds 1.0
    self.assertEqual(self.tree.GetAtomNodePath(self.mol, 0, params), [M, 0, 1])

  def testCharges(self):
    raw = [-0.45, -0.18, -0.65, 0.12, 0.12, 0.12, 0.08, 0.08, 0.31]  # deepest node with a value
    std = [0.03, 0.02, 0.05, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02]
    deficit = 0.0 - sum(raw)
    details = self.tree.GetPartialChargesDetails(self.mol)
    for got, want in zip(details["raw"], raw):
      self.assertAlmostEqual(got, want)
    for got, want in zip(details["std"], std):
      self.assertAlmostEqual(got, want)
    self.assertEqual(details["matchDepth"], [4, 3, 2, 3, 3, 3, 2, 2, 3])
    N = rdDASHTree.Normalization
    for normalization, want in ((N.NONE, raw), (N.SYMMETRIC, [r + deficit / 9 for r in raw]),
                                (N.STD_WEIGHTED,
                                 [r + deficit * s / sum(std) for r, s in zip(raw, std)])):
      options = rdDASHTree.NormalizationOptions()
      options.normalization = normalization
      for got, w in zip(self.tree.GetPartialCharges(self.mol, options), want):
        self.assertAlmostEqual(got, w)
    # a batch agrees with the single-molecule call, and None yields an empty row
    batch = self.tree.GetPartialChargesBatch([self.mol, None, self.mol], numThreads=0)
    self.assertEqual(batch[1], [])
    self.assertEqual(batch[0], self.tree.GetPartialCharges(self.mol))
    self.assertEqual(batch[0], batch[2])
    # molecules that nothing but the generator holds stay alive for the batch
    self.assertEqual(self.tree.GetPartialChargesBatch(Chem.Mol(self.mol) for _ in range(2)),
                     [batch[0], batch[0]])

  def testNormalizedMolProperty(self):
    # charges are the normalised assignment with the formal charge as target
    charges = self.tree.GetPartialCharges(self.mol)
    self.assertEqual(self.tree.GetNormalizedMolProperty(self.mol, 0.0), charges)
    self.assertAlmostEqual(sum(self.tree.GetNormalizedMolProperty(self.mol, 1.0)), 1.0, places=12)
    details = self.tree.GetPartialChargesDetails(self.mol)
    batch = self.tree.GetNormalizedMolPropertyBatch([self.mol, None, self.mol], [0.0, 0.0, 1.0])
    self.assertEqual(batch[0], charges)
    self.assertEqual(batch[1], [])
    self.assertAlmostEqual(sum(batch[2]), 1.0, places=12)
    with self.assertRaises(ValueError):  # one target per molecule
      self.tree.GetNormalizedMolPropertyBatch([self.mol], [0.0, 1.0])
    options = rdDASHTree.NormalizationOptions()
    options.normalization = rdDASHTree.Normalization.NONE
    options.valueProperty = "std"
    self.assertEqual(self.tree.GetNormalizedMolProperty(self.mol, 1.0, options), details["std"])
    # a column without a deviation column still normalises, except by weight
    options.stdProperty = "no_such_column"
    options.normalization = rdDASHTree.Normalization.SYMMETRIC
    self.assertAlmostEqual(sum(self.tree.GetNormalizedMolProperty(self.mol, 1.0, options)), 1.0,
                           places=12)
    options.normalization = rdDASHTree.Normalization.STD_WEIGHTED
    with self.assertRaises(ValueError):
      self.tree.GetNormalizedMolProperty(self.mol, 1.0, options)
    options.stdProperty = "std"
    options.defaultStdValue = 0.0
    with self.assertRaises(ValueError):  # would divide by the deviation total
      self.tree.GetNormalizedMolProperty(self.mol, 1.0, options)

  def testRefusals(self):
    with self.assertRaises(ValueError):  # methane's carbon is not a DASH class
      self.tree.GetPartialCharges(Chem.AddHs(Chem.MolFromSmiles("C")))
    with self.assertRaises(ValueError):
      rdDASHTree.DASHTree(FIXTURE, ["result", "definitely_absent"])
    with self.assertRaises(OSError):
      rdDASHTree.DASHTree('no_such_file_at_all.dash')
    with self.assertRaises(AttributeError):  # a typo is not silently stored
      rdDASHTree.DASHParams().maxdepth = 3

  def testNodes(self):
    root = self.tree.GetRoot(O)
    self.assertEqual((root.GetBranch(), root.GetId(), root.GetAtomFeatureIndex()), (O, 0, O))
    self.assertEqual((root.GetConAtom(), root.GetConType()), (-1, -1))
    self.assertEqual(len(root), 2)
    first, second = root
    self.assertEqual((first.GetAtomFeatureIndex(), first.GetConAtom(), first.GetConType()),
                     (H, 0, 1))
    self.assertAlmostEqual(first.GetAttention(), 11.0)
    self.assertTrue(first.IsStopNode())
    self.assertFalse(second.IsStopNode())
    self.assertAlmostEqual(second.GetValue("result"), -0.99)
    self.assertEqual(root[-1].GetId(), second.GetId())
    self.assertEqual(first[0].GetId(), 3)
    self.assertTrue(math.isnan(self.tree.GetNode(M, 3).GetValue("result")))
    with self.assertRaises(IndexError):
      root[2]
    with self.assertRaises(IndexError):
      root.GetChild(2)
    with self.assertRaises(IndexError):
      self.tree.GetNode(4, 0)
    with self.assertRaises(ValueError):
      root.GetValue("no_such_property")
    # a node keeps the tree it came from alive, here a temporary
    orphan = rdDASHTree.DASHTree(FIXTURE, ['result']).GetRoot(O)[1]
    self.assertAlmostEqual(orphan.GetValue('result'), -0.99)


class TestTies(unittest.TestCase):
  """Atoms that fit the same node, on the pruned tree in Data/DASHTree."""

  @classmethod
  def setUpClass(cls):
    cls.tree = rdDASHTree.DASHTree(SHIPPED_TREE, ['result', 'std'])
    cls.legacy = rdDASHTree.DASHParams()
    cls.legacy.legacyTieBreaking = True

  def classValues(self, mol, params=None):
    """raw values by symmetry class, which any atom order agrees on"""
    values = self.tree.GetMolProperty(mol, 'result', params or rdDASHTree.DASHParams())
    classes = {}
    for atom, rank in enumerate(Chem.CanonicalRankAtoms(mol, breakTies=False)):
      classes.setdefault(rank, []).append(values[atom])
    return {rank: sorted(v) for rank, v in classes.items()}

  def testAtomOrderDoesNotMatter(self):
    # the two CH2 next to 5-aminoisoindoline's ring NH look alike for five
    # levels, and under the legacy rule the SMILES decides which comes first
    first, second = (molWithHs(s) for s in ('Nc1ccc2CNCc2c1', 'Nc1ccc2c(c1)CNC2'))
    self.assertEqual(self.classValues(first), self.classValues(second))
    self.assertNotEqual(self.classValues(first, self.legacy), self.classValues(second, self.legacy))
    params = Chem.SmilesParserParams()
    params.removeHs = False
    sdf = os.path.join(TEST_DATA_DIR, 'dash_ref_mols.sdf')
    for mol in list(Chem.SDMolSupplier(sdf, removeHs=False))[:20]:
      expected = self.classValues(mol)
      for smiles in Chem.MolToRandomSmilesVect(mol, 3, randomSeed=42):
        self.assertEqual(self.classValues(Chem.MolFromSmiles(smiles, params)), expected, smiles)

  def testTheBetterAttendedChildWins(self):
    # the legacy descents of the two SMILES part at two children of one node;
    # both now take the earlier child, which has the higher attention
    first, second = (molWithHs(s) for s in ('Nc1ccc2CNCc2c1', 'Nc1ccc2c(c1)CNC2'))
    paths = [
      list(self.tree.GetAtomNodePath(first, 6, self.legacy)),
      list(self.tree.GetAtomNodePath(second, 8, self.legacy))
    ]
    level = next(i for i, (a, b) in enumerate(zip(*paths)) if a != b)
    earlier, later = sorted(paths, key=lambda path: path[level])
    self.assertGreaterEqual(
      self.tree.GetNode(earlier[0], earlier[level]).GetAttention(),
      self.tree.GetNode(later[0], later[level]).GetAttention())
    self.assertEqual(list(self.tree.GetAtomNodePath(first, 6)), earlier)
    self.assertEqual(list(self.tree.GetAtomNodePath(second, 8)), earlier)


class TestDASHTree(unittest.TestCase):

  @classmethod
  def setUpClass(cls):
    cls.tree = rdDASHTree.DASHTree(SHIPPED_TREE, ['result', 'std'])

  def testIntrospection(self):
    self.assertGreater(self.tree.GetNumNodes(), 0)
    self.assertEqual(sorted(self.tree.GetPropertyNames()), ['result', 'std'])
    self.assertIn('result', self.tree.GetAvailablePropertyNames())
    self.assertEqual(self.tree.GetFilename(), SHIPPED_TREE)

  def testUnnormalisedIsRaw(self):
    mol = molWithHs('CC(=O)Nc1ccc(O)cc1')
    options = rdDASHTree.NormalizationOptions()
    options.normalization = rdDASHTree.Normalization.NONE
    details = self.tree.GetPartialChargesDetails(mol, options)
    self.assertEqual(list(details['charges']), list(details['raw']))
    self.assertEqual(list(details['raw']), list(self.tree.GetMolProperty(mol, 'result')))
    self.assertTrue(all(s > 0 for s in details['std']))
    self.assertTrue(all(2 <= d <= 17 for d in details['matchDepth']))

  def testSymmetryEquivalentAtoms(self):
    charges = self.tree.GetPartialCharges(molWithHs('c1ccccc1'))
    self.assertEqual(len(set(charges[:6])), 1)
    self.assertEqual(len(set(charges[6:])), 1)

  def testMaxDepth(self):
    mol = molWithHs('CC(C)Cc1ccc(cc1)C(C)C(=O)O')
    previous = 0
    for maxDepth in (0, 1, 2, 3, 5, 16):
      params = rdDASHTree.DASHParams()
      params.maxDepth = maxDepth
      deepest = 0
      for i in range(mol.GetNumAtoms()):
        path = self.tree.GetAtomNodePath(mol, i, params)
        # a hydrogen spends its first node reaching the heavy atom, so two
        # nodes are always possible however small maxDepth is
        self.assertLessEqual(len(path), 1 + max(2, maxDepth))
        deepest = max(deepest, len(path))
      self.assertGreaterEqual(deepest, previous)
      previous = deepest

  def testStopFlagAgreesWithTheAttentions(self):
    # the default path takes a precomputed flag; nudging the threshold off its
    # default forces the arithmetic, and the two must agree
    for smiles in ('CC(C)Cc1ccc(cc1)C(C)C(=O)O', 'CN(C)CCOC(c1ccccc1)c1ccccc1'):
      mol = molWithHs(smiles)
      slow = rdDASHTree.NormalizationOptions()
      slow.params.attentionThreshold = 10.000000001
      self.assertEqual(list(self.tree.GetPartialCharges(mol)),
                       list(self.tree.GetPartialCharges(mol, slow)))

  def testAtomWithNoDataIsRefused(self):
    # this bridging [N+] classifies into a branch holding one node with no value
    mol = molWithHs('COC(=O)C1=CO[C@@H](C)[C@H]2C=[N+]3CCc4c([nH]c5ccccc45)[C@@H]3C[C@H]12')
    with self.assertRaises(ValueError):
      self.tree.GetPartialCharges(mol)
    self.assertTrue(math.isnan(self.tree.GetAtomProperty(mol, 11, 'result')))

  def testStalePropertyCacheIsRefused(self):
    mol = Chem.RWMol()
    mol.AddAtom(Chem.Atom(6))
    mol.AddAtom(Chem.Atom(8))
    mol.AddBond(0, 1, Chem.BondType.SINGLE)
    with self.assertRaises(ValueError):
      self.tree.GetPartialCharges(mol)

  def testOutOfRangeAtomIndex(self):
    mol = molWithHs('CCO')
    with self.assertRaises(IndexError):
      self.tree.GetAtomNodePath(mol, mol.GetNumAtoms())
    with self.assertRaises(IndexError):
      self.tree.GetAtomProperty(mol, mol.GetNumAtoms(), 'result')
    with self.assertRaises(IndexError):
      self.tree.GetMatchedSubstructure(mol, mol.GetNumAtoms())


class TestBatch(unittest.TestCase):

  @classmethod
  def setUpClass(cls):
    cls.tree = rdDASHTree.DASHTree(SHIPPED_TREE, ['result', 'std'])
    cls.mols = [
      molWithHs(s)
      for s in ('CCO', 'CC(=O)Nc1ccc(O)cc1', 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
                'CN1C=NC2=C1C(=O)N(C)C(=O)N2C', 'CC(=O)Oc1ccccc1C(=O)O', 'c1ccc2c(c1)[nH]c1ccccc12',
                'OC(=O)CCc1c[nH]c2ccccc12', 'CN(C)CCOC(c1ccccc1)c1ccccc1')
    ]

  def testEveryThreadCountAgrees(self):
    single = [list(self.tree.GetPartialCharges(m)) for m in self.mols]
    for numThreads in (1, 2, 4, 8, 0, -1):
      batch = self.tree.GetPartialChargesBatch(self.mols, rdDASHTree.NormalizationOptions(),
                                               numThreads)
      self.assertEqual([list(r) for r in batch], single, msg=str(numThreads))

  def testEmptyBatch(self):
    self.assertEqual(
      len(self.tree.GetPartialChargesBatch([], rdDASHTree.NormalizationOptions(), 4)), 0)

  def testPropertyBatch(self):
    batch = self.tree.GetMolPropertyBatch(self.mols, 'result', rdDASHTree.DASHParams(), 4)
    self.assertEqual(len(batch), len(self.mols))
    for row, mol in zip(batch, self.mols):
      self.assertEqual(list(row), list(self.tree.GetMolProperty(mol, 'result')))

  def testAcceptsAnyPythonSequence(self):
    asTuple = self.tree.GetPartialChargesBatch(tuple(self.mols))
    asList = self.tree.GetPartialChargesBatch(self.mols)
    self.assertEqual([list(r) for r in asTuple], [list(r) for r in asList])


@needsPublishedTree
class TestReferenceAgreement(unittest.TestCase):
  """Exact agreement with the values the DASH-tree python package produces."""

  def testAgainstStoredReference(self):
    sdf = os.path.join(TEST_DATA_DIR, 'dash_ref_mols.sdf')
    values = os.path.join(TEST_DATA_DIR, 'dash_ref_values.txt')

    tree = rdDASHTree.DASHTree(PUBLISHED_TREE, ['result', 'std'])
    mols = [m for m in Chem.SDMolSupplier(sdf, removeHs=False) if m is not None]
    self.assertGreater(len(mols), 0)

    cache = {}
    for i, mol in enumerate(mols):
      options = rdDASHTree.NormalizationOptions()
      # the package takes the first of several atoms that fit a node
      options.params.legacyTieBreaking = True
      options.normalization = rdDASHTree.Normalization.NONE
      details = tree.GetPartialChargesDetails(mol, options)
      entry = dict(details)
      for name, norm in (('symmetric', rdDASHTree.Normalization.SYMMETRIC),
                         ('std_weighted', rdDASHTree.Normalization.STD_WEIGHTED)):
        options.normalization = norm
        entry[name] = tree.GetPartialCharges(mol, options)
      cache[i] = entry

    nChecked = 0
    with open(values) as f:
      for line in f:
        if line.startswith('#'):
          continue
        parts = line.split()
        if not parts:
          continue
        molIdx, atomIdx = int(parts[0]), int(parts[1])
        rawResult, stdUsed = float(parts[2]), float(parts[3])
        depth = int(parts[4])
        qNone, qSym, qStd = float(parts[5]), float(parts[6]), float(parts[7])
        entry = cache[molIdx]
        # the stored values are table lookups, so these must be exact
        self.assertEqual(entry['raw'][atomIdx], rawResult)
        self.assertEqual(entry['std'][atomIdx], stdUsed)
        self.assertEqual(entry['matchDepth'][atomIdx], depth)
        self.assertEqual(entry['charges'][atomIdx], qNone)
        # the normalised ones go through a division; allow the last bit
        self.assertAlmostEqual(entry['symmetric'][atomIdx], qSym, delta=1e-15)
        self.assertAlmostEqual(entry['std_weighted'][atomIdx], qStd, delta=1e-15)
        nChecked += 1
    self.assertGreater(nChecked, 0)


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
