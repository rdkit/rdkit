#
#  Copyright (C) 2026 Marc Lehner and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""Tests for rdkit.Chem.DASHTree.

A tree is a few hundred megabytes of published data and is not in the
repository, so these run only when DASH_TREE_FILE points at a converted
container and are skipped otherwise.
"""

import os
import unittest

from rdkit import Chem
from rdkit.Chem import DASHTree, rdDASHTree

treeFile = os.environ.get("DASH_TREE_FILE")


@unittest.skipUnless(treeFile and os.path.isfile(treeFile),
                     "set DASH_TREE_FILE to a converted container")
class TestCase(unittest.TestCase):

  def setUp(self):
    self.tree = rdDASHTree.DASHTree(treeFile, ["result", "std"])
    self.mol = Chem.AddHs(Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1"))

  def testNodes(self):
    path = self.tree.GetAtomNodePath(self.mol, 0)
    root = self.tree.GetRoot(path[0])
    self.assertEqual(root.GetId(), path[1])
    self.assertEqual(root.GetAtomFeatureIndex(), path[0])
    self.assertEqual((root.GetConAtom(), root.GetConType()), (-1, -1))
    # every node of the path is a child of the one before it
    node = root
    for nodeId in path[2:]:
      children = [node.GetChild(i).GetId() for i in range(node.GetNumChildren())]
      self.assertIn(nodeId, children)
      node = self.tree.GetNode(path[0], nodeId)
      self.assertGreaterEqual(node.GetConAtom(), 0)
    with self.assertRaises(ValueError):
      node.GetChild(node.GetNumChildren())
    # a node is a sequence of its children
    self.assertEqual(len(root), root.GetNumChildren())
    self.assertEqual([c.GetId() for c in root],
                     [root.GetChild(i).GetId() for i in range(len(root))])
    self.assertEqual(root[-1].GetId(), root.GetChild(len(root) - 1).GetId())
    with self.assertRaises(IndexError):
      root[len(root)]
    # the deepest node carrying a value is what the atom is assigned
    values = [self.tree.GetNode(path[0], i).GetValue("result") for i in path[1:]]
    deepest = next(v for v in reversed(values) if v == v)
    self.assertEqual(deepest, self.tree.GetAtomProperty(self.mol, 0, "result"))

  def testPatternMatchesItsOwnAtom(self):
    for atom in self.mol.GetAtoms():
      smarts = DASHTree.GetAtomMatchSmarts(self.tree, self.mol, atom.GetIdx())
      query = Chem.MolFromSmarts(smarts)
      self.assertIsNotNone(query, smarts)
      mapped = next(a.GetIdx() for a in query.GetAtoms() if a.GetAtomMapNum() == 1)
      landings = {m[mapped] for m in self.mol.GetSubstructMatches(query, uniquify=False)}
      self.assertIn(atom.GetIdx(), landings, smarts)

  def testMatchedSubstructure(self):
    path = self.tree.GetAtomNodePath(self.mol, 0)
    atoms = self.tree.GetMatchedSubstructure(self.mol, 0)
    self.assertEqual(len(atoms), len(path) - 1)
    self.assertEqual(atoms[0], 0)
    # the matched atoms are, position for position, one of the ways the
    # pattern fits the molecule
    query = DASHTree.NodePathToQueryMol(self.tree, path, foldHydrogens=False)
    fits = self.mol.GetSubstructMatches(query, uniquify=False, maxMatches=100000)
    self.assertIn(tuple(atoms), fits)


if __name__ == '__main__':
  unittest.main()
