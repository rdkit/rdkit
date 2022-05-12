#
# Created by Gareth Jones on 5/30/2020.
#
# Copyright 2020-2022 Schrodinger, Inc and other RDKit contributors
#  @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
""" A rough test of Tautomer Queries 

Rough in that only basic functionality is evaluated.
"""

from rdkit import Chem, DataStructs
from rdkit.Chem import rdTautomerQuery
from unittest import TestCase, main
import os
import pickle


class TautomerQueryTestCase(TestCase):

  def test_basic(self):
    mol = Chem.MolFromSmiles("O=C1CCCCC1")
    tautomer_query = rdTautomerQuery.TautomerQuery(mol)
    self.assertEqual(2, len(tautomer_query.GetTautomers()))
    modified_bonds = tautomer_query.GetModifiedBonds()
    self.assertEqual(3, len(modified_bonds))
    modified_atoms = tautomer_query.GetModifiedAtoms()
    self.assertEqual(3, len(modified_atoms))

    target = Chem.MolFromSmiles("OC1=CCCC(CC)C1")
    self.assertTrue(target.HasSubstructMatch(tautomer_query.GetTemplateMolecule()))

    match = tautomer_query.IsSubstructOf(target)
    self.assertTrue(match)
    self.assertTrue(target.HasSubstructMatch(tautomer_query.GetTemplateMolecule()))
    #self.assertTrue(target.HasSubstructMatch(tautomer_query.GetTautomers()[1]))
    match = tautomer_query.GetSubstructMatch(target)
    self.assertEqual(7, len(match))

    matches = tautomer_query.GetSubstructMatches(target)
    self.assertEqual(1, len(matches))

    tautomer_matches = tautomer_query.GetSubstructMatchesWithTautomers(target)
    self.assertEqual(1, len(tautomer_matches))
    self.assertEqual(len(tautomer_matches[0][0]), 7)
    matching_tautomer = tautomer_matches[0][1]
    self.assertTrue(target.HasSubstructMatch(matching_tautomer))

    tautomer_fingerprint = tautomer_query.PatternFingerprintTemplate()
    target_fingerprint = rdTautomerQuery.PatternFingerprintTautomerTarget(target)
    matching_fingerprints = DataStructs.AllProbeBitsMatch(tautomer_fingerprint, target_fingerprint)
    self.assertTrue(matching_fingerprints)

    file = os.environ['RDBASE'] + '/Code/GraphMol/MolStandardize/test_data/tautomerTransforms.in'
    tautomer_query2 = rdTautomerQuery.TautomerQuery(mol, file)
    self.assertEqual(2, len(tautomer_query.GetTautomers()))

  def test_parameter_searches(self):
    mol = Chem.MolFromSmiles("O=C1CCCCC1")
    tautomer_query = rdTautomerQuery.TautomerQuery(mol)
    target = Chem.MolFromSmiles("O=C1CCCC(CC)C1")
    params = Chem.SubstructMatchParameters()
    params.uniquify = False
    tautomer_matches = tautomer_query.GetSubstructMatchesWithTautomers(target, params)
    self.assertEqual(2, len(tautomer_matches))
    matches = tautomer_query.GetSubstructMatches(target, params)
    self.assertEqual(2, len(matches))
    match = tautomer_query.GetSubstructMatch(target, params)
    self.assertEqual(7, len(match))
    is_match = tautomer_query.IsSubstructOf(target, params)
    self.assertTrue(is_match)

  def test_types(self):
    mol = Chem.MolFromSmiles("O=C1CCCCC1")
    tautomer_query = rdTautomerQuery.TautomerQuery(mol)
    target = Chem.MolFromSmiles("OC1=CCCC(CC)C1")
    try:
      self.assertTrue(target.HasSubstructMatch(tautomer_query.GetTautomers()[1]))
    except Exception:
      self.fail("Boost type error")

  def test_simple(self):
    mol = Chem.MolFromSmiles("CC=O")
    tautomer_query = rdTautomerQuery.TautomerQuery(mol)
    target = Chem.MolFromSmiles("OC(C)=O")
    tautomers = tautomer_query.GetTautomers()
    self.assertTrue(target.HasSubstructMatch(mol))
    self.assertTrue(tautomer_query.IsSubstructOf(target))

  def test_serialization(self):
    mol = Chem.MolFromSmiles("O=C1CCCCC1")
    base_tautomer_query = rdTautomerQuery.TautomerQuery(mol)
    for tautomer_query in (pickle.loads(pickle.dumps(base_tautomer_query)),
                           rdTautomerQuery.TautomerQuery(base_tautomer_query.ToBinary())):
      self.assertEqual(2, len(tautomer_query.GetTautomers()))
      modified_bonds = tautomer_query.GetModifiedBonds()
      self.assertEqual(3, len(modified_bonds))
      modified_atoms = tautomer_query.GetModifiedAtoms()
      self.assertEqual(3, len(modified_atoms))

      target = Chem.MolFromSmiles("OC1=CCCC(CC)C1")
      self.assertTrue(target.HasSubstructMatch(tautomer_query.GetTemplateMolecule()))


if __name__ == "__main__":
  main()
