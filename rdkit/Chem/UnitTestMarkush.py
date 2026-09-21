#  Copyright (C) 2026 RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.

import unittest

from rdkit import Chem
from rdkit.Chem import Markush


class TestCase(unittest.TestCase):

  def testGenericScopeUsesGenericMatchersByDefault(self):
    query = Chem.MolFromSmarts('OC* |$;;ARY$|')
    Chem.SetGenericQueriesFromProperties(query)

    self.assertTrue(
        Markush.IsInMarkushScope(query, Chem.MolFromSmiles('c1ccccc1CO')))
    self.assertFalse(
        Markush.IsInMarkushScope(query, Chem.MolFromSmiles('C1CCCCC1CO')))

  def testExplicitParametersOverrideDefaultMatchingPolicy(self):
    query = Chem.MolFromSmarts('OC* |$;;ARY$|')
    Chem.SetGenericQueriesFromProperties(query)
    parameters = Chem.SubstructMatchParameters()

    self.assertTrue(
        Markush.IsInMarkushScope(query, Chem.MolFromSmiles('C1CCCCC1CO'),
                                 parameters))

  def testEnumerationFiltersAndDeduplicatesInInputOrder(self):
    query = Chem.MolFromSmarts('OC* |$;;ARY$|')
    Chem.SetGenericQueriesFromProperties(query)
    candidates = [
        Chem.MolFromSmiles('c1ccccc1CO'),
        Chem.MolFromSmiles('C1CCCCC1CO'),
        Chem.MolFromSmiles('OCC1=CC=CC=C1'),
    ]

    result = Markush.EnumerateMarkush(query, candidates)

    self.assertEqual(len(result), 1)
    self.assertEqual(Chem.MolToSmiles(result[0]), 'OCc1ccccc1')

  def testEnumerationHonorsExplicitParameters(self):
    query = Chem.MolFromSmarts('OC* |$;;ARY$|')
    Chem.SetGenericQueriesFromProperties(query)
    parameters = Chem.SubstructMatchParameters()

    result = Markush.EnumerateMarkush(
        query, (Chem.MolFromSmiles('C1CCCCC1CO'), ), parameters)

    self.assertEqual([Chem.MolToSmiles(mol) for mol in result], ['OCC1CCCCC1'])

  def testEnumerationHandlesEmptyAndNonMatchingCandidateSets(self):
    query = Chem.MolFromSmarts('OC* |$;;ARY$|')
    Chem.SetGenericQueriesFromProperties(query)

    self.assertEqual(Markush.EnumerateMarkush(query, ()), ())
    self.assertEqual(
        Markush.EnumerateMarkush(query, (Chem.MolFromSmiles('C1CCCCC1CO'), )),
        ())

  def testMakeFormulaCoversInputsAndRemovesDuplicates(self):
    ethanol = Chem.MolFromSmiles('CCO')
    benzene = Chem.MolFromSmiles('c1ccccc1')
    formula = Markush.MakeMarkushFormula(
        (ethanol, benzene, Chem.MolFromSmiles('OCC')))

    self.assertEqual(len(formula.queries), 2)
    self.assertIsInstance(formula.queries, tuple)
    self.assertTrue(Markush.IsInMarkushScope(formula, ethanol))
    self.assertTrue(Markush.IsInMarkushScope(formula, benzene))
    self.assertEqual(
        [Chem.MolToSmiles(mol)
         for mol in Markush.EnumerateMarkush(formula, (benzene, ethanol))],
        ['c1ccccc1', 'CCO'])

  def testInvalidInputsAreRejected(self):
    molecule = Chem.MolFromSmiles('CC')
    with self.assertRaisesRegex(ValueError, '^a Markush formula needs at least one query$'):
      Markush.MarkushFormula(())
    with self.assertRaisesRegex(ValueError, '^Markush queries cannot be None$'):
      Markush.MarkushFormula((None, ))
    with self.assertRaisesRegex(ValueError, '^molecule cannot be None$'):
      Markush.IsInMarkushScope(molecule, None)
    with self.assertRaisesRegex(ValueError, '^candidate molecules cannot be None$'):
      Markush.EnumerateMarkush(molecule, (None, ))
    with self.assertRaisesRegex(ValueError, '^molecules cannot contain None$'):
      Markush.MakeMarkushFormula((None, ))
    with self.assertRaisesRegex(ValueError, '^a Markush formula needs at least one query$'):
      Markush.MakeMarkushFormula(())


if __name__ == '__main__':
  unittest.main()
