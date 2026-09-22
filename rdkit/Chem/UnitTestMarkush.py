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

  @staticmethod
  def _generic_query(smarts):
    query = Chem.MolFromSmarts(smarts)
    Chem.SetGenericQueriesFromProperties(query)
    return query

  def testGenericScopeUsesGenericMatchersByDefault(self):
    query = self._generic_query('OC* |$;;ARY$|')

    self.assertTrue(
        Markush.IsInMarkushScope(query, Chem.MolFromSmiles('c1ccccc1CO')))
    self.assertFalse(
        Markush.IsInMarkushScope(query, Chem.MolFromSmiles('C1CCCCC1CO')))

  def testExplicitParametersOverrideDefaultMatchingPolicy(self):
    query = self._generic_query('OC* |$;;ARY$|')
    parameters = Chem.SubstructMatchParameters()

    self.assertTrue(
        Markush.IsInMarkushScope(query, Chem.MolFromSmiles('C1CCCCC1CO'),
                                 parameters))

  def testEnumerationFiltersAndDeduplicatesInInputOrder(self):
    query = self._generic_query('OC* |$;;ARY$|')
    candidates = [
        Chem.MolFromSmiles('c1ccccc1CO'),
        Chem.MolFromSmiles('C1CCCCC1CO'),
        Chem.MolFromSmiles('OCC1=CC=CC=C1'),
    ]

    result = Markush.EnumerateMarkush(query, candidates)

    self.assertEqual(len(result), 1)
    self.assertEqual(Chem.MolToSmiles(result[0]), 'OCc1ccccc1')

  def testEnumerationHonorsExplicitParameters(self):
    query = self._generic_query('OC* |$;;ARY$|')
    parameters = Chem.SubstructMatchParameters()

    result = Markush.EnumerateMarkush(
        query, (Chem.MolFromSmiles('C1CCCCC1CO'), ), parameters)

    self.assertEqual([Chem.MolToSmiles(mol) for mol in result], ['OCC1CCCCC1'])

  def testEnumerationHandlesEmptyAndNonMatchingCandidateSets(self):
    query = self._generic_query('OC* |$;;ARY$|')

    self.assertEqual(Markush.EnumerateMarkush(query, ()), ())
    self.assertEqual(
        Markush.EnumerateMarkush(query, (Chem.MolFromSmiles('C1CCCCC1CO'), )),
        ())

  def testDocumentedGenericGroupExamplesAndEnumeration(self):
    # These Beilstein/Reaxys generic-group meanings are documented in the RDKit
    # Book, "Generic (Markush) queries in substructure matching" section.
    cases = (
        ('O* |$;ALK$|', 'CCO', 'O=C=O'),  # ethanol: alkyl; carbon dioxide: not alkyl
        ('C* |$;AEL$|', 'CC=C', 'CCC'),  # propene: alkenyl; propane: not alkenyl
        ('C* |$;AYL$|', 'CC#C', 'CCC'),  # propyne: alkynyl; propane: not alkynyl
        ('C* |$;AOX$|', 'COC', 'CCC'),  # dimethyl ether: alkoxy; propane: not alkoxy
        ('C* |$;CAL$|', 'CC1CCCCC1', 'Cc1ccccc1'),
        # ethylcyclohexane; ethylbenzene
        ('C* |$;CEL$|', 'CC1=CC=CC=C1', 'CC1CCCCC1'),
        # ethylbenzene; ethylcyclohexane
        ('C* |$;HAR$|', 'Cc1ccncc1', 'Cc1ccccc1'),
        # methylpyridine; ethylbenzene
    )
    for smarts, matching_smiles, nonmatching_smiles in cases:
      with self.subTest(smarts=smarts):
        query = self._generic_query(smarts)
        matching = Chem.MolFromSmiles(matching_smiles)
        nonmatching = Chem.MolFromSmiles(nonmatching_smiles)

        self.assertTrue(Markush.IsInMarkushScope(query, matching))
        self.assertFalse(Markush.IsInMarkushScope(query, nonmatching))
        self.assertEqual(
            [Chem.MolToSmiles(mol)
             for mol in Markush.EnumerateMarkush(query,
                                                   (nonmatching, matching))],
            [Chem.MolToSmiles(matching)])

  def testMakeFormulaCoversInputsAndRemovesDuplicates(self):
    ethanol = Chem.MolFromSmiles('CCO')
    benzene = Chem.MolFromSmiles('c1ccccc1')
    formula = Markush.MakeMarkushFormula(
        (ethanol, benzene, Chem.MolFromSmiles('OCC')))

    self.assertEqual(len(formula.queries), 2)
    self.assertIsInstance(formula.queries, tuple)
    self.assertTrue(Markush.IsInMarkushScope(formula, ethanol))
    self.assertTrue(Markush.IsInMarkushScope(formula, benzene))
    self.assertFalse(Markush.IsInMarkushScope(formula, Chem.MolFromSmiles('CCCO')))
    self.assertEqual(
        [Chem.MolToSmiles(mol)
         for mol in Markush.EnumerateMarkush(formula, (benzene, ethanol))],
        ['c1ccccc1', 'CCO'])

  def testMadeFormulaPreservesChirality(self):
    first = Chem.MolFromSmiles('C[C@H](F)Cl')
    second = Chem.MolFromSmiles('C[C@@H](F)Cl')
    formula = Markush.MakeMarkushFormula((first, ))

    self.assertTrue(Markush.IsInMarkushScope(formula, first))
    self.assertFalse(Markush.IsInMarkushScope(formula, second))

  def testFormulaIdentitiesCannotBeSuppliedByCallers(self):
    query = Chem.MolFromSmarts('CCO')

    with self.assertRaises(TypeError):
      Markush.MarkushFormula((query, ), identities=('CCO', ))

  def testDefaultGenericMatchingPreservesChirality(self):
    query = Chem.MolFromSmarts('C[C@H](F)Cl')
    first = Chem.MolFromSmiles('C[C@H](F)Cl')
    second = Chem.MolFromSmiles('C[C@@H](F)Cl')

    self.assertTrue(Markush.IsInMarkushScope(query, first))
    self.assertFalse(Markush.IsInMarkushScope(query, second))

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
    with self.assertRaisesRegex(ValueError, '^Markush queries cannot be None$'):
      Markush.EnumerateMarkush(None, ())


if __name__ == '__main__':
  unittest.main()
