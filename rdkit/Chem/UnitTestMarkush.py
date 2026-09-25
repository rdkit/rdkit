#  Copyright (C) 2026 RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.

import time
import unittest

from rdkit import Chem
from rdkit.Chem import Markush


class TestCase(unittest.TestCase):

  @staticmethod
  def _generic_query(smarts):
    query = Chem.MolFromSmarts(smarts)
    Chem.SetGenericQueriesFromProperties(query)
    return query

  @staticmethod
  def _branched_scaffold(index):
    return 'C' + ''.join(
        'C(C)' if index >> bit & 1 else 'C' for bit in range(7))

  @classmethod
  def _documented_generic_formula_cases(cls):
    # The group codes come from the RDKit Book's Beilstein/Reaxys generic
    # group table. Every tail below is a member of the associated group.
    tails_by_group = {
        'ALK': ('C', 'CC', 'CCC', 'C(C)C', 'C(C)(C)C'),
        'AEL': ('C=C', 'C=CC', 'C=CCC', 'C=CCCC', 'C=CCCCC'),
        'AYL': ('C#C', 'C#CC', 'C#CCC', 'C#CCCC', 'C#CCCCC'),
        'AOX': ('OC', 'OCC', 'OCCC', 'OCCCC', 'OCCCCC'),
        'CAL': ('C1CC1', 'C1CCC1', 'C1CCCC1', 'C1CCCCC1', 'C1CCCCCC1'),
        'CEL': ('c1ccccc1', 'C1=CCCC1', 'C1=CCCCC1', 'C1=CCCCCC1',
                'C1=CCCCCCC1'),
        'HAR': ('c1ccncc1', 'c1ccoc1', 'c1ccsc1', 'c1cn[nH]c1',
                'c1ncc[nH]1'),
    }
    for group, tails in tails_by_group.items():
      for index in range(16):
        scaffold = cls._branched_scaffold(index)
        atom_count = Chem.MolFromSmarts(scaffold + '*').GetNumAtoms() - 1
        smarts = scaffold + '* |$' + ';' * atom_count + group + '$|'
        yield smarts, tuple(scaffold + tail for tail in tails)

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

  def testFormulaUsesAllDisjunctiveAlternatives(self):
    formula = Markush.MarkushFormula(
        (Chem.MolFromSmarts('CCN'), Chem.MolFromSmarts('CCO')))
    amine = Chem.MolFromSmiles('CCN')
    alcohol = Chem.MolFromSmiles('CCO')
    neither = Chem.MolFromSmiles('CCC')

    self.assertTrue(Markush.IsInMarkushScope(formula, amine))
    self.assertTrue(Markush.IsInMarkushScope(formula, alcohol))
    self.assertFalse(Markush.IsInMarkushScope(formula, neither))
    self.assertEqual(
        Markush.EnumerateMarkush(formula, (neither, alcohol, amine)),
        (alcohol, amine))

  def testExplicitParametersAreNotMutated(self):
    query = self._generic_query('OC* |$;;ARY$|')
    parameters = Chem.SubstructMatchParameters()

    Markush.IsInMarkushScope(query, Chem.MolFromSmiles('C1CCCCC1CO'), parameters)

    self.assertFalse(parameters.useGenericMatchers)
    self.assertFalse(parameters.useChirality)

  def testInvalidParametersAreRejectedBeforeMatchingOrEnumeration(self):
    molecule = Chem.MolFromSmiles('CCO')
    formula = Markush.MakeMarkushFormula((molecule, ))

    with self.assertRaisesRegex(TypeError,
                                '^params must be SubstructMatchParameters or None$'):
      Markush.IsInMarkushScope(formula, molecule, params='not parameters')
    with self.assertRaisesRegex(TypeError,
                                '^params must be SubstructMatchParameters or None$'):
      Markush.EnumerateMarkush(formula, (), params='not parameters')

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

  def testOneHundredDocumentedMarkushFormulaeAndTheirCompounds(self):
    formula_count = 0
    compound_count = 0
    formulae = set()
    identities = set()
    max_enumeration_microseconds = 0
    for smarts, candidate_smiles in self._documented_generic_formula_cases():
      with self.subTest(smarts=smarts):
        query = self._generic_query(smarts)
        candidates = tuple(Chem.MolFromSmiles(smiles) for smiles in candidate_smiles)
        self.assertNotIn(None, candidates)

        started = time.perf_counter_ns()
        enumerated = Markush.EnumerateMarkush(query, candidates)
        elapsed_microseconds = (time.perf_counter_ns() - started) // 1000
        max_enumeration_microseconds = max(max_enumeration_microseconds,
                                           elapsed_microseconds)

        self.assertEqual(enumerated, candidates)
        self.assertTrue(all(Markush.IsInMarkushScope(query, molecule)
                            for molecule in candidates))
        exact_formula = Markush.MakeMarkushFormula(enumerated)
        self.assertEqual(Markush.EnumerateMarkush(exact_formula, candidates),
                         candidates)
        formula_count += 1
        compound_count += len(candidates)
        formulae.add(smarts)
        identities.update(Chem.MolToSmiles(molecule) for molecule in candidates)

    self.assertGreaterEqual(formula_count, 100)
    self.assertGreaterEqual(len(formulae), 100)
    self.assertGreaterEqual(compound_count, 500)
    self.assertGreaterEqual(len(identities), 500)
    # Each individual five-member finite library is deliberately small enough
    # to be enumerated within 100,000 microseconds on the test host.
    self.assertLess(max_enumeration_microseconds, 100000)

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

  def testDirectQueriesAreSubstructureQueriesButMadeFormulaeAreExact(self):
    ethanol_query = Chem.MolFromSmarts('CCO')
    propanol = Chem.MolFromSmiles('CCCO')

    self.assertTrue(Markush.IsInMarkushScope(ethanol_query, propanol))
    formula = Markush.MakeMarkushFormula((Chem.MolFromSmiles('CCO'), ))
    self.assertFalse(Markush.IsInMarkushScope(formula, propanol))

  def testMakeFormulaTreatsIsotopesAndStereoisomersAsDistinctIdentities(self):
    ethanol = Chem.MolFromSmiles('CCO')
    labeled_ethanol = Chem.MolFromSmiles('[13CH3]CO')
    first_enantiomer = Chem.MolFromSmiles('C[C@H](F)Cl')
    second_enantiomer = Chem.MolFromSmiles('C[C@@H](F)Cl')
    formula = Markush.MakeMarkushFormula(
        (ethanol, labeled_ethanol, first_enantiomer, second_enantiomer))

    self.assertEqual(len(formula.queries), 4)
    self.assertEqual(
        Markush.EnumerateMarkush(
            formula,
            (ethanol, labeled_ethanol, first_enantiomer, second_enantiomer)),
        (ethanol, labeled_ethanol, first_enantiomer, second_enantiomer))

  def testMakeFormulaSnapshotsTheInputChemicalIdentity(self):
    molecule = Chem.MolFromSmiles('CCO')
    original = Chem.Mol(molecule)
    formula = Markush.MakeMarkushFormula((molecule, ))
    molecule.GetAtomWithIdx(0).SetAtomicNum(7)

    self.assertTrue(Markush.IsInMarkushScope(formula, original))
    self.assertFalse(Markush.IsInMarkushScope(formula, molecule))

  def testEnumerationPreservesFirstEquivalentCandidateObject(self):
    query = Chem.MolFromSmarts('CCO')
    first = Chem.MolFromSmiles('OCC')
    duplicate = Chem.MolFromSmiles('CCO')

    result = Markush.EnumerateMarkush(query, (first, duplicate))

    self.assertEqual(len(result), 1)
    self.assertIs(result[0], first)

  def testGeneratorInputsAreAccepted(self):
    ethanol = Chem.MolFromSmiles('CCO')
    benzene = Chem.MolFromSmiles('c1ccccc1')
    formula = Markush.MakeMarkushFormula(molecule for molecule in (ethanol, benzene))

    result = Markush.EnumerateMarkush(
        formula, (molecule for molecule in (benzene, ethanol)))

    self.assertEqual(result, (benzene, ethanol))

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
    with self.assertRaisesRegex(TypeError, '^Markush queries must be iterable$'):
      Markush.MarkushFormula(None)
    with self.assertRaisesRegex(ValueError, '^a Markush formula needs at least one query$'):
      Markush.MarkushFormula(())
    with self.assertRaisesRegex(ValueError, '^Markush queries cannot be None$'):
      Markush.MarkushFormula((None, ))
    with self.assertRaisesRegex(TypeError, '^Markush queries must be RDKit molecules$'):
      Markush.MarkushFormula(('CC', ))
    with self.assertRaisesRegex(ValueError, '^molecule cannot be None$'):
      Markush.IsInMarkushScope(molecule, None)
    with self.assertRaisesRegex(TypeError, '^molecule must be an RDKit molecule$'):
      Markush.IsInMarkushScope(molecule, 'CC')
    with self.assertRaisesRegex(ValueError, '^candidate molecules cannot be None$'):
      Markush.EnumerateMarkush(molecule, (None, ))
    with self.assertRaisesRegex(ValueError, '^candidate molecules cannot be None$'):
      Markush.EnumerateMarkush(molecule, (molecule, None))
    with self.assertRaisesRegex(TypeError, '^candidate molecules must be RDKit molecules$'):
      Markush.EnumerateMarkush(molecule, ('CC', ))
    with self.assertRaisesRegex(TypeError, '^candidates must be iterable$'):
      Markush.EnumerateMarkush(molecule, None)
    with self.assertRaisesRegex(ValueError, '^molecules cannot contain None$'):
      Markush.MakeMarkushFormula((None, ))
    with self.assertRaisesRegex(TypeError, '^molecules must be RDKit molecules$'):
      Markush.MakeMarkushFormula(('CC', ))
    with self.assertRaisesRegex(TypeError, '^molecules must be iterable$'):
      Markush.MakeMarkushFormula(None)
    with self.assertRaisesRegex(ValueError, '^a Markush formula needs at least one query$'):
      Markush.MakeMarkushFormula(())
    with self.assertRaisesRegex(ValueError, '^Markush queries cannot be None$'):
      Markush.EnumerateMarkush(None, ())


if __name__ == '__main__':
  unittest.main()
