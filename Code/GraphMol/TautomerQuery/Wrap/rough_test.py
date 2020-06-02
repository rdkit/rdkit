#
# Created by Gareth Jones on 5/30/2020.
#
# Copyright 2020 Schrodinger, Inc
#

""" A rough test of Tautomer Queries 

Rough in that only basic functionality is evaluated.
"""

from rdkit import Chem, DataStructs
from rdkit.Chem import rdTautomerQuery
from unittest import TestCase, main

class TautomerQueryTestCase(TestCase):
    
    def test_basic(self):
        mol = Chem.MolFromSmiles("O=C1CCCCC1")
        tautomer_query = rdTautomerQuery.TautomerQuery(mol)

        target = Chem.MolFromSmiles("OC1=CCCC(CC)C1")
        match = tautomer_query.IsSubstructOf(target)
        self.assertTrue(match)
        match = tautomer_query.GetSubstructMatch(target);
        self.assertEqual(7, len(match))

        matches = tautomer_query.GetSubstructMatches(target);
        self.assertEqual(1, len(matches))

        tautomer_matches = tautomer_query.GetSubstructMatchesWithTautomers(target)
        self.assertEqual(1, len(tautomer_matches))
        self.assertEqual(len(tautomer_matches[0][0]), 7)
        matching_tautomer = tautomer_matches[0][1];
        self.assertTrue(target.HasSubstructMatch(matching_tautomer))

        tautomer_fingerprint = tautomer_query.PatternFingerprintTemplate()
        target_fingerprint = rdTautomerQuery.PatternFingerprintTautomerTarget(target)
        matching_fingerprints = DataStructs.AllProbeBitsMatch(tautomer_fingerprint, target_fingerprint)
        self.assertTrue(matching_fingerprints)
        
        del tautomer_fingerprint
        del target_fingerprint




if __name__ == "__main__":
    main();