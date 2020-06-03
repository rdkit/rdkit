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
import os

class TautomerQueryTestCase(TestCase):
    
    def test_basic(self):
        mol = Chem.MolFromSmiles("O=C1CCCCC1")
        tautomer_query = rdTautomerQuery.TautomerQuery(mol)
        self.assertEqual(2, len(tautomer_query.GetTautomers()))
        modified_bonds = tautomer_query.GetModifiedBonds()
        self.assertEqual(3, len(modified_bonds));
        modified_atoms = tautomer_query.GetModifiedAtoms()
        self.assertEqual(3, len(modified_atoms));

        target = Chem.MolFromSmiles("OC1=CCCC(CC)C1")
        self.assertTrue(target.HasSubstructMatch(tautomer_query.GetTemplateMolecule()))

        match = tautomer_query.IsSubstructOf(target)
        self.assertTrue(match)
        self.assertTrue(target.HasSubstructMatch(tautomer_query.GetTemplateMolecule()))
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

        file = os.environ['RDBASE']+'/Data/MolStandardize/tautomerTransforms.in';
        tautomer_query2 = rdTautomerQuery.TautomerQuery(mol, file)
        self.assertEqual(2, len(tautomer_query.GetTautomers()))




if __name__ == "__main__":
    main();