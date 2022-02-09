"""
unit tests for the MolStandardize module
tests include:
reorder_tautomers
"""

import unittest

from rdkit import Chem
from rdkit.Chem import MolStandardize
from rdkit.Chem.MolStandardize import rdMolStandardize

class TestCase(unittest.TestCase):

    def testBasic(self):
        m = Chem.MolFromSmiles('Oc1c(cccc3)c3nc2ccncc12')
        enumerator = rdMolStandardize.TautomerEnumerator()
        canon = enumerator.Canonicalize(m)
        reord = MolStandardize.ReorderTautomers(m)[0]
        canonSmile = Chem.MolToSmiles(canon)
        reordSmile = Chem.MolToSmiles(reord)
        self.assertEquals(canonSmile, reordSmile)

    def testLength(self):
        m = Chem.MolFromSmiles('Oc1c(cccc3)c3nc2ccncc12')
        enumerator = rdMolStandardize.TautomerEnumerator()
        tauts = enumerator.Enumerate(m)
        reordtauts = MolStandardize.ReorderTautomers(m)
        self.assertEquals(len(reordtauts), len(tauts))
