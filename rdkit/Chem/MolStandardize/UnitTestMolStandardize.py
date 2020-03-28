"""
unit tests for the MolStandardize module
tests include:
reorder_tautomers
"""


import unittest


from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import Draw


from Chem.MolStandardize import reorder_tautomers


class TestCase(unittest.TestCase):

    def testBasic(self):
        enumerator = rdMolStandardize.TautomerEnumerator()
        m = Chem.MolFromSmiles('Oc1c(cccc3)c3nc2ccncc12')
        canon = enumerator.Canonicalize(m)
        reord = reorder_tautomers(m)[0]
        canonSmile = Chem.MolToSmiles(canon)
        reordSmile = Chem.MolToSmiles(reord)
        self.assertEquals(canonSmile, reordSmile)

    



