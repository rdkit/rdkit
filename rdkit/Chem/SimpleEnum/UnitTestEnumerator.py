#
# Created by Peter Gedeck, January 2017
#

import doctest
import os
import unittest

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.SimpleEnum import Enumerator


def load_tests(loader, tests, ignore):
  """ Add the Doctests from the module """
  tests.addTests(doctest.DocTestSuite(Enumerator, optionflags=doctest.ELLIPSIS))
  return tests


class TestCase(unittest.TestCase):

  def test_EnumerateReactionException(self):
    testFile = os.sep.join(
      [os.path.dirname(os.path.abspath(__file__)), 'test_data', 'boronic1.rxn'])
    rxn = AllChem.ReactionFromRxnFile(testFile)
    rxn.Initialize()
    reacts1 = ['Brc1ccccc1', 'Brc1ncccc1', 'Brc1cnccc1']
    reacts1 = [Chem.MolFromSmiles(x) for x in reacts1]

    self.assertRaises(ValueError, Enumerator.EnumerateReaction, rxn, (reacts1, ))

  def test_UniquifyingEnhancedStereo(self):
    rxn = AllChem.ReactionFromSmarts('[C:1]-[O:2]>>[C:1]-[N:2]')
    rxn.Initialize()
    reacts1 = ['C[C@@H](F)CO','C[C@@H](F)CO |a:1|','C[C@@H](F)CO |o1:1|','C[C@@H](F)CO |&1:1|']
    reacts1 = [Chem.MolFromSmiles(x) for x in reacts1]

    prods = list(Enumerator.EnumerateReaction(rxn,(reacts1,),uniqueProductsOnly=True))
    ps = Chem.SmilesWriteParams()
    print([Chem.MolToCXSmiles(tpl[0],ps,Chem.CXSmilesFields.CX_ENHANCEDSTEREO) for tpl in prods])
    self.assertEqual(len(prods),len(reacts1))


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
