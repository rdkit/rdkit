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


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
