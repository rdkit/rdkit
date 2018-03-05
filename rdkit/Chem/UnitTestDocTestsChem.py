from __future__ import print_function
import unittest
import doctest
from rdkit import Chem
from rdkit.Chem import MCS, FragmentMatcher, MACCSkeys, Descriptors, TemplateAlign
from rdkit.Chem import Recap, BRICS, AllChem, PropertyMol, SaltRemover, EnumerateHeterocycles, EnumerateStereoisomers

def load_tests(loader, tests, ignore):  # pylint: disable=unused-argument
  """ Add the Doctests from the module """
  tests.addTests(doctest.DocTestSuite(EnumerateStereoisomers, optionflags=doctest.ELLIPSIS))
  tests.addTests(doctest.DocTestSuite(EnumerateHeterocycles, optionflags=doctest.ELLIPSIS))
  tests.addTests(doctest.DocTestSuite(MCS, optionflags=doctest.ELLIPSIS))
  tests.addTests(doctest.DocTestSuite(FragmentMatcher, optionflags=doctest.ELLIPSIS))
  tests.addTests(doctest.DocTestSuite(MACCSkeys, optionflags=doctest.ELLIPSIS))
  tests.addTests(doctest.DocTestSuite(Descriptors, optionflags=doctest.ELLIPSIS))
  tests.addTests(doctest.DocTestSuite(TemplateAlign, optionflags=doctest.ELLIPSIS))
  tests.addTests(doctest.DocTestSuite(Recap, optionflags=doctest.ELLIPSIS))
  tests.addTests(doctest.DocTestSuite(BRICS, optionflags=doctest.ELLIPSIS))
  tests.addTests(doctest.DocTestSuite(AllChem, optionflags=doctest.ELLIPSIS))
  tests.addTests(doctest.DocTestSuite(PropertyMol, optionflags=doctest.ELLIPSIS))
  tests.addTests(doctest.DocTestSuite(SaltRemover, optionflags=doctest.ELLIPSIS))
  tests.addTests(doctest.DocTestSuite(Chem, optionflags=doctest.ELLIPSIS))
  return tests

if __name__ == '__main__':  # pragma: nocover
  unittest.main()
