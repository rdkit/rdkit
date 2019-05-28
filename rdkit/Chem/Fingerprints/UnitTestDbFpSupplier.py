# $Id$
#
#  Copyright (C) 2007  greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#


import doctest
import unittest

from rdkit import RDConfig
from rdkit.Chem.Fingerprints import DbFpSupplier
from rdkit.Dbase.DbConnection import DbConnect


def load_tests(loader, tests, ignore):
  """ Add the Doctests from the module """
  tests.addTests(doctest.DocTestSuite(DbFpSupplier, optionflags=doctest.ELLIPSIS))
  return tests


class TestCase(unittest.TestCase):

  def test_ForwardDbFpSupplier(self):
    # Additional tests to complete code coverage
    conn = DbConnect(RDConfig.RDTestDatabase, 'simple_combined')

    self.assertRaises(ValueError, DbFpSupplier.ForwardDbFpSupplier, conn.GetData(),
                      fpColName='typo')

    suppl = DbFpSupplier.ForwardDbFpSupplier(conn.GetData(), fpColName='AutoFragmentFp')
    self.assertIn('ID', suppl.GetColumnNames())
