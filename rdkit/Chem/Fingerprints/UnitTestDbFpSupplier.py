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
from __future__ import print_function

import doctest

from rdkit.Chem.Fingerprints import DbFpSupplier


def load_tests(loader, tests, ignore):  # pylint: disable=unused-argument
  """ Add the Doctests from the module """
  tests.addTests(doctest.DocTestSuite(DbFpSupplier, optionflags=doctest.ELLIPSIS))
  return tests
