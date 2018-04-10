# $Id$
#
#  Copyright (C) 2006  greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import doctest
import unittest

from rdkit.Chem.FeatMaps import FeatMapPoint


def load_tests(loader, tests, ignore):
  """ Add the Doctests from the module """
  tests.addTests(doctest.DocTestSuite(FeatMapPoint, optionflags=doctest.ELLIPSIS))
  return tests


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
