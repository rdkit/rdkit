# $Id$
#
#  Copyright (C) 2004-2006  Greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""unit testing code for the AnalyzeComposite functionality

"""
import os
import unittest

from rdkit import RDConfig
from rdkit.ML import AnalyzeComposite
from rdkit.six.moves import cPickle as pickle  # @UnresolvedImport


class TestCase(unittest.TestCase):

  def setUp(self):
    self.baseDir = os.path.join(RDConfig.RDCodeDir, 'ML', 'test_data')

  def test1_Issue163(self):
    name1 = os.path.join(self.baseDir, 'humanoral.1.pkl')
    try:
      with open(name1, 'rb') as pklF:
        c1 = pickle.load(pklF)
    except Exception:  # pragma: nocover
      c1 = None
    self.assertTrue(c1)
    name2 = os.path.join(self.baseDir, 'humanoral.2.pkl')
    try:
      with open(name2, 'rb') as pklF:
        c2 = pickle.load(pklF)
    except Exception:  # pragma: nocover
      c2 = None
    self.assertTrue(c2)

    try:
      res = AnalyzeComposite.ProcessIt([c1, c2], verbose=-1)
    except Exception:  # pragma: nocover
      import traceback
      traceback.print_exc()
      ok = 0
    else:
      ok = 1
    self.assertTrue(ok)

    self.assertTrue(res[0][0] == 'BALABANJ')
    self.assertTrue(res[1][0] == 'BERTZCT')
    self.assertTrue(res[-1][0] == 'FR_ALLYLIC_OXID')
    for entry in res:
      self.assertTrue(len(entry) == 5)


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
