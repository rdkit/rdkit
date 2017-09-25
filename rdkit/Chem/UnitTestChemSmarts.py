# $Id$
#
#  Copyright (C) 2003-2006  Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""basic unit testing code for the wrapper of the SMARTS matcher

"""
from rdkit import RDConfig
import unittest
import os.path
from rdkit import Chem


class TestCase(unittest.TestCase):

  def setUp(self):
    #print '\n%s: '%self.shortDescription(),
    fName = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', 'quinone.mol')
    self.m = Chem.MolFromMolFile(fName)
    assert self.m.GetNumAtoms() == 8, 'bad nAtoms'

  def testMatch(self):
    " testing smarts match "
    p = Chem.MolFromSmarts('CC(=O)C')
    matches = self.m.GetSubstructMatches(p)
    assert len(matches) == 2, 'bad UMapList: %s' % (str(res))
    for match in matches:
      assert len(match) == 4, 'bad match: %s' % (str(match))

  def testOrder(self):
    " testing atom order in smarts match "
    p = Chem.MolFromSmarts('CC(=O)C')
    matches = self.m.GetSubstructMatches(p)
    m = matches[0]
    atomList = [self.m.GetAtomWithIdx(x).GetSymbol() for x in m]
    assert atomList == ['C', 'C', 'O', 'C'], 'bad atom ordering: %s' % str(atomList)


if __name__ == '__main__':
  unittest.main()
