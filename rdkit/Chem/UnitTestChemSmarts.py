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
import os.path
import unittest

from rdkit import Chem, RDConfig


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

  def testBondRingCount(self):
    " testing @{n} bond ring count SMARTS "
    # parsing
    for sma in ['*@{2}*', '*@{2-}*', '*@{-2}*', '*@{1-3}*', '*@{0}*']:
      q = Chem.MolFromSmarts(sma)
      assert q is not None, 'failed to parse %s' % sma

    # naphthalene: one fusion bond in 2 rings
    m = Chem.MolFromSmiles('c1ccc2ccccc2c1')
    q = Chem.MolFromSmarts('*@{2}*')
    matches = m.GetSubstructMatches(q)
    assert len(matches) == 1, 'naphthalene @{2} match count'

    q2 = Chem.MolFromSmarts('*@{2-}*')
    matches2 = m.GetSubstructMatches(q2)
    assert len(matches2) == 1, 'naphthalene @{2-} match count'

    q3 = Chem.MolFromSmarts('*@{1}*')
    matches3 = m.GetSubstructMatches(q3)
    assert len(matches3) == 10, 'naphthalene @{1} match count'

    # biphenylene: two fusion bonds
    bip = Chem.MolFromSmiles('c1ccc2c(c1)c1ccccc12')
    q = Chem.MolFromSmarts('*@{2-}*')
    matches = bip.GetSubstructMatches(q)
    assert len(matches) == 2, 'biphenylene @{2-} match count'

    # cubane: all 12 bonds in 2 rings
    cub = Chem.MolFromSmiles('C12C3C4C1C5C2C3C45')
    q = Chem.MolFromSmarts('*@{2}*')
    matches = cub.GetSubstructMatches(q)
    assert len(matches) == 12, 'cubane @{2} match count'

    # cyclohexane: all bonds in 1 ring
    cyc = Chem.MolFromSmiles('C1CCCCC1')
    q = Chem.MolFromSmarts('*@{1}*')
    matches = cyc.GetSubstructMatches(q)
    assert len(matches) == 6, 'cyclohexane @{1} match count'

    # ethane: bond in 0 rings
    eth = Chem.MolFromSmiles('CC')
    q = Chem.MolFromSmarts('*@{0}*')
    matches = eth.GetSubstructMatches(q)
    assert len(matches) == 1, 'ethane @{0} match count'

    # round-trip SMARTS
    for sma in ['*@{2}*', '*@{2-}*', '*@{-2}*', '*@{1-3}*']:
      q = Chem.MolFromSmarts(sma)
      out = Chem.MolToSmarts(q)
      q2 = Chem.MolFromSmarts(out)
      assert Chem.MolToSmarts(q2) == out, 'round-trip failed for %s' % sma

    # negation
    m = Chem.MolFromSmiles('c1ccc2ccccc2c1')
    q = Chem.MolFromSmarts('*!@{2}*')
    matches = m.GetSubstructMatches(q)
    assert len(matches) == 10, 'negation match count'

if __name__ == '__main__':
  unittest.main()
