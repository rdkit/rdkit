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

  def testBondRingCountQueryQueryMatches(self):
    " testing query-query matching with not queries on bond ring counts "
    params = Chem.SubstructMatchParameters()
    params.useQueryQueryMatches = True

    # --- @{0} vs @{0} ---
    target = Chem.MolFromSmarts('*@{0}*')
    query = Chem.MolFromSmarts('*@{0}*')
    Chem.GetSymmSSSR(target)
    Chem.GetSymmSSSR(query)
    matches = target.GetSubstructMatches(query, params=params)
    assert len(matches) == 1, 'positive @{0} should match positive @{0}'

    # --- !@{0} vs @{0} (should NOT match) ---
    target = Chem.MolFromSmarts('*!@{0}*')
    query = Chem.MolFromSmarts('*@{0}*')
    Chem.GetSymmSSSR(target)
    Chem.GetSymmSSSR(query)
    matches = target.GetSubstructMatches(query, params=params)
    assert len(matches) == 0, 'not @{0} should NOT match positive @{0}'

    # --- !@{0} vs !@{0} (should match) ---
    target = Chem.MolFromSmarts('*!@{0}*')
    query = Chem.MolFromSmarts('*!@{0}*')
    Chem.GetSymmSSSR(target)
    Chem.GetSymmSSSR(query)
    matches = target.GetSubstructMatches(query, params=params)
    assert len(matches) == 1, 'not @{0} should match not @{0}'

    # --- @{0} vs !@{0} (should NOT match) ---
    target = Chem.MolFromSmarts('*@{0}*')
    query = Chem.MolFromSmarts('*!@{0}*')
    Chem.GetSymmSSSR(target)
    Chem.GetSymmSSSR(query)
    matches = target.GetSubstructMatches(query, params=params)
    assert len(matches) == 0, 'positive @{0} should NOT match not @{0}'

    # --- @{2} vs @{2} ---
    target = Chem.MolFromSmarts('*@{2}*')
    query = Chem.MolFromSmarts('*@{2}*')
    Chem.GetSymmSSSR(target)
    Chem.GetSymmSSSR(query)
    matches = target.GetSubstructMatches(query, params=params)
    assert len(matches) == 1, 'positive @{2} should match positive @{2}'

    # --- !@{2} vs @{2} (should NOT match) ---
    target = Chem.MolFromSmarts('*!@{2}*')
    query = Chem.MolFromSmarts('*@{2}*')
    Chem.GetSymmSSSR(target)
    Chem.GetSymmSSSR(query)
    matches = target.GetSubstructMatches(query, params=params)
    assert len(matches) == 0, 'not @{2} should NOT match positive @{2}'

    # --- !@{2} vs !@{2} (should match) ---
    target = Chem.MolFromSmarts('*!@{2}*')
    query = Chem.MolFromSmarts('*!@{2}*')
    Chem.GetSymmSSSR(target)
    Chem.GetSymmSSSR(query)
    matches = target.GetSubstructMatches(query, params=params)
    assert len(matches) == 1, 'not @{2} should match not @{2}'

    # --- @{2} vs !@{2} (should NOT match) ---
    target = Chem.MolFromSmarts('*@{2}*')
    query = Chem.MolFromSmarts('*!@{2}*')
    Chem.GetSymmSSSR(target)
    Chem.GetSymmSSSR(query)
    matches = target.GetSubstructMatches(query, params=params)
    assert len(matches) == 0, 'positive @{2} should NOT match not @{2}'

    # --- @{1} vs @{2} (should NOT match, different ring counts) ---
    target = Chem.MolFromSmarts('*@{1}*')
    query = Chem.MolFromSmarts('*@{2}*')
    Chem.GetSymmSSSR(target)
    Chem.GetSymmSSSR(query)
    matches = target.GetSubstructMatches(query, params=params)
    assert len(matches) == 0, 'positive @{1} should NOT match positive @{2}'

    # --- !@{1} vs !@{2} (should NOT match, different negated ring counts) ---
    target = Chem.MolFromSmarts('*!@{1}*')
    query = Chem.MolFromSmarts('*!@{2}*')
    Chem.GetSymmSSSR(target)
    Chem.GetSymmSSSR(query)
    matches = target.GetSubstructMatches(query, params=params)
    assert len(matches) == 0, 'not @{1} should NOT match not @{2}'

    # --- @{1-2} vs @{1-2} (should match, same range) ---
    target = Chem.MolFromSmarts('*@{1-2}*')
    query = Chem.MolFromSmarts('*@{1-2}*')
    Chem.GetSymmSSSR(target)
    Chem.GetSymmSSSR(query)
    matches = target.GetSubstructMatches(query, params=params)
    assert len(matches) == 1, 'range @{1-2} should match range @{1-2}'

    # --- @{1-2} vs @{3-4} (should NOT match, non-overlapping ranges) ---
    target = Chem.MolFromSmarts('*@{1-2}*')
    query = Chem.MolFromSmarts('*@{3-4}*')
    Chem.GetSymmSSSR(target)
    Chem.GetSymmSSSR(query)
    matches = target.GetSubstructMatches(query, params=params)
    assert len(matches) == 0, 'range @{1-2} should NOT match range @{3-4}'

    # --- range vs exact: both directions should NOT match ---
    # queriesMatch correctly returns false for range-vs-exact (different query types)
    # and there is no fallback, so no match in either direction
    target = Chem.MolFromSmarts('*@{1-3}*')
    query = Chem.MolFromSmarts('*@{2}*')
    Chem.GetSymmSSSR(target)
    Chem.GetSymmSSSR(query)
    matches = target.GetSubstructMatches(query, params=params)
    assert len(matches) == 0, 'range @{1-3} should NOT match exact @{2}'

    target = Chem.MolFromSmarts('*@{2}*')
    query = Chem.MolFromSmarts('*@{1-3}*')
    Chem.GetSymmSSSR(target)
    Chem.GetSymmSSSR(query)
    matches = target.GetSubstructMatches(query, params=params)
    assert len(matches) == 0, 'exact @{2} should NOT match range @{1-3}'

if __name__ == '__main__':
  unittest.main()
