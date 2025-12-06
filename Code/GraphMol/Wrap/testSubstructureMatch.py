#
#  Copyright (C) 2025 Greg Landrum
#         All Rights Reserved
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

#
import unittest

from rdkit import Chem


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def testExtraAtomMatch(self):
    """Test setting extra atom matching functions """
    mol = Chem.MolFromSmiles('CCCC')
    mol.GetAtomWithIdx(1).SetIntProp('foo', 1)
    mol.GetAtomWithIdx(2).SetIntProp('foo', 2)
    patt = Chem.MolFromSmiles('CC')
    patt.GetAtomWithIdx(0).SetIntProp('bar', 2)
    patt.GetAtomWithIdx(1).SetIntProp('bar', 1)

    def atomCheck(qatom, ratom):
      if not qatom.HasProp('bar') or not ratom.HasProp('foo'):
        return False
      return qatom.GetIntProp('bar') == ratom.GetIntProp('foo')

    matches = mol.GetSubstructMatches(patt)
    self.assertEqual(len(matches), 3)

    params = Chem.SubstructMatchParameters()
    params.setExtraAtomCheckFunc(atomCheck)

    matches = mol.GetSubstructMatches(patt, params)
    self.assertEqual(len(matches), 1)
    self.assertEqual(matches[0], (2, 1))

  def testExtraBondMatch(self):
    """Test setting extra bond matching functions """
    mol = Chem.MolFromSmiles('CCCC')
    mol.GetBondWithIdx(1).SetIntProp('foo', 1)
    patt = Chem.MolFromSmiles('CC')
    patt.GetBondWithIdx(0).SetIntProp('bar', 1)

    def bondCheck(qbond, rbond):
      if not qbond.HasProp('bar') or not rbond.HasProp('foo'):
        return False
      return qbond.GetIntProp('bar') == rbond.GetIntProp('foo')

    matches = mol.GetSubstructMatches(patt)
    self.assertEqual(len(matches), 3)

    params = Chem.SubstructMatchParameters()
    params.setExtraBondCheckFunc(bondCheck)

    matches = mol.GetSubstructMatches(patt, params)
    self.assertEqual(len(matches), 1)
    self.assertEqual(matches[0], (1, 2))

  def testCoordsFunctor(self):
    m = Chem.MolFromSmiles('CCCC |(0,0,0;1,0,0;2,0,0;3,0,0)|')
    patt = Chem.MolFromSmiles('CC |(3,0,0.1;2,0,0)|')

    params = Chem.SubstructMatchParameters()
    coordMatcher = Chem.AtomCoordsMatcher(tol=.11)
    params.setExtraAtomCheckFunc(coordMatcher)
    matches = m.GetSubstructMatches(patt, params)
    self.assertEqual(len(matches), 1)
    self.assertEqual(matches[0], (3, 2))
    
    # lifetime management 1:
    coordMatcher = None
    matches = m.GetSubstructMatches(patt, params)
    self.assertEqual(len(matches), 1)
    self.assertEqual(matches[0], (3, 2))
    
    # lifetime management 2:
    params.setExtraAtomCheckFunc(Chem.AtomCoordsMatcher(tol=.11))
    matches = m.GetSubstructMatches(patt, params)
    self.assertEqual(len(matches), 1)
    self.assertEqual(matches[0], (3, 2))


if __name__ == '__main__':
  unittest.main()
