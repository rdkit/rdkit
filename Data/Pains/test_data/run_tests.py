#
# Copyright (C) 2015 Greg Landrum
#   This file is part of the RDKit.
#   The contents are covered by the terms of the BSD license
#   which is included in the file license.txt, found at the root
#   of the RDKit source tree.
#

import csv
import os
import unittest

from rdkit import Chem, RDConfig


class TestCase(unittest.TestCase):

  def setUp(self):
    self.basePath = os.path.join(RDConfig.RDDataDir, 'Pains')
    self.painsFile = os.path.join(self.basePath, 'wehi_pains.csv')
    with open(self.painsFile, 'r') as inf:
      self.painsDefs = [x for x in csv.reader(inf)]
    self.matchers = [Chem.MolFromSmarts(x[0], mergeHs=True) for x in self.painsDefs]

  def test1(self):
    " molecules that we know should match "
    with open(os.path.join(self.basePath, 'test_data', 'test_set3.txt'), 'r') as inf:
      testData = [x.strip().split() for x in inf if x[0] != '#']
    for line in testData:
      self.assertEqual(len(line), 5)
      id_ = int(line[0])
      m = Chem.MolFromSmiles(line[2])
      self.assertTrue(m is not None)
      self.assertTrue(m.HasSubstructMatch(self.matchers[id_]))
      self.assertTrue(Chem.AddHs(m).HasSubstructMatch(self.matchers[id_]))


if __name__ == '__main__':
  unittest.main()
