#
#  Copyright (C) 2026 Greg Landrum and other RDKit contributors
#         All Rights Reserved
#
import os
import unittest

from rdkit import Chem


def lists_to_tuple(lst):
  return tuple(tuple(x) for x in lst)


def ats_to_tuple(lst):
  return tuple(x.GetIdx() for x in lst)


class TestCase(unittest.TestCase):

  def test_recursiveatoms(self):
    m = Chem.MolFromSmiles('CC=C')
    sma = '[C;!$(C=C)]'
    p = Chem.MolFromSmarts(sma)
    tpl = lists_to_tuple(m.GetSubstructMatches(p))
    self.assertEqual(tpl, ((0, ), ))

    qat = Chem.AtomFromSmarts(sma)
    self.assertEqual(ats_to_tuple(m.GetAtomsMatchingQuery(qat)), (0, ))


if __name__ == '__main__':
  import os
  if "RDTESTCASE" in os.environ:
    suite = unittest.TestSuite()
    testcases = os.environ["RDTESTCASE"]
    for name in testcases.split(':'):
      suite.addTest(TestCase(name))

    runner = unittest.TextTestRunner()
    runner.run(suite)
  else:
    unittest.main()
