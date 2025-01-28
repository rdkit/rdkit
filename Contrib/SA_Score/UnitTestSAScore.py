import os.path
import unittest

import sascorer

from rdkit import Chem, RDConfig

print(sascorer.__file__)


class TestCase(unittest.TestCase):

  def test1(self):
    with open('data/zim.100.txt') as f:
      testData = [x.strip().split('\t') for x in f]
    testData.pop(0)
    for row in testData:
      smi = row[0]
      m = Chem.MolFromSmiles(smi)
      tgt = float(row[2])
      val = sascorer.calculateScore(m)
      self.assertAlmostEqual(tgt, val, 3)


if __name__ == '__main__':
  import getopt
  import re
  import sys
  doLong = 0
  if len(sys.argv) > 1:
    args, extras = getopt.getopt(sys.argv[1:], 'l')
    for arg, val in args:
      if arg == '-l':
        doLong = 1
      sys.argv.remove('-l')
  if doLong:
    for methName in dir(TestCase):
      if re.match('_test', methName):
        newName = re.sub('_test', 'test', methName)
        exec('TestCase.%s = TestCase.%s' % (newName, methName))

  unittest.main()
