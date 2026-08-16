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

  def testGithub8251(self):
    # the smoothing above 8 used to take log(sascore - 8), which diverges at
    # the branch boundary and reported this molecule as 4.36, i.e. easier to
    # make than aspirin
    smi = ('C[N+]12CCCC1c1ccc[n+](c1)[Zn+2]21[O]C(=O)C(=O)[O]1.'
           'O=C1[OH+][Zn+2]2([OH+]C1=O)[OH+]C(=O)C(=O)[OH+]2')
    m = Chem.MolFromSmiles(smi)
    self.assertAlmostEqual(sascorer.calculateScore(m), 8.026, 3)

  def testEmptyMolecule(self):
    self.assertIsNone(sascorer.calculateScore(Chem.MolFromSmiles('')))

  def testFragmentTableSources(self):
    # the shipped binary table and the pickle it was generated from have to
    # give the same answers
    mols = [Chem.MolFromSmiles(smi) for smi in ('CC(=O)Oc1ccccc1C(=O)O',
                                                'C[C@H](N)C(=O)O',
                                                'C1CCC2(CC1)CCCCC2')]
    sascorer.readFragmentScores()
    fromBin = [sascorer.calculateScore(m) for m in mols]

    sascorer.readFragmentScores(
      os.path.join(RDConfig.RDContribDir, 'SA_Score', 'fpscores.pkl.gz'))
    fromPickle = [sascorer.calculateScore(m) for m in mols]

    sascorer.readFragmentScores()
    self.assertEqual(fromBin, fromPickle)


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
