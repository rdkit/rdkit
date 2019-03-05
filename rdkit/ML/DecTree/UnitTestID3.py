#
#  Copyright (C) 2000  greg Landrum
#
""" unit tests for the ID3 implementation """


import io
import unittest

from rdkit import RDConfig
from rdkit.ML.Data import MLData
from rdkit.ML.DecTree import ID3
import pickle


class ID3TestCase(unittest.TestCase):

  def setUp(self):
    self.basicTreeName = RDConfig.RDCodeDir + '/ML/DecTree/test_data/BasicTree.pkl'
    self.multiTreeName = RDConfig.RDCodeDir + '/ML/DecTree/test_data/MultiTree.pkl'

  def _setupBasicTree(self):
    examples = [[0, 0, 0, 0, 0], [0, 0, 0, 1, 0], [1, 0, 0, 0, 1], [2, 1, 0, 0, 1], [2, 2, 1, 0, 1],
                [2, 2, 1, 1, 0], [1, 2, 1, 1, 1], [0, 1, 0, 0, 0], [0, 2, 1, 0, 1], [2, 1, 1, 0, 1],
                [0, 1, 1, 1, 1], [1, 1, 0, 1, 1], [1, 0, 1, 0, 1], [2, 1, 0, 1, 0]]

    data = MLData.MLQuantDataSet(examples)
    attrs = list(range(0, data.GetNVars()))
    t1 = ID3.ID3Boot(data.GetAllData(), attrs, data.GetNPossibleVals())
    self.t1 = t1
    self.examples = examples

  def testBasicTree(self):
    # " testing basic tree growth "
    self._setupBasicTree()
    with open(self.basicTreeName, 'r') as inTFile:
      buf = inTFile.read().replace('\r\n', '\n').encode('utf-8')
      inTFile.close()
    with io.BytesIO(buf) as inFile:
      t2 = pickle.load(inFile)
    assert self.t1 == t2, 'Incorrect tree generated.'

  def _setupMultiTree(self):
    examples = [[0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 2], [0, 1, 1, 2], [1, 0, 0, 2], [1, 0, 1, 2],
                [1, 1, 0, 2], [1, 1, 1, 0]]

    data = MLData.MLQuantDataSet(examples)
    attrs = list(range(0, data.GetNVars()))
    t1 = ID3.ID3Boot(data.GetAllData(), attrs, data.GetNPossibleVals())
    self.t1 = t1
    self.examples = examples

  def testMultiTree(self):
    # " testing multivalued tree growth "
    self._setupMultiTree()
    with open(self.multiTreeName, 'r') as inTFile:
      buf = inTFile.read().replace('\r\n', '\n').encode('utf-8')
      inTFile.close()
    with io.BytesIO(buf) as inFile:
      t2 = pickle.load(inFile)
    assert self.t1 == t2, 'Incorrect tree generated.'

  def testClassify(self):
    # " testing basic tree classification "
    self._setupBasicTree()
    self.assertEqual(
      self.t1.ClassifyExample(self.examples[0]), self.examples[0][-1],
      'BasicExample 0 misclassified')
    self.assertEqual(
      self.t1.ClassifyExample(self.examples[1]), self.examples[1][-1],
      'BasicExample 1 misclassified')
    self.assertEqual(
      self.t1.ClassifyExample(self.examples[6]), self.examples[6][-1],
      'BasicExample 6 misclassified')
    self._setupMultiTree()
    self.assertEqual(
      self.t1.ClassifyExample(self.examples[0]), self.examples[0][-1],
      'MultiExample 0 misclassified')
    self.assertEqual(
      self.t1.ClassifyExample(self.examples[1]), self.examples[1][-1],
      'MultiExample 1 misclassified')
    self.assertEqual(
      self.t1.ClassifyExample(self.examples[6]), self.examples[6][-1],
      'MultiExample 6 misclassified')

  # ------------- force python in the ID3 code
  def _setupPyBasicTree(self):
    from rdkit.ML.InfoTheory import entropy
    ID3.entropy.InfoEntropy = entropy.PyInfoEntropy
    ID3.entropy.InfoGain = entropy.PyInfoGain

    examples = [[0, 0, 0, 0, 0], [0, 0, 0, 1, 0], [1, 0, 0, 0, 1], [2, 1, 0, 0, 1], [2, 2, 1, 0, 1],
                [2, 2, 1, 1, 0], [1, 2, 1, 1, 1], [0, 1, 0, 0, 0], [0, 2, 1, 0, 1], [2, 1, 1, 0, 1],
                [0, 1, 1, 1, 1], [1, 1, 0, 1, 1], [1, 0, 1, 0, 1], [2, 1, 0, 1, 0]]

    data = MLData.MLQuantDataSet(examples)
    attrs = list(range(0, data.GetNVars()))
    t1 = ID3.ID3Boot(data.GetAllData(), attrs, data.GetNPossibleVals())
    self.t1 = t1
    self.examples = examples

  def testPyBasicTree(self):
    # " testing basic tree growth (python entropy code) "
    self._setupPyBasicTree()
    with open(self.basicTreeName, 'r') as inTFile:
      buf = inTFile.read().replace('\r\n', '\n').encode('utf-8')
      inTFile.close()
    with io.BytesIO(buf) as inFile:
      t2 = pickle.load(inFile)
    assert self.t1 == t2, 'Incorrect tree generated.'

  def _setupPyMultiTree(self):
    from rdkit.ML.InfoTheory import entropy
    ID3.entropy.InfoEntropy = entropy.PyInfoEntropy
    ID3.entropy.InfoGain = entropy.PyInfoGain

    examples = [[0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 2], [0, 1, 1, 2], [1, 0, 0, 2], [1, 0, 1, 2],
                [1, 1, 0, 2], [1, 1, 1, 0]]

    data = MLData.MLQuantDataSet(examples)
    attrs = list(range(0, data.GetNVars()))
    t1 = ID3.ID3Boot(data.GetAllData(), attrs, data.GetNPossibleVals())
    self.t1 = t1
    self.examples = examples

  def testPyMultiTree(self):
    # " testing multivalued tree growth (python entropy code) "
    self._setupPyMultiTree()
    with open(self.multiTreeName, 'r') as inTFile:
      buf = inTFile.read().replace('\r\n', '\n').encode('utf-8')
      inTFile.close()
    with io.BytesIO(buf) as inFile:
      t2 = pickle.load(inFile)
    assert self.t1 == t2, 'Incorrect tree generated.'

  def testPyClassify(self):
    # " testing tree classification (python entropy code) "
    self._setupPyBasicTree()
    self.assertEqual(
      self.t1.ClassifyExample(self.examples[0]), self.examples[0][-1],
      'BasicExample 0 misclassified')
    self.assertEqual(
      self.t1.ClassifyExample(self.examples[1]), self.examples[1][-1],
      'BasicExample 1 misclassified')
    self.assertEqual(
      self.t1.ClassifyExample(self.examples[6]), self.examples[6][-1],
      'BasicExample 6 misclassified')
    self._setupMultiTree()
    self.assertEqual(
      self.t1.ClassifyExample(self.examples[0]), self.examples[0][-1],
      'MultiExample 0 misclassified')
    self.assertEqual(
      self.t1.ClassifyExample(self.examples[1]), self.examples[1][-1],
      'MultiExample 1 misclassified')
    self.assertEqual(
      self.t1.ClassifyExample(self.examples[6]), self.examples[6][-1],
      'MultiExample 6 misclassified')


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
