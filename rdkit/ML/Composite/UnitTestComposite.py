# $Id$
#
#  Copyright (C) 2001-2008  greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" unit testing code for composite models

"""
import unittest
import io
import pickle
from rdkit.ML.Composite import Composite
from rdkit.ML.DecTree.DecTree import DecTreeNode as Node
from rdkit import RDConfig


class TestCase(unittest.TestCase):

  def setUp(self):
    with open(RDConfig.RDCodeDir + '/ML/Composite/test_data/ferro.pkl', 'r') as pklTF:
      buf = pklTF.read().replace('\r\n', '\n').encode('utf-8')
      pklTF.close()
    with io.BytesIO(buf) as pklF:
      self.examples = pickle.load(pklF)
    self.varNames = ['composition', 'max_atomic', 'has3d', 'has4d', 'has5d', 'elconc', 'atvol',
                     'isferro']
    self.qBounds = [[], [1.89, 3.53], [], [], [], [0.55, 0.73], [11.81, 14.52], []]
    self.nPoss = [0, 3, 2, 2, 2, 3, 3, 2]
    self.attrs = list(range(1, len(self.varNames) - 1))
    from rdkit.ML.Data import DataUtils
    DataUtils.InitRandomNumbers((23, 43))

  def testQuantize(self):
    # testing data quantization
    qBounds = [[], [1, 2, 3]]
    examples = [['foo', 0], ['foo', 1.5], ['foo', 5.5], ['foo', 2.5]]
    answers = [['foo', 0], ['foo', 1], ['foo', 3], ['foo', 2]]
    nPoss = [0, 4]
    composite = Composite.Composite()
    composite.SetQuantBounds(qBounds, nPoss)
    for i in range(len(examples)):
      qEx = composite.QuantizeExample(examples[i])
      self.assertEqual(qEx, answers[i])

  def testTreeGrow(self):
    # testing tree-based composite
    with open(RDConfig.RDCodeDir + '/ML/Composite/test_data/composite_base.pkl', 'r') as pklTF:
      buf = pklTF.read().replace('\r\n', '\n').encode('utf-8')
      pklTF.close()
    with io.BytesIO(buf) as pklF:
      self.refCompos = pickle.load(pklF)

    composite = Composite.Composite()
    composite._varNames = self.varNames
    composite.SetQuantBounds(self.qBounds, self.nPoss)
    from rdkit.ML.DecTree import CrossValidate
    driver = CrossValidate.CrossValidationDriver
    pruner = None
    composite.Grow(self.examples, self.attrs, [], buildDriver=driver, pruner=pruner, nTries=100,
                   silent=1)
    composite.AverageErrors()
    composite.SortModels(sortOnError=False)
    self.assertEqual(composite.countList, sorted(composite.countList))
    self.assertNotEqual(composite.errList, sorted(composite.errList))
    composite.SortModels()
    self.assertNotEqual(composite.countList, sorted(composite.countList))
    self.assertEqual(composite.errList, sorted(composite.errList))

    # with open(RDConfig.RDCodeDir+'/ML/Composite/test_data/composite_base.pkl','wb') as pklF:
    #   pickle.dump(composite,pklF)

    self.treeComposite = composite
    self.assertEqual(len(composite), len(self.refCompos))
    for i in range(len(composite)):
      t1, c1, e1 = composite[i]
      t2, c2, e2 = self.refCompos[i]
      self.assertEqual(e1, e2)
      # we used to check for equality here, but since there are redundant errors,
      #  that's non-trivial.
      # assert t1 == t2, 'tree mismatch'
      # assert c1 == c2, 'count mismatch'
    s = str(composite)
    self.assertIn('Composite', s)
    self.assertIn('Model', s)
    self.assertIn('error', s)

  def testErrorEstimate(self):
    # testing out-of-bag error estimates

    compos = Composite.Composite()
    compos.SetQuantBounds([(0.5, ), (0.5, ), (0.5, ), []], [2, 2, 2, 2])
    compos.SetDescriptorNames(('D0', 'D1', 'D2', 'Act'))
    compos.SetInputOrder(('ID', 'D0', 'D1', 'D2', 'Act'))
    data = [['A', 0, 0, 0, 0], ['B', 1, 0, 0, 1], ['C', 0, 1, 0, 0], ['D', 1, 1, 1, 1]]

    #
    #  Build and validate three simple trees:
    #
    t1 = Node(None, 'D0', 0)
    n = Node(t1, 'D1', 1)
    t1.AddChildNode(n)
    n.AddChildNode(Node(n, '0', 0, isTerminal=1))
    n.AddChildNode(Node(n, '1', 1, isTerminal=1))
    n = Node(t1, 'D2', 2)
    t1.AddChildNode(n)
    n.AddChildNode(Node(n, '1', 1, isTerminal=1))
    n.AddChildNode(Node(n, '0', 0, isTerminal=1))
    assert t1.ClassifyExample(data[0][1:]) == 0
    assert t1.ClassifyExample(data[1][1:]) == 1
    assert t1.ClassifyExample(data[2][1:]) == 1
    assert t1.ClassifyExample(data[3][1:]) == 0
    t1._trainIndices = (0, 1)
    compos.AddModel(t1, .5)

    t2 = Node(None, 'D1', 1)
    n = Node(t2, 'D0', 0)
    t2.AddChildNode(n)
    n.AddChildNode(Node(n, '0', 0, isTerminal=1))
    n.AddChildNode(Node(n, '1', 1, isTerminal=1))
    n = Node(t2, 'D2', 2)
    t2.AddChildNode(n)
    n.AddChildNode(Node(n, '0', 0, isTerminal=1))
    n.AddChildNode(Node(n, '1', 1, isTerminal=1))
    assert t2.ClassifyExample(data[0][1:]) == 0
    assert t2.ClassifyExample(data[1][1:]) == 1
    assert t2.ClassifyExample(data[2][1:]) == 0
    assert t2.ClassifyExample(data[3][1:]) == 1
    t2._trainIndices = (1, 2)
    compos.AddModel(t2, 0.0)

    t3 = Node(None, 'D0', 0)
    n = Node(t3, 'D2', 2)
    t3.AddChildNode(n)
    n.AddChildNode(Node(n, '0', 0, isTerminal=1))
    n.AddChildNode(Node(n, '1', 1, isTerminal=1))
    n = Node(t3, 'D1', 1)
    t3.AddChildNode(n)
    n.AddChildNode(Node(n, '0', 0, isTerminal=1))
    n.AddChildNode(Node(n, '1', 1, isTerminal=1))
    assert t3.ClassifyExample(data[0][1:]) == 0
    assert t3.ClassifyExample(data[1][1:]) == 0
    assert t3.ClassifyExample(data[2][1:]) == 0
    assert t3.ClassifyExample(data[3][1:]) == 1
    t3._trainIndices = (2, 3)
    compos.AddModel(t3, 0.25)

    #
    #  validate the composite itself:
    #
    pred, conf = compos.ClassifyExample(data[0])
    assert pred == 0
    assert conf == 1.0
    pred, conf = compos.ClassifyExample(data[1])
    assert pred == 1
    assert conf == 2. / 3.
    pred, conf = compos.ClassifyExample(data[2])
    assert pred == 0
    assert conf == 2. / 3.
    pred, conf = compos.ClassifyExample(data[3])
    assert pred == 1
    assert conf == 2. / 3.

    self.assertEqual(compos.GetVoteDetails(), [0, 1, 1])
    self.assertEqual(compos.GetInputOrder(), [1, 2, 3, 4])

    #
    #  now test the out-of-bag calculation:
    #
    pred, conf = compos.ClassifyExample(data[0], onlyModels=(1, 2))
    assert pred == 0
    assert conf == 1.0
    pred, conf = compos.ClassifyExample(data[1], onlyModels=(2, ))
    assert pred == 0
    assert conf == 1.0
    pred, conf = compos.ClassifyExample(data[2], onlyModels=(0, ))
    assert pred == 1
    assert conf == 1.0
    pred, conf = compos.ClassifyExample(data[3], onlyModels=(0, 1))
    assert pred == 0
    assert conf == 0.5

    compos.ClearModelExamples()

  def test_exceptions(self):
    compos = Composite.Composite()
    compos.SetQuantBounds([(0.5, ), (0.5, ), (0.5, ), []], [2, 2, 2, 2])
    compos.SetDescriptorNames(('ID', 'D0', 'D1', 'D2', 'Act'))
    compos.SetInputOrder(('ID', 'D2', 'D1', 'D0', 'Act'))
    self.assertEqual(compos._mapOrder, [0, 3, 2, 1, 4])
    # Probes caught exception for ID
    compos.SetInputOrder(('D2', 'D1', 'D0', 'Act'))
    self.assertEqual(compos._mapOrder, [0, 2, 1, 0, 3])
    # Probes caught exception for Act
    compos.SetInputOrder(('ID', 'D2', 'D1', 'D0'))
    self.assertEqual(compos._mapOrder, [0, 3, 2, 1, -1])

    self.assertRaises(ValueError, compos.SetInputOrder, ('Unknown', 'D0'))


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
