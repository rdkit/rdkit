#
#  Copyright (C) 2000  greg Landrum
#
""" unit tests for the Neural network trainer implementation

   this basically works out **all** of the network code

"""


import unittest

from rdkit.ML.Neural.ActFuncs import Sigmoid, TanH
from rdkit.ML.Neural.NetNode import NetNode
from rdkit.ML.Neural.Network import Network


class TestCaseActFuncs(unittest.TestCase):

  def test_Sigmoid(self):
    f = Sigmoid()
    self.assertAlmostEqual(f(0), 0.5)
    self.assertAlmostEqual(f(0), f.Eval(0))
    self.assertAlmostEqual(f.Deriv(0), 0.25)
    self.assertAlmostEqual(f(1), 1.0 - f(-1))
    self.assertAlmostEqual(f(2), 1.0 - f(-2))
    self.assertAlmostEqual(f.Deriv(1), f.Deriv(-1))
    self.assertAlmostEqual(f.Deriv(2), f.Deriv(-2))
    self.assertLess(f(1), f(2))
    self.assertLess(f.Deriv(2), f.Deriv(1))
    self.assertAlmostEqual(f.Deriv(1), f.DerivFromVal(f(1)))

  def test_TanH(self):
    f = TanH()
    self.assertAlmostEqual(f(0), 0.0)
    self.assertAlmostEqual(f(0), f.Eval(0))
    self.assertAlmostEqual(f.Deriv(0), 1.0)
    self.assertAlmostEqual(f(1), -f(-1))
    self.assertAlmostEqual(f(2), -f(-2))
    self.assertAlmostEqual(f.Deriv(1), f.Deriv(-1))
    self.assertAlmostEqual(f.Deriv(2), f.Deriv(-2))
    self.assertLess(f(1), f(2))
    self.assertLess(f.Deriv(2), f.Deriv(1))
    self.assertAlmostEqual(f.Deriv(1), f.DerivFromVal(f(1)))


class TestCaseNetNode(unittest.TestCase):

  def test_NetNode(self):
    # A node without input always returns 1
    nodeList = [None] * 2
    node = NetNode(0, nodeList)
    nodeList[0] = node
    valVect = [None] * 2
    self.assertEqual(node.Eval(valVect), 1)
    self.assertEqual(valVect, [1, None])

    node = NetNode(1, nodeList, inputNodes=[0], weights=[0.1])
    self.assertRaises(AssertionError, node.SetWeights, [0, 1])
    self.assertRaises(AssertionError, node.SetInputs, [0, 1])


class TestCaseNetwork(unittest.TestCase):

  def test_Network(self):
    nodeCounts = [2, 2, 1, 2]
    net = Network(nodeCounts)
    self.assertEqual(net.GetNumNodes(), 7)
    self.assertEqual(len(net.GetAllNodes()), 7)
    self.assertEqual(net.GetInputNodeList(), [0, 1])
    self.assertEqual(net.GetHiddenLayerNodeList(0), [2, 3])
    self.assertEqual(net.GetHiddenLayerNodeList(1), [4])
    self.assertEqual(net.GetOutputNodeList(), [5, 6])

    # We get a representation of the network
    s = str(net)
    self.assertIn('Network', s)


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
