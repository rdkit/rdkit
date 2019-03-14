#
#  Copyright (C) 2001  greg Landrum
#
""" unit testing code for compound descriptors

"""


import unittest

from rdkit.ML.Descriptors import Parser
import math


class TestCase(unittest.TestCase):

  def setUp(self):
    self.piece1 = [['d1', 'd2', 's1'], ['d1', 'd2', 's1']]
    self.aDict = {'Fe': {'d1': 1, 'd2': 2, 's1': 'abc'}, 'Pt': {'d1': 10, 'd2': 20, 's1': 'def'}}
    self.pDict = {'d1': 100., 'd2': 200.}
    self.compos = [('Fe', 1), ('Pt', 1)]
    self.cExprs = ["SUM($1)", "SUM($1)+SUM($2)", "MEAN($1)", "DEV($2)", "MAX($1)", "MIN($2)",
                   "SUM($1)/$a", 'HAS($3,"def")', 'HAS($3,"xyz")', 'HAS($3)', "sqrt($a+($b+$b))"]
    self.results = [11., 33., 5.5, 9., 10., 2., 0.11, 1, 0, -666, math.sqrt(500)]
    self.tol = 0.0001

  def testSingleCalcs(self):
    " testing calculation of a single descriptor "
    for i in range(len(self.cExprs)):
      cExpr = self.cExprs[i]
      argVect = self.piece1 + [cExpr]
      res = Parser.CalcSingleCompoundDescriptor(self.compos, argVect, self.aDict, self.pDict)
      self.assertAlmostEqual(res, self.results[i], 2)

  def testMultipleCalcs(self):
    " testing calculation of multiple descriptors "
    for i in range(len(self.cExprs)):
      cExpr = self.cExprs[i]
      argVect = self.piece1 + [cExpr]
      res = Parser.CalcMultipleCompoundsDescriptor([self.compos, self.compos], argVect, self.aDict,
                                                   [self.pDict, self.pDict])
      self.assertAlmostEqual(res[0], self.results[i], delta=self.tol,
                             msg='Expression {0} failed'.format(cExpr))
      self.assertAlmostEqual(res[1], self.results[i], delta=self.tol,
                             msg='Expression {0} failed'.format(cExpr))

  def _test_exampleCode(self):
    Parser._exampleCode()


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
