# $Id$
#
#  Copyright (C) 2004-2005  Rational Discovery LLC
#   All Rights Reserved
#
import os
import unittest

from rdkit import RDConfig
from rdkit.DataStructs import ExplicitBitVect
from rdkit.ML.Data import DataUtils
from rdkit.ML.NaiveBayes import CrossValidate
from rdkit.ML.NaiveBayes.ClassificationModel import NaiveBayesClassifier


class TestCase(unittest.TestCase):

  def setUp(self):
    DataUtils.InitRandomNumbers((25, 25))

  def test1NaiveBayes(self):
    fName = os.path.join(RDConfig.RDCodeDir, 'ML', 'NaiveBayes', 'test_data', 'stddata.csv')
    data = DataUtils.TextFileToData(fName)

    examples = data.GetNamedData()

    nvars = data.GetNVars()
    attrs = list(range(1, nvars + 1))
    npvals = [0] + [3] * nvars + [2]
    qBounds = [0] + [2] * nvars + [0]
    mod, err = CrossValidate.CrossValidationDriver(examples, attrs, npvals, qBounds, silent=True)
    self.assertAlmostEqual(mod._classProbs[0], 0.5000, 4)
    self.assertAlmostEqual(mod._classProbs[1], 0.5000, 4)
    self.assertAlmostEqual(mod._QBoundVals[1][0], -0.0360, 4)
    self.assertAlmostEqual(mod._QBoundVals[1][1], 0.114)
    self.assertAlmostEqual(mod._QBoundVals[2][0], -0.7022, 4)
    self.assertAlmostEqual(mod._QBoundVals[2][1], -0.16635, 4)
    self.assertAlmostEqual(mod._QBoundVals[3][0], -0.3659, 4)
    self.assertAlmostEqual(mod._QBoundVals[3][1], 0.4305, 4)
    self.assertAlmostEqual(err, 0.2121, 4)

    mod, err = CrossValidate.CrossValidationDriver(examples, attrs, npvals, qBounds, silent=True,
                                                   calcTotalError=True)
    self.assertAlmostEqual(mod._classProbs[0], 0.515151, 4)
    self.assertAlmostEqual(mod._classProbs[1], 0.484848, 4)
    self.assertAlmostEqual(mod._QBoundVals[1][0], -0.40315, 4)
    self.assertAlmostEqual(mod._QBoundVals[1][1], 0.114)
    self.assertAlmostEqual(mod._QBoundVals[2][0], -0.62185, 4)
    self.assertAlmostEqual(mod._QBoundVals[2][1], -0.19965, 4)
    self.assertAlmostEqual(mod._QBoundVals[3][0], 0.4305, 4)
    self.assertAlmostEqual(mod._QBoundVals[3][1], 0.80305, 4)
    self.assertAlmostEqual(err, 0.14563, 4)

    mod, err = CrossValidate.CrossValidationDriver(examples, attrs, npvals, qBounds, silent=True,
                                                   replacementSelection=True)
    self.assertAlmostEqual(mod._classProbs[0], 0.5131578, 4)
    self.assertAlmostEqual(mod._classProbs[1], 0.4868421, 4)
    self.assertAlmostEqual(mod._QBoundVals[1][0], -0.036, 4)
    self.assertAlmostEqual(mod._QBoundVals[1][1], 0.93465, 4)
    self.assertAlmostEqual(mod._QBoundVals[2][0], -0.6696, 4)
    self.assertAlmostEqual(mod._QBoundVals[2][1], -0.19965, 4)
    self.assertAlmostEqual(mod._QBoundVals[3][0], -1.06785, 4)
    self.assertAlmostEqual(mod._QBoundVals[3][1], 0.4305, 4)
    self.assertAlmostEqual(err, 0.3, 4)

  def test2NaiveBayes(self):
    fName = os.path.join(RDConfig.RDCodeDir, 'ML', 'NaiveBayes', 'test_data', 'stddata.csv')
    data = DataUtils.TextFileToData(fName)
    examples = data.GetNamedData()

    nvars = data.GetNVars()
    attrs = list(range(1, nvars + 1))
    npvals = [0] + [3] * nvars + [2]
    qBounds = [0] + [2] * nvars + [0]
    mod, err = CrossValidate.CrossValidationDriver(examples, attrs, npvals, qBounds,
                                                   mEstimateVal=20.0, silent=True)
    self.assertTrue(isinstance(mod, NaiveBayesClassifier))
    self.assertAlmostEqual(err, 0.1818, 4)

    self.assertEqual(mod.GetName(), '')
    mod.SetName('modelName')
    self.assertEqual(mod.GetName(), 'modelName')
    mod.NameModel(None)
    self.assertEqual(mod.GetName(), 'NaiveBayesClassifier')

    self.assertGreater(len(mod.GetExamples()), 0)
    self.assertGreater(len(mod.GetTrainingExamples()), 0)
    self.assertEqual(sorted(mod.GetTrainingExamples() + mod.GetExamples()), sorted(examples))

  def test3(self):
    examples = [
      ['a', 1, 0, 1, 0, 1],
      ['b', 1, 0, 0, 0, 1],
      ['c', 1, 0, 1, 0, 0],
      ['d', 0, 1, 1, 0, 0],
      ['e', 0, 1, 1, 1, 0],
    ]

    nvars = len(examples[0]) - 2
    attrs = list(range(1, nvars + 1))
    npvals = [0] + [2] * nvars + [2]
    qBounds = [0] + [0] * nvars + [0]
    mdl = CrossValidate.makeNBClassificationModel(examples, attrs, npvals, qBounds)
    nWrong = 0
    for eg in examples:
      p = mdl.ClassifyExample(eg)
      if p != eg[-1]:
        nWrong += 1
    self.assertEqual(nWrong, 1)

    bitEx = []
    for eg in examples:
      newEg = [eg[0], None, eg[-1]]
      bv = ExplicitBitVect(nvars)
      for i in range(nvars):
        if eg[i + 1]:
          bv.SetBit(i)
      newEg[1] = bv
      bitEx.append(newEg)

    attrs = list(range(nvars))
    mdl2 = CrossValidate.makeNBClassificationModel(bitEx, attrs, npvals, qBounds, useSigs=True)
    nWrong = 0
    for eg in bitEx:
      p = mdl2.ClassifyExample(eg)
      if p != eg[-1]:
        nWrong += 1
    self.assertEqual(nWrong, 1)

    # now compare:
    for i in range(len(bitEx)):
      eg = examples[i]
      p1 = mdl.ClassifyExample(eg)
      bitEg = bitEx[i]
      p2 = mdl2.ClassifyExample(bitEg)
      self.assertEqual(p1, p2)
      v1 = mdl.GetClassificationDetails()
      v2 = mdl.GetClassificationDetails()
      self.assertAlmostEqual(v1, v2, 4)

  def test4(self):
    examples = [
      ['a', 1, 0, 1, 0, 1],
      ['b', 1, 0, 0, 0, 1],
      ['c', 1, 0, 1, 0, 0],
      ['d', 0, 1, 1, 0, 0],
      ['e', 0, 1, 1, 1, 0],
    ]

    nvars = len(examples[0]) - 2
    origNVars = nvars
    nvars = 10
    npvals = [0] + [2] * nvars + [2]
    qBounds = [0] + [0] * nvars + [0]

    bitEx = []
    for eg in examples:
      newEg = [eg[0], None, eg[-1]]
      bv = ExplicitBitVect(nvars)
      for i in range(origNVars):
        if eg[i + 1]:
          bv.SetBit(i)

      # this bit will yield perfect accuracy if
      #  the attrs argument isn't being used properly:
      if eg[-1]:
        bv.SetBit(origNVars)
      newEg[1] = bv
      bitEx.append(newEg)

    attrs = list(range(origNVars))
    mdl2 = CrossValidate.makeNBClassificationModel(bitEx, attrs, npvals, qBounds, useSigs=True)
    nWrong = 0
    for eg in bitEx:
      p = mdl2.ClassifyExample(eg)
      if p != eg[-1]:
        nWrong += 1
    self.assertEqual(nWrong, 1)

  def _test5(self):  # disabled because CMIM was removed # pragma: nocover
    examples = [
      ['a', 1, 0, 1, 0, 1, 1, 0, 1],
      ['b', 1, 0, 0, 0, 1, 0, 0, 1],
      ['c', 1, 0, 1, 0, 1, 1, 0, 0],
      ['d', 0, 1, 1, 0, 1, 0, 0, 0],
      ['e', 0, 1, 1, 1, 0, 1, 0, 0],
    ]

    nvars = len(examples[0]) - 2
    npvals = [0] + [2] * nvars + [2]
    qBounds = [0] + [0] * nvars + [0]

    bitEx = []
    for eg in examples:
      newEg = [eg[0], None, eg[-1]]
      bv = ExplicitBitVect(nvars)
      for i in range(nvars):
        if eg[i + 1]:
          bv.SetBit(i)

      # this bit will yield perfect accuracy if
      #  the attrs argument isn't being used properly:
      newEg[1] = bv
      bitEx.append(newEg)

    attrs = list(range(nvars))
    mdl2 = CrossValidate.makeNBClassificationModel(bitEx, attrs, npvals, qBounds, useSigs=True,
                                                   useCMIM=2)
    nWrong = 0
    for eg in bitEx:
      p = mdl2.ClassifyExample(eg)
      if p != eg[-1]:
        nWrong += 1
    self.assertEqual(nWrong, 1)


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
