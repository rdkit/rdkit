#
#  Copyright (C) 2001  greg Landrum
#
""" unit testing code for the descriptor COM server

"""


import unittest

from rdkit import RDConfig
import numpy as np

try:
  from win32com.client import Dispatch
except ImportError:
  Dispatch = None


@unittest.skipIf(Dispatch is None, 'Windows test')
class TestCase(unittest.TestCase):

  def setUp(self):
    print('\n%s: ' % self.shortDescription(), end='')

  def testConnectToCOMServer(self):
    # " testing connection "
    Dispatch('RD.DescCalc')

  def testLoadCalculator(self):
    # " testing load "
    c = Dispatch('RD.DescCalc')
    c.LoadCalculator(RDConfig.RDCodeDir + '/ml/descriptors/test_data/ferro.dsc')

  def testNames(self):
    # " testing GetDescriptorNames "
    c = Dispatch('RD.DescCalc')
    c.LoadCalculator(RDConfig.RDCodeDir + '/ml/descriptors/test_data/ferro.dsc')
    names = c.GetDescriptorNames()
    expectedNames = ('MAX_DED', 'has3d', 'has4d', 'has5d', 'elconc', 'atvol')
    assert names == expectedNames, 'GetDescriptorNames failed (%s != %s)' % (repr(names),
                                                                             repr(expectedNames))

  def testCalc(self):
    # " testing descriptor calculation "
    argV = ['CrPt3', 'fcc', 'AuCu3', 58.09549962, 1, 4, 0.228898, 8.876, 1]
    nameV = ['Compound', 'Structure', 'Structure_Type', 'Volume', 'Z', 'Atoms_per_Formula_Unit',
             'Hardness', 'RawDOS_Ef', 'IsFerromagnetic']
    c = Dispatch('RD.DescCalc')
    c.LoadCalculator(RDConfig.RDCodeDir + '/ml/descriptors/test_data/ferro.dsc')
    descVect = np.array(c.CalcDescriptors(argV, nameV))
    expected = np.array((3.67481803894, 1, 0, 1, 0.619669341609, 14.523874905))
    diffV = abs(descVect - expected)
    assert max(diffV) < 0.0001, 'bad descriptors: %s, %s' % (str(expected), str(descVect))


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
