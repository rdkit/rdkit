#
#  Copyright (C) 2001  greg Landrum
#

# unit testing code for the composite model COM server

from rdkit import RDConfig
import unittest
try:
  from win32com.client import Dispatch
except ImportError:
  Dispatch = None


@unittest.skipIf(Dispatch is None, 'Test for Windows only')
class TestCase(unittest.TestCase):

  def testConnect(self):
    # connecting to COM server
    ok = 1
    try:
      c = Dispatch('RD.Composite')
    except Exception:
      ok = 0
    assert ok and c is not None, 'connection to COM server failed'

  def testLoad(self):
    # loading a composite
    c = Dispatch('RD.Composite')
    ok = 1
    try:
      c.LoadComposite(RDConfig.RDCodeDir + '/ml/composite/test_data/composite_base.pkl')
    except Exception:
      ok = 0
    assert ok, 'LoadComposite failed'

  def testNames(self):
    # testing descriptor names
    c = Dispatch('RD.Composite')
    c.LoadComposite(RDConfig.RDCodeDir + '/ml/composite/test_data/composite_base.pkl')
    names = c.GetDescriptorNames()
    expectedNames = ('composition', 'max_atomic', 'has3d', 'has4d', 'has5d', 'elconc', 'atvol',
                     'isferro')
    assert names == expectedNames, 'GetDescriptorNames failed'

  def testInputOrder(self):
    # testing input order
    c = Dispatch('RD.Composite')
    c.LoadComposite(RDConfig.RDCodeDir + '/ml/composite/test_data/composite_base.pkl')
    names = c.GetDescriptorNames()
    ok = 1
    try:
      c.SetInputOrder(names)
    except Exception:
      ok = 0
    assert ok, 'SetInputOrder failed'

  def testClassify(self):
    # testing classification
    argV = ['CrPt3', 'fcc', 'AuCu3', 58.09549962, 36, 4, 0.228898, 2.219, 1, 3.67481803894, 1, 0, 1,
            0.619669341609, 14.523874905]
    nameV = ['composition', 'Structure', 'Structure_Type', 'Volume', 'Electrons_Per_Unit',
             'Atoms_Per_Unit', 'Hardness', 'DOS_Ef', 'isferro', 'max_atomic', 'has3d', 'has4d',
             'has5d', 'elconc', 'atvol']
    c = Dispatch('RD.Composite')
    c.LoadComposite(RDConfig.RDCodeDir + '/ml/composite/test_data/composite_base.pkl')
    c.SetInputOrder(nameV)
    res = c.ClassifyExample(argV)
    expected = [1, 1.0]
    assert res[0] == expected[0], 'bad prediction'
    assert res[1] == expected[1], 'bad confidence'


def TestSuite():  # pragma: nocover
  suite = unittest.TestSuite()
  suite.addTest(TestCase('testConnect'))
  suite.addTest(TestCase('testLoad'))
  suite.addTest(TestCase('testNames'))
  suite.addTest(TestCase('testInputOrder'))
  suite.addTest(TestCase('testClassify'))
  return suite


if __name__ == '__main__':  # pragma: nocover
  suite = TestSuite()
  unittest.TextTestRunner().run(suite)
