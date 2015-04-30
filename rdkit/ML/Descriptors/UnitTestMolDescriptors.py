#
#  Copyright (C) 2002-2008  greg Landrum and Rational Discovery LLC
#

""" unit testing code for molecular descriptor calculators

"""
import unittest,os.path
import io
from rdkit.six.moves import cPickle
from rdkit import RDConfig
from rdkit.ML.Descriptors import MoleculeDescriptors
import numpy
from rdkit import Chem

class TestCase(unittest.TestCase):
  def setUp(self):
    self.descs = ['MolLogP','Chi1v']
    self.vers= ('1.1.0','1.0.0')
    self.calc = MoleculeDescriptors.MolecularDescriptorCalculator(self.descs)
    self.testD = [
      ('CCOC',    (0.6527, 1.40403)),
      ('CC=O',    (0.2052, 0.81305)),
      ('CCC(=O)O',(0.481, 1.48839))]

  def testGetNames(self):
    self.assertEqual(self.calc.GetDescriptorNames(),tuple(self.descs))
    
  def _testVals(self,calc,testD):
    for smi,vals in testD:
      mol = Chem.MolFromSmiles(smi)
      ans = numpy.array(vals)
      res = numpy.array(calc.CalcDescriptors(mol))
      self.assertTrue(max(abs(res-ans))<1e-4,'bad descriptor values for SMILES %s (%s)'%(smi,str(res)))
    
  def testCalcVals(self):
    self._testVals(self.calc,self.testD)

  def testSaveState(self):
    fName = os.path.join(RDConfig.RDCodeDir,'ML/Descriptors/test_data','molcalc.dsc')
    with open(fName,'r') as inTF:
      buf = inTF.read().replace('\r\n', '\n').encode('utf-8')
      inTF.close()
    with io.BytesIO(buf) as inF:
      calc = cPickle.load(inF)
    self.assertEqual(calc.GetDescriptorNames(),tuple(self.descs))
    self.assertEqual(calc.GetDescriptorVersions(),tuple(self.vers))
    self._testVals(calc,self.testD)
    
    
if __name__ == '__main__':
  unittest.main()
