#
#  Copyright (C) 2002-2008  greg Landrum and Rational Discovery LLC
#
""" unit testing code for molecular descriptor calculators

"""
import os.path
import pickle
import unittest
from io import BytesIO, StringIO

import numpy

from rdkit import Chem, RDConfig
from rdkit.ML.Descriptors import Descriptors, MoleculeDescriptors
from rdkit.TestRunner import redirect_stdout


class TestCase(unittest.TestCase):

    def setUp(self):
        self.descs = ['MolLogP', 'Chi1v']
        self.vers = ('1.1.0', '1.0.0')
        self.calc = MoleculeDescriptors.MolecularDescriptorCalculator(self.descs)
        self.testD = [('CCOC', (0.6527, 1.40403)), ('CC=O', (0.2052, 0.81305)), ('CCC(=O)O',
                                                                                 (0.481, 1.48839))]

    def testGetNames(self):
        self.assertEqual(self.calc.GetDescriptorNames(), tuple(self.descs))

    def _testVals(self, calc, testD):
        for smi, vals in testD:
            mol = Chem.MolFromSmiles(smi)
            ans = numpy.array(vals)
            res = numpy.array(calc.CalcDescriptors(mol))
            self.assertTrue(
              max(abs(res - ans)) < 1e-4, 'bad descriptor values for SMILES %s (%s)' % (smi, str(res)))

    def testCalcVals(self):
        self._testVals(self.calc, self.testD)

    def testSaveState(self):
        fName = os.path.join(RDConfig.RDCodeDir, 'ML/Descriptors/test_data', 'molcalc.dsc')
        with open(fName, 'r') as inTF:
            buf = inTF.read().replace('\r\n', '\n').encode('utf-8')
            inTF.close()
        inF = BytesIO(buf)
        calc = pickle.load(inF)
        self.assertEqual(calc.GetDescriptorNames(), tuple(self.descs))
        self.assertEqual(calc.GetDescriptorVersions(), tuple(self.vers))
        self._testVals(calc, self.testD)

        f = StringIO()
        with redirect_stdout(f):
            calc.ShowDescriptors()
        s = f.getvalue()
        for name in calc.GetDescriptorNames():
            self.assertIn(name, s)

        self.assertIn('Wildman-Crippen LogP value', calc.GetDescriptorSummaries())
        self.assertIn('N/A', calc.GetDescriptorSummaries())

        funcs = calc.GetDescriptorFuncs()
        self.assertEqual(len(funcs), len(self.descs))
        for f in funcs:
            self.assertTrue(callable(f))


class TestDescriptors(unittest.TestCase):

    def test_DescriptorCalculator(self):
        calc = Descriptors.DescriptorCalculator()
        self.assertRaises(NotImplementedError, calc.ShowDescriptors)
        self.assertRaises(NotImplementedError, calc.GetDescriptorNames)
        self.assertRaises(NotImplementedError, calc.CalcDescriptors, None)

        calc.simpleList = ['simple1', 'simple2']
        calc.compoundList = ['cmpd1', 'cmpd2']
        f = StringIO()
        with redirect_stdout(f):
            calc.ShowDescriptors()
        s = f.getvalue()
        for name in calc.simpleList:
            self.assertIn(name, s)
        for name in calc.compoundList:
            self.assertIn(name, s)

    def test_github3511(self):
        mol = Chem.MolFromSmiles('C')
        descriptors = [name for name, _ in Chem.Descriptors.descList]
        calculator = MoleculeDescriptors.MolecularDescriptorCalculator(descriptors)
        calculator.CalcDescriptors(mol)

        # This should not raise a pickling exception
        pickle.dumps(mol)




if __name__ == '__main__':  # pragma: nocover
    unittest.main()
