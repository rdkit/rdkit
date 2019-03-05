#
#  Copyright (C) 2002-2008  greg Landrum and Rational Discovery LLC
#
""" unit tests for the model and descriptor packager """
import os
import random
import unittest
from xml.dom import minidom
from xml.etree import ElementTree as ET

from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import Descriptors
from rdkit.ML.Composite import Composite
from rdkit.ML.Data import DataUtils
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator
from rdkit.ML.ModelPackage import Packager, PackageUtils
from rdkit.ML.ModelPackage.Packager import ModelPackage
from io import BytesIO
import pickle


def feq(a, b, tol=1e-4):
    return abs(a - b) <= tol


class TestCase(unittest.TestCase):

    def setUp(self):
        self.dataDir = os.path.join(RDConfig.RDCodeDir, 'ML/ModelPackage/test_data')
        self.testD = [
          # NOTE: the confidences here can be twitchy due to changes in descriptors:
          ('Fc1ccc(NC(=O)c2cccnc2Oc3cccc(c3)C(F)(F)F)c(F)c1', 0, 0.8),
          # (r'CN/1(=C\C=C(/C=C1)\C\2=C\C=N(C)(Cl)\C=C2)Cl',0,0.70),
          (r'NS(=O)(=O)c1cc(ccc1Cl)C2(O)NC(=O)c3ccccc32', 1, 0.70),
        ]

    def _loadPackage(self):
        with open(os.path.join(self.dataDir, 'Jan9_build3_pkg.pkl'), 'r') as pkgTF:
            buf = pkgTF.read().replace('\r\n', '\n').encode('utf-8')
            pkgTF.close()
        io = BytesIO(buf)
        pkg = pickle.load(io)
        return pkg

    def _verify(self, pkg, testD):
        for smi, pred, conf in testD:
            m = Chem.MolFromSmiles(smi)
            self.assertTrue(m is not None, 'SMILES: %s failed\n' % (smi))
            p, c = pkg.Classify(m)
            assert p == pred, 'bad prediction (%d) for smiles %s' % (p, smi)
            assert feq(c, conf), 'bad confidence (%f) for smiles %s' % (c, smi)

    def _verify2(self, pkg, testD):
        for smi, pred, conf in testD:
            m = Chem.MolFromSmiles(smi)
            self.assertTrue(m is not None, 'SMILES: %s failed\n' % (smi))
            p, c = pkg.Classify(m)
            assert p == pred, 'bad prediction (%d) for smiles %s' % (p, smi)
            assert feq(c, conf), 'bad confidence (%f) for smiles %s' % (c, smi)
            p, c = pkg.Classify(m)
            assert p == pred, 'bad prediction (%d) for smiles %s' % (p, smi)
            assert feq(c, conf), 'bad confidence (%f) for smiles %s' % (c, smi)

    def testBuild(self):
        # """ tests building and screening a packager """
        with open(os.path.join(self.dataDir, 'Jan9_build3_calc.dsc'), 'r') as calcTF:
            buf = calcTF.read().replace('\r\n', '\n').encode('utf-8')
            calcTF.close()
        calc = pickle.load(BytesIO(buf))
        with open(os.path.join(self.dataDir, 'Jan9_build3_model.pkl'), 'rb') as modelF:
            model = pickle.load(modelF)
        pkg = Packager.ModelPackage(descCalc=calc, model=model)
        self._verify(pkg, self.testD)

    def testLoad(self):
        # """ tests loading and screening a packager """
        pkg = self._loadPackage()
        self._verify(pkg, self.testD)

    def testLoad2(self):
        # """ tests loading and screening a packager 2 """
        pkg = self._loadPackage()
        self._verify2(pkg, self.testD)

    def testPerm1(self):
        # """ tests the descriptor remapping stuff in a packager """
        pkg = self._loadPackage()
        calc = pkg.GetCalculator()
        names = calc.GetDescriptorNames()
        ref = {}
        DataUtils.InitRandomNumbers((23, 42))
        for smi, _, _ in self.testD:
            for desc in names:
                fn = getattr(Descriptors, desc, lambda x: 777)
                m = Chem.MolFromSmiles(smi)
                ref[desc] = fn(m)

            for _ in range(5):
                perm = list(names)
                random.shuffle(perm, random=random.random)

                m = Chem.MolFromSmiles(smi)
                for desc in perm:
                    fn = getattr(Descriptors, desc, lambda x: 777)
                    val = fn(m)
                    assert feq(val, ref[desc], 1e-4), '%s: %s(%s): %f!=%f' % (str(perm), smi, desc, val,
                                                                              ref[desc])

    def testPerm2(self):
        # """ tests the descriptor remapping stuff in a packager """
        pkg = self._loadPackage()
        calc = pkg.GetCalculator()
        names = calc.GetDescriptorNames()
        DataUtils.InitRandomNumbers((23, 42))
        perm = list(names)
        random.shuffle(perm, random=random.random)
        calc.simpleList = perm
        calc.descriptorNames = perm
        pkg.Init()
        self._verify(pkg, self.testD)

    def test_ModelPackage(self):
        pkg = self._loadPackage()

        self.assertTrue(isinstance(pkg.GetCalculator(), MolecularDescriptorCalculator))
        pkg.SetCalculator('calculator')
        self.assertEqual(pkg.GetCalculator(), 'calculator')

        self.assertTrue(isinstance(pkg.GetModel(), Composite.Composite))
        pkg.SetModel('model')
        self.assertEqual(pkg.GetModel(), 'model')

        self.assertEqual(pkg.GetDataset(), None)
        pkg.SetDataset('dataset')
        self.assertEqual(pkg.GetDataset(), 'dataset')

        self.assertEqual(pkg.GetNotes(), 'General purpose model built from PhysProp data')
        pkg.SetNotes('notes')
        self.assertEqual(pkg.GetNotes(), 'notes')

        # Here seems to be a difference between Python 2 and 3. The next assert works in Python 3,
        # but fails in Python 2
        # self.assertFalse(hasattr(pkg, '_supplementalData'))
        self.assertEqual(pkg.GetSupplementalData(), [])
        self.assertTrue(hasattr(pkg, '_supplementalData'))

        delattr(pkg, '_supplementalData')
        pkg.AddSupplementalData('supp1')
        self.assertTrue(hasattr(pkg, '_supplementalData'))
        self.assertEqual(pkg.GetSupplementalData(), ['supp1'])
        pkg.AddSupplementalData('supp2')
        self.assertEqual(pkg.GetSupplementalData(), ['supp1', 'supp2'])

        pkg = ModelPackage()
        self.assertFalse(pkg._initialized)
        pkg.Init()
        self.assertFalse(pkg._initialized)

    def test_PackageUtils(self):
        pkg = self._loadPackage()
        xml = PackageUtils.PackageToXml(
          pkg, dataPerformance=[('label', ['accuracy', 'avgCorrect', 'avgIncorrect']), ],
          recommendedThreshold=0.2, classDescriptions=[('a', 'texta'), ('b', 'textb')],
          modelType='model type', modelOrganism='model organism')
        s = prettyXML(xml.getroot())
        self.assertIn('<RDModelInfo>', s)


def prettyXML(xml):
    s = ET.tostring(xml, encoding='utf-8')
    tree = minidom.parseString(s)
    return tree.toprettyxml(indent=' ')


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
