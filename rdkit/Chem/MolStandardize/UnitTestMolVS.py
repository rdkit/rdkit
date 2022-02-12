#!/usr/env/bin python
# -*- coding: utf-8 -*-
"""
unit testing for MolVS steps for PubChem Substances
tests include
molvs.standardize_smiles
Standardizer().normalize
Standardizer().disconnect_metals
Standardizer().reionize
molvs.standardize.canonicalize_tautomer_smiles
molvs.validate.Validator()
Standardizer().fragment_parent
"""
import gzip
import os.path
import unittest
from collections import namedtuple

import molvs
from molvs import Standardizer, validate
from rdkit import Chem, RDConfig

doLong = False
TestData = namedtuple('TestData', 'lineNo,smiles,mol,expected')

class TestCase(unittest.TestCase):
    dataPCS_standardize_smiles100k = os.path.join(RDConfig.RDBaseDir, 'rdkit', 'Chem', 'MolStandardize', 'test_data', '100kPCS_standardize_sm.csv.gz')
    dataPCS_standardize_smiles1k = os.path.join(RDConfig.RDBaseDir,'rdkit', 'Chem', 'MolStandardize', 'test_data', '1kPCS_standardize_sm.csv.gz')
    dataPCS_nomralized1k = os.path.join(RDConfig.RDBaseDir,'rdkit', 'Chem', 'MolStandardize', 'test_data', '1kPCS_normalized.csv.gz')
    dataPCS_nomralized100k = os.path.join(RDConfig.RDBaseDir, 'rdkit', 'Chem', 'MolStandardize', 'test_data', '100kPCS_normalized.csv.gz')
    dataPCS_metal100k = os.path.join(RDConfig.RDBaseDir,'rdkit', 'Chem', 'MolStandardize', 'test_data', '100kPCS_metals.csv.gz')
    dataPCS_metal1k = os.path.join(RDConfig.RDBaseDir,'rdkit', 'Chem', 'MolStandardize', 'test_data', '1kPCS_metals.csv.gz')
    dataPCS_reionize100k = os.path.join(RDConfig.RDBaseDir,'rdkit', 'Chem', 'MolStandardize', 'test_data', '100kPCS_reionize.csv.gz')
    dataPCS_reionize1k = os.path.join(RDConfig.RDBaseDir,'rdkit', 'Chem', 'MolStandardize', 'test_data', '1kPCS_reionize.csv.gz')
    dataPCS_fragment100k = os.path.join(RDConfig.RDBaseDir,'rdkit', 'Chem', 'MolStandardize', 'test_data', '100kPCS_fragment.csv.gz')
    dataPCS_fragmnet1k = os.path.join(RDConfig.RDBaseDir,'rdkit', 'Chem', 'MolStandardize', 'test_data', '1kPCS_fragment.csv.gz')
    dataPCS_tautomer100k = os.path.join(RDConfig.RDBaseDir,'rdkit', 'Chem', 'MolStandardize', 'test_data', '100kPCS_tautomer.csv.gz')
    dataPCS_tautomer1k = os.path.join(RDConfig.RDBaseDir,'rdkit', 'Chem', 'MolStandardize', 'test_data', '1kPCS_tautomer.csv.gz')
    dataPCS_validate100k = os.path.join(RDConfig.RDBaseDir,'rdkit', 'Chem', 'MolStandardize', 'test_data', '100kPCS_validate.csv.gz')
    dataPCS_validate1k = os.path.join(RDConfig.RDBaseDir,'rdkit', 'Chem', 'MolStandardize', 'test_data', '1kPCS_validate.csv.gz')

    @classmethod
    def readPCSdata(cls, datafile):
        for data in readPCStestData(datafile):
            yield data

    def testStandardizeSmShort(self):
        for data in self.readPCSdata(self.dataPCS_standardize_smiles1k):
            ss = molvs.standardize_smiles(data.smiles)
            self.assertEqual(ss, data.expected)

    def testStandardizeSmLong(self):
        if not doLong:
            raise unittest.SkipTest('long test')
        for data in self.readPCSdata(self.dataPCS_standardize_smiles100k):
            try:
                ss = molvs.standardize_smiles(data.smiles)
            except Exception:
                raise AssertionError(f'Line {data.lineNo}: MolVS standardization failed for SMILES {data.smiles}')
            self.assertEqual(ss, data.expected)

    def testNormalizeShort(self):
        for data in self.readPCSdata(self.dataPCS_nomralized1k):
            n = Standardizer()
            nm = n.normalize(data.mol)
            ns = Chem.MolToSmiles(nm)
            self.assertEqual(ns, data.expected)

    def testNormalizeLong(self):
        if not doLong:
            raise unittest.SkipTest('long test')
        for data in self.readPCSdata(self.dataPCS_nomralized100k):
            try:
                n = Standardizer()
                nm = n.normalize(data.mol)
                ns = Chem.MolToSmiles(nm)
            except Exception:
                raise AssertionError(f'Line {data.lineNo}: MolVS normalization failed for SMILES {data.smiles}')
            self.assertEqual(ns, data.expected)

    def testMetalShort(self):
        for data in self.readPCSdata(self.dataPCS_metal1k):
            n = Standardizer()
            nm = n.disconnect_metals(data.mol)
            ns = Chem.MolToSmiles(nm)
            self.assertEqual(ns, data.expected)

    def testMetalLong(self):
        if not doLong:
            raise unittest.SkipTest('long test')
        for data in self.readPCSdata(self.dataPCS_metal100k):
            try:
                n = Standardizer()
                nm = n.disconnect_metals(data.mol)
                ns = Chem.MolToSmiles(nm)
            except Exception:
                raise AssertionError(f'Line {data.lineNo}: MolVS normalization failed for SMILES {data.smiles}')
            self.assertEqual(ns, data.expected)

    def testReionizeShort(self):
        for data in self.readPCSdata(self.dataPCS_reionize1k):
            n = Standardizer()
            nm = n.reionize(data.mol)
            ns = Chem.MolToSmiles(nm)
            self.assertEqual(ns, data.expected)

    def testReionizeLong(self):
        if not doLong:
            raise unittest.SkipTest('long test')
        for data in self.readPCSdata(self.dataPCS_reionize100k):
            try:
                n = Standardizer()
                nm = n.reionize(data.mol)
                ns = Chem.MolToSmiles(nm)
            except Exception:
                raise AssertionError(f'Line {data.lineNo}: MolVS normalization failed for SMILES {data.smiles}')
            self.assertEqual(ns, data.expected)

    def testFragmentShort(self):
        for data in self.readPCSdata(self.dataPCS_fragmnet1k):
            s = Standardizer()
            frag = s.fragment_parent(data.mol)
            ns = Chem.MolToSmiles(frag)
            self.assertEqual(ns, data.expected)

    def testFragmentLong(self):
        if not doLong:
            raise unittest.SkipTest('long test')
        for data in self.readPCSdata(self.dataPCS_fragment100k):
            try:
                s = Standardizer()
                frag = s.fragment_parent(data.mol)
                ns = Chem.MolToSmiles(frag)
            except Exception:
                raise AssertionError(f'Line {data.lineNo}: MolVS normalization failed for SMILES {data.smiles}')
            self.assertEqual(ns, data.expected)

    def testTautomerShort(self):
        for data in self.readPCSdata(self.dataPCS_tautomer1k):
            canon_taut = molvs.standardize.canonicalize_tautomer_smiles(data.smiles)
            self.assertEqual(canon_taut, data.expected)

    def testTautomerLong(self):
        if not doLong:
            raise unittest.SkipTest('long test')
        for data in self.readPCSdata(self.dataPCS_tautomer100k):
            try:
                canon_taut = molvs.standardize.canonicalize_tautomer_smiles(data.smiles)
            except Exception:
                raise AssertionError(f'Line {data.lineNo}: MolVS normalization failed for SMILES {data.smiles}')
            self.assertEqual(canon_taut, data.expected)

    def testValidateShort(self):
        for data in self.readPCSdata(self.dataPCS_validate1k):
            v = validate.Validator()
            stdout = v(data.mol)
            if stdout:
                stdout = ' '.join(stdout)
            else:
                stdout = str(stdout)
            self.assertEqual(stdout, data.expected)

    def testValidateLong(self):
        if not doLong:
            raise unittest.SkipTest('long test')
        for i, data in enumerate(self.readPCSdata(self.dataPCS_validate100k)):
            try:
                v = validate.Validator()
                stdout = v(data.mol)
                if stdout:
                    stdout = ','.join(stdout)
                else:
                    stdout = str(stdout)
                # print(i, stdout)
            except Exception:
                raise AssertionError(
                    f'Line {data.lineNo}: MolVS normalization failed for SMILES {data.smiles}')
            self.assertEqual(stdout, data.expected)

def readPCStestData(filename):
    """ Read test data for MolVS from file """
    with gzip.open(filename, 'rt') as f:
        for lineNo, line in enumerate(f, 1):
            if line[0] == '#':
                continue
            smiles = line.strip().split(',')[0]
            expected = ",".join(line.strip().split(',')[1:])
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                raise AssertionError(f'molecule construction failed on line {lineNo}')
            yield TestData(lineNo, smiles, mol, expected)

if __name__ == '__main__':
    import argparse
    import sys

    parser = argparse.ArgumentParser()
    parser.add_argument('-l', default=False, action='store_true', dest='doLong')
    args = parser.parse_args()
    doLong = args.doLong

    # Remove the -l flag if present so that it isn't interpreted by unittest.main()
    if '-l' in sys.argv:
        sys.argv.remove('-l')
    unittest.main()
