
import unittest
import os
import io

import pickle

from rdkit import Chem
from rdkit.Chem import rdPartialCharges
from rdkit import RDConfig


def feq(v1, v2, tol2=1e-4):
    return abs(v1 - v2) <= tol2


class TestCase(unittest.TestCase):

    def setUp(self):
        pass

    def test0HalgrenSet(self):
        smiSup = Chem.SmilesMolSupplier(
          os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'PartialCharges', 'Wrap', 'test_data',
                       'halgren.smi'), delimiter='\t')

        # parse the original file
        with open(
            os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'PartialCharges', 'Wrap', 'test_data',
                         'halgren_out.txt'), 'r') as infil:
            lines = infil.readlines()

        tab = Chem.GetPeriodicTable()

        olst = []
        for mol in smiSup:
            rdPartialCharges.ComputeGasteigerCharges(mol)
            tstr = "Molecule: "
            tstr += mol.GetProp("_Name")
            olst.append(tstr)
            for i in range(mol.GetNumAtoms()):
                at = mol.GetAtomWithIdx(i)
                en = tab.GetElementSymbol(at.GetAtomicNum())
                chg = float(at.GetProp("_GasteigerCharge"))
                tstr = "%i %s %6.4f" % (i, en, chg)
                olst.append(tstr)

        i = 0
        for line in lines:
            self.assertTrue(line.strip() == olst[i])
            i += 1

    def test1PPDataset(self):
        fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'PartialCharges', 'Wrap',
                             'test_data', 'PP_descrs_regress.2.csv')
        infil = open(fileN, 'r')
        lines = infil.readlines()
        infil.close()

        infile = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'PartialCharges', 'Wrap',
                              'test_data', 'PP_combi_charges.pkl')
        with open(infile, 'r') as cchtFile:
            buf = cchtFile.read().replace('\r\n', '\n').encode('utf-8')
            cchtFile.close()
        with io.BytesIO(buf) as cchFile:
            combiCharges = pickle.load(cchFile)

        for lin in lines:
            if (lin[0] == '#'):
                continue
            tlst = lin.strip().split(',')
            smi = tlst[0]
            rdmol = Chem.MolFromSmiles(smi)
            rdPartialCharges.ComputeGasteigerCharges(rdmol)

            nat = rdmol.GetNumAtoms()
            failed = False
            for ai in range(nat):
                rdch = float(rdmol.GetAtomWithIdx(ai).GetProp('_GasteigerCharge'))
                if not feq(rdch, combiCharges[smi][ai], 1.e-2):
                    failed = True
                    print(smi, ai, rdch, combiCharges[smi][ai])
            if failed:
                rdmol.Debug()
            self.assertFalse(failed)

    def test2Params(self):
        """ tests handling of Issue187 """
        m1 = Chem.MolFromSmiles('C(=O)[O-]')
        rdPartialCharges.ComputeGasteigerCharges(m1)

        m2 = Chem.MolFromSmiles('C(=O)[O-].[Na+]')
        rdPartialCharges.ComputeGasteigerCharges(m2)

        for i in range(m1.GetNumAtoms()):
            c1 = float(m1.GetAtomWithIdx(i).GetProp('_GasteigerCharge'))
            c2 = float(m2.GetAtomWithIdx(i).GetProp('_GasteigerCharge'))
            self.assertTrue(feq(c1, c2, 1e-4))

    def test3Params(self):
        """ tests handling of Issue187 """
        m2 = Chem.MolFromSmiles('C(=O)[O-].[Na+]')
        with self.assertRaisesRegex(Exception, ""):
            rdPartialCharges.ComputeGasteigerCharges(m2, 12, 1)

    def testGithubIssue20(self):
        """ tests handling of Github issue 20 """
        m1 = Chem.MolFromSmiles('CB(O)O')
        rdPartialCharges.ComputeGasteigerCharges(m1)
        chgs = [-0.030, 0.448, -0.427, -0.427]
        for i in range(m1.GetNumAtoms()):
            c1 = float(m1.GetAtomWithIdx(i).GetProp('_GasteigerCharge'))
            self.assertAlmostEqual(c1, chgs[i], 3)

    def testGithubIssue577(self):
        """ tests handling of Github issue 577 """
        m1 = Chem.MolFromSmiles('CCO')
        from locale import setlocale, LC_NUMERIC
        try:
            setlocale(LC_NUMERIC, "de_DE")
        except Exception:
            # can't set the required locale, might as well just return
            return
        try:
            rdPartialCharges.ComputeGasteigerCharges(m1)
            for at in m1.GetAtoms():
                float(at.GetProp('_GasteigerCharge'))
        finally:
            setlocale(LC_NUMERIC, "C")
        rdPartialCharges.ComputeGasteigerCharges(m1)
        for at in m1.GetAtoms():
            float(at.GetProp('_GasteigerCharge'))
    def testGithub2480(self):
        with self.assertRaisesRegex(Exception, "^Python argument types"):
            rdPartialCharges.ComputeGasteigerCharges(None)


if __name__ == '__main__':
    unittest.main()
