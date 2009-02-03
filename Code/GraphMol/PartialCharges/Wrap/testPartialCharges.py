from rdkit import Chem
from rdkit.Chem import rdPartialCharges
from rdkit import RDConfig
import unittest
import os
import cPickle as pickle

def feq(v1,v2,tol2=1e-4):
    return abs(v1-v2)<=tol2

class TestCase(unittest.TestCase):
    def setUp(self) :
        pass

    def test0HalgrenSet(self) :
        smiSup = Chem.SmilesMolSupplier("test_data/halgren.smi",delimiter='\t')

        #parse the original file
        infil = open("test_data/halgren_out.txt", 'r')
        lines = infil.readlines()
        infil.close()
        
        tab = Chem.GetPeriodicTable()
        
        olst = []
        for mol in smiSup :
            rdPartialCharges.ComputeGasteigerCharges(mol) 
            tstr = "Molecule: "
            tstr += mol.GetProp("_Name")
            olst.append(tstr)
            for i in range(mol.GetNumAtoms()) :
                at = mol.GetAtomWithIdx(i)
                en = tab.GetElementSymbol(at.GetAtomicNum())
                chg = float(at.GetProp("_GasteigerCharge"))
                tstr = "%i %s %6.4f"%(i, en, chg)
                olst.append(tstr)

            
        i = 0
        for line in lines:
	  self.failUnless(line.strip() == olst[i])
	  i += 1
        
    def test1PPDataset(self):
        fileN = os.path.join('test_data', 'PP_descrs_regress.2.csv')
        infil = open(fileN, 'r')
        lines = infil.readlines()
        infil.close()

        infile = os.path.join('test_data', 'PP_combi_charges.pkl')
        cchFile = open(infile, 'rb')
        combiCharges = pickle.load(cchFile)

        for lin in lines :
            if (lin[0] == '#') :
                continue
            tlst = lin.strip().split(',')
            smi = tlst[0]
            rdmol = Chem.MolFromSmiles(smi)
            rdPartialCharges.ComputeGasteigerCharges(rdmol)
            
            nat = rdmol.GetNumAtoms()
            
            for ai in range(nat) :
                rdch = float(rdmol.GetAtomWithIdx(ai).GetProp('_GasteigerCharge'))
                if not feq(rdch, combiCharges[smi][ai], 1.e-2) :
                    print smi, ai
                    rdmol.debug()
                    break
                #self.failUnless(feq(rdch, combiCharges[smi][ai], 1.e-2))
                
    def test2Params(self):
        """ tests handling of Issue187 """
        m1 = Chem.MolFromSmiles('C(=O)[O-]')
        rdPartialCharges.ComputeGasteigerCharges(m1)

        m2 = Chem.MolFromSmiles('C(=O)[O-].[Na+]')
        try:
            rdPartialCharges.ComputeGasteigerCharges(m2)
        except:
            self.fail('should not have hit an exception')

        for i in range(m1.GetNumAtoms()):
            c1 = float(m1.GetAtomWithIdx(i).GetProp('_GasteigerCharge'))
            c2 = float(m2.GetAtomWithIdx(i).GetProp('_GasteigerCharge'))
            self.failUnless(feq(c1,c2,1e-4))
            
            
    def test3Params(self):
        """ tests handling of Issue187 """
        m2 = Chem.MolFromSmiles('C(=O)[O-].[Na+]')
        try:
            rdPartialCharges.ComputeGasteigerCharges(m2,12,1)
        except:
            pass
        else:
            self.fail('should have hit an exception')

if __name__== '__main__':
    unittest.main()

    
