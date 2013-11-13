#!/usr/bin/env python

import rdkit
import rdkit.Chem.rdShape
#print rdkit.Chem.rdShape.Align.__doc__
from rdkit import Chem
from rdkit import RDConfig
import unittest, os,time

class TestCase(unittest.TestCase) :
    def test_shape_reference(self):
        #mol= Chem.MolFromMolFile(os.path.join(RDConfig.RDBaseDir,'Code/GraphMol/Depictor','test_data/7UPJ_xtal.mol'))
        f1 = os.path.join(RDConfig.RDBaseDir,'Code','Shape','test_data','reference.sdf')
        f2 = os.path.join(RDConfig.RDBaseDir,'Code','Shape','test_data','database.sdf')
        d1 = os.path.join(RDConfig.RDBaseDir,'Code','Shape','test_data','output_rdkit_1.tab')
        #d2 = os.path.join(RDConfig.RDBaseDir,'Code','Shape','test_data','output_rdkit_2.tab')


        mol1 = Chem.SDMolSupplier(f1).next()
        supp = Chem.SDMolSupplier(f1)
        mol2 = supp.next()
        score = rdkit.Chem.rdShape.Align(mol1, mol2, maxIter=10)
        self.failUnlessAlmostEqual(score,1.0,1)
        
        supp = Chem.SDMolSupplier(f2)
        t1=time.time()
        tanimoto = [rdkit.Chem.rdShape.Align(mol1, m, maxIter=0) for m in supp]
        t2=time.time()
        print "Alignment time: %.2f"%(t2-t1)
        
        df1 = [x.strip().split('\t') for x in file(d1).readlines()]
        taniCol = df1[0].index('Shape-it::Tanimoto')
        for i,tani in enumerate(tanimoto):
            self.failUnlessAlmostEqual(float(df1[i+1][taniCol]),tani,3)

if __name__ == '__main__':
  unittest.main()
