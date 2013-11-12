#!/usr/bin/env python

import rdkit
import rdkit.Chem.rdShape
#print rdkit.Chem.rdShape.Align.__doc__
from rdkit import Chem
from rdkit import RDConfig
import unittest, os

class TestCase(unittest.TestCase) :
    def test_shape_reference(self):
        #mol= Chem.MolFromMolFile(os.path.join(RDConfig.RDBaseDir,'Code/GraphMol/Depictor','test_data/7UPJ_xtal.mol'))
        f1 = os.path.join(RDConfig.RDBaseDir,'Code','Shape','test_data','reference.sdf')
        f2 = os.path.join(RDConfig.RDBaseDir,'Code','Shape','test_data','database.sdf')
        d1 = os.path.join(RDConfig.RDBaseDir,'Code','Shape','test_data','output_rdkit_1.tab')
        d2 = os.path.join(RDConfig.RDBaseDir,'Code','Shape','test_data','output_rdkit_2.tab')


        mol1 = Chem.SDMolSupplier(f1).next()
        supp = Chem.SDMolSupplier(f1)
        mol2 = supp.next()
        score = rdkit.Chem.rdShape.Align(mol1, mol2, maxIter=10)
        print score
        
        supp = Chem.SDMolSupplier(f2)
        tanimoto = [rdkit.Chem.rdShape.Align(mol1, m, maxIter=0) for m in supp]
        
        
        import pandas as pd
        df1 = pd.read_csv(d1, sep="\t")
        df2 = pd.read_csv(d2, sep="\t")
        df = pd.DataFrame(tanimoto, columns=["Shape-it::Tanimoto",])
        
        import numpy as np
        np.testing.assert_array_almost_equal(df["Shape-it::Tanimoto"], df1["Shape-it::Tanimoto"], decimal=3)
        np.testing.assert_array_almost_equal(df["Shape-it::Tanimoto"], df2["Shape-it::Tanimoto"], decimal=3)

if __name__ == '__main__':
  unittest.main()
