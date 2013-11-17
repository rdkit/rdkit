#!/usr/bin/env python

import rdkit
import rdkit.Chem.rdShape
#print rdkit.Chem.rdShape.Align.__doc__
from rdkit import Chem
from rdkit import RDConfig
import unittest, os,time

class TestCase(unittest.TestCase) :
    def test_shape_reference(self):
        f1 = os.path.join(RDConfig.RDBaseDir,'Code','Shape','test_data','ace.sdf')
        d1 = os.path.join(RDConfig.RDBaseDir,'Code','Shape','test_data','output_ace_1.tab')

        mol1 = Chem.SDMolSupplier(f1).next()
        mol2 = Chem.SDMolSupplier(f1).next()
        score = rdkit.Chem.rdShape.Align(mol1, mol2, maxIter=10)
        self.failUnlessAlmostEqual(score,1.0,1)

        
        #df1 = file(d1,'w+')
        df1 = [x.strip().split() for x in file(d1,'r')]
        scores = [(x,float(y)) for x,y in df1]
        supp = Chem.SDMolSupplier(f1)
        tanimoto=[]
        w = Chem.SDWriter('/tmp/shapeout.sdf')
        for i,m in enumerate(supp):
            tani = rdkit.Chem.rdShape.Align(mol1, m, maxIter=0)
            tanimoto.append(tani)
            ##print >>df1,'%d %.4f'%(i,tani)
            print i,scores[i][1],'%.4f'%tani,'%.4f'%(scores[i][1]-tani)
            self.failUnlessAlmostEqual(tani,scores[i][1],1)
            m.SetProp('ShapeTani',str(tani))
            w.write(m)


if __name__ == '__main__':
  unittest.main()
