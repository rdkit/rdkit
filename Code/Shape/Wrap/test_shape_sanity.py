
from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.Chem.rdShape
import numpy.testing
import unittest, os
from rdkit import RDConfig

class TestCase(unittest.TestCase) :
    def test_shape_sanity(self):
        rod = "C#CC#CC#CC#C"
        small_rod = "C#CC#C"
        dot = "C"
        
        rod, small_rod, dot = [Chem.MolFromSmiles(smi) for smi in [rod, small_rod, dot]]
        
        [AllChem.EmbedMolecule(mol) for mol in [rod, small_rod, dot ]]
        [AllChem.UFFOptimizeMolecule(mol) for mol in [rod, small_rod, dot ]]
        
        #print >>file(os.path.join(RDConfig.RDBaseDir,'Code','Shape','test_data', 'rod.mol'),'w'),Chem.MolToMolBlock(rod)
        #print >>file(os.path.join(RDConfig.RDBaseDir,'Code','Shape','test_data', 'small_rod.mol'),'w'),Chem.MolToMolBlock(small_rod)
        #print >>file(os.path.join(RDConfig.RDBaseDir,'Code','Shape','test_data', 'dot.mol'),'w'),Chem.MolToMolBlock(dot)
        
        score = rdkit.Chem.rdShape.Align(rod, small_rod, maxIter=10)
        print "rod-small_rod", score
        self.failUnlessAlmostEqual(0.510,score,2)
        
        
        score1 = rdkit.Chem.rdShape.Align(rod, dot, maxIter=10)
        print "rod-dot", score1
        self.failUnlessAlmostEqual(0.185,score1,2)
        score2 = rdkit.Chem.rdShape.Align(small_rod, dot, maxIter=10)
        print "small_rod-dot", score2
        self.failUnlessAlmostEqual(score1*2,score2,1)
        
        score = rdkit.Chem.rdShape.Align(dot, dot, maxIter=10)
        print "dot-dot", score
        self.failUnlessAlmostEqual(1.0,score,2)

if __name__ == '__main__':
  unittest.main()
