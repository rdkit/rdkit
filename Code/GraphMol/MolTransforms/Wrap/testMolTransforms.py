from rdkit import RDConfig
import os,sys
import unittest
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Geometry import rdGeometry as geom
from rdkit.Chem import rdMolTransforms as rdmt

def feq(v1, v2, tol=1.0e-4):
    return abs(v1-v2) < tol

def ptEq(pt, tpt, tol=0.0001):
    pt -= tpt
    return feq(pt.Length(), 0.0)

class TestCase(unittest.TestCase):
    def setUp(self) :
        pass

    def test1Canonicalization(self) :
        mol = Chem.MolFromSmiles("C")
        conf = Chem.Conformer(1)
        conf.SetAtomPosition(0, (4.0, 5.0, 6.0))
        mol.AddConformer(conf,1)

        conf = mol.GetConformer()
        pt  = rdmt.ComputeCentroid(conf)
        self.failUnless(ptEq(pt, geom.Point3D(4.0, 5.0, 6.0)))

        fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','MolTransforms',
                                            'test_data','1oir.mol')
        m = Chem.MolFromMolFile(fileN)
        cpt = rdmt.ComputeCentroid(m.GetConformer())
        trans = rdmt.ComputeCanonicalTransform(m.GetConformer(), cpt)
        trans2 = rdmt.ComputeCanonicalTransform(m.GetConformer())
        for i in range(4):
            for j in range(4) :
                self.failUnless(feq(trans[i,j], trans2[i, j]))
        rdmt.TransformConformer(m.GetConformer(), trans2)
        m2 = Chem.MolFromMolFile(fileN)
        rdmt.CanonicalizeConformer(m2.GetConformer())
        nats = m.GetNumAtoms()
        cnf1 = m.GetConformer()
        cnf2 = m2.GetConformer()
        for i in range(nats):
            p1 = list(cnf1.GetAtomPosition(i))
            p2 = list(cnf2.GetAtomPosition(i))
            self.failUnless(feq(p1[0],p2[0]))
            self.failUnless(feq(p1[1],p2[1]))
            self.failUnless(feq(p1[2],p2[2]))

        m3 = Chem.MolFromMolFile(fileN)
        rdmt.CanonicalizeMol(m3)
        cnf1 = m.GetConformer()
        cnf2 = m3.GetConformer()
        for i in range(nats):
            p1 = list(cnf1.GetAtomPosition(i))
            p2 = list(cnf2.GetAtomPosition(i))
            self.failUnless(feq(p1[0],p2[0]))
            self.failUnless(feq(p1[1],p2[1]))
            self.failUnless(feq(p1[2],p2[2]))
            
if __name__=="__main__":
    unittest.main()
