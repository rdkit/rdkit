from rdkit import RDConfig
import os,sys,math
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

    def testGetSetBondLength(self):
        file = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'MolTransforms',
                             'test_data', '3-cyclohexylpyridine.mol')

        m = Chem.MolFromMolFile(file, True, False)
        conf = m.GetConformer()
        dist = rdmt.GetBondLength(conf, 0, 19)
        self.failUnlessAlmostEqual(dist, 1.36, 2)
        rdmt.SetBondLength(conf, 0, 19, 2.5)
        dist = rdmt.GetBondLength(conf, 0, 19)
        self.failUnlessAlmostEqual(dist, 2.5, 1)
        rdmt.SetBondLength(conf, 19, 0, 3.0)
        dist = rdmt.GetBondLength(conf, 0, 19)
        self.failUnlessAlmostEqual(dist, 3.0, 1)

    def testGetSetAngle(self):
        file = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'MolTransforms',
                             'test_data', '3-cyclohexylpyridine.mol')

        m = Chem.MolFromMolFile(file, True, False)
        conf = m.GetConformer()
        angle = rdmt.GetAngleDeg(conf, 0, 19, 21)
        self.failUnlessAlmostEqual(angle, 109.7, 1)
        rdmt.SetAngleDeg(conf, 0, 19, 21, 125.0)
        angle = rdmt.GetAngleDeg(conf, 0, 19, 21);
        self.failUnlessAlmostEqual(angle, 125.0, 1)
        rdmt.SetAngleRad(conf, 21, 19, 0, math.pi / 2.)
        angle = rdmt.GetAngleRad(conf, 0, 19, 21)
        self.failUnlessAlmostEqual(angle, math.pi / 2., 1)
        angle = rdmt.GetAngleDeg(conf, 0, 19, 21)
        self.failUnlessAlmostEqual(angle, 90.0, 1)

    def testGetSetDihedral(self):
        file = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'MolTransforms',
                             'test_data', '3-cyclohexylpyridine.mol')

        m = Chem.MolFromMolFile(file, True, False)
        conf = m.GetConformer()
        dihedral = rdmt.GetDihedralDeg(conf, 0, 19, 21, 24)
        self.failUnlessAlmostEqual(dihedral, 176.05, 2)
        rdmt.SetDihedralDeg(conf, 8, 0, 19, 21, 65.0)
        dihedral = rdmt.GetDihedralDeg(conf, 8, 0, 19, 21)
        self.failUnlessAlmostEqual(dihedral, 65.0, 1)
        rdmt.SetDihedralDeg(conf, 8, 0, 19, 21, -130.0)
        dihedral = rdmt.GetDihedralDeg(conf, 8, 0, 19, 21)
        self.failUnlessAlmostEqual(dihedral, -130.0, 1)
        rdmt.SetDihedralRad(conf, 21, 19, 0, 8, -2. / 3. * math.pi)
        dihedral = rdmt.GetDihedralRad(conf, 8, 0, 19, 21)
        self.failUnlessAlmostEqual(dihedral, -2. / 3. * math.pi, 1)
        dihedral = rdmt.GetDihedralDeg(conf, 8, 0, 19, 21)
        self.failUnlessAlmostEqual(dihedral, -120.0, 1)
            
if __name__=="__main__":
    unittest.main()
