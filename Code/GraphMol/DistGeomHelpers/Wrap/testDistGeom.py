from __future__ import print_function
import unittest
import os,copy
import math
import numpy

from rdkit.six.moves import cPickle as pickle
from rdkit.six import next
from rdkit import Chem
from rdkit.Chem import rdDistGeom,ChemicalForceFields,rdMolAlign
from rdkit import RDConfig
from rdkit.Geometry import rdGeometry as geom
from rdkit.RDLogger import logger
logger=logger()

def feq(v1, v2, tol=1.e-4) :
    return abs(v1-v2) < tol

def lstEq(l1, l2, tol=1.0e-4) :
    ln = len(l1)
    if (ln != len(l2) ) :
        return 0

    for i in range(ln) :
        if abs(l1[i] - l2[i]) > tol :
            return 0
    return 1

def compareWithOld(smilesFile, sdFile) :
    smiSup = Chem.SmilesMolSupplier(smilesFile, ",", 0, -1)
    sdsup = Chem.SDMolSupplier(sdFile)
    im = 0
    for mol in smiSup :
        cid = rdDistGeom.EmbedMolecule(mol, 10,1)
        omol = sdsup[im]
        assert cid == 0
        conf = mol.GetConformer(0)
        oconf = omol.GetConformer()
        nat = mol.GetNumAtoms()
        for i in range(nat) :
            #atm = mol.GetAtomWithIdx(i)
            #oatm = omol.GetAtomWithIdx(i)
            pos = conf.GetAtomPosition(i)
            opos = oconf.GetAtomPosition(i)
            if not lstEq(pos, opos):
                return 0
        im += 1
    return 1

def compareMatrices(bm1, bm2, map, tol=1.0e-5) :
    N = numpy.shape(bm1)[0]
    for i in range(1,N):
        for j in range(i):
            l, m = map[i], map[j]
            if (l < m) :
                l, m = m,l
            if (abs(bm1[l,m] - bm2[i,j]) > tol):
                return 0
                
            if (abs(bm1[m,l] - bm2[j,i]) > tol):
                return 0
                
    return 1

def compareOrder(smi1, smi2, tol=1.0e-5) :
    m1 = Chem.MolFromSmiles(smi1)
    m2 = Chem.MolFromSmiles(smi2)
    bm1 = rdDistGeom.GetMoleculeBoundsMatrix(m1)
    bm2 = rdDistGeom.GetMoleculeBoundsMatrix(m2)
    map = m1.GetSubstructMatch(m2)
    return compareMatrices(bm1, bm2,map, tol)

def computeDist(lst1, lst2):
    res = 0.0
    for i, val in enumerate(lst1):
        res += (val - lst2[i])*(val - lst2[i])
    res = math.sqrt(res)
    return res

def computeChiralVol(pt1, pt2, pt3, pt4):
    v1 = pt1 - pt4
    v2 = pt2 - pt4
    v3 = pt3 - pt4
    cp = v2.CrossProduct(v3)
    vol = v1.DotProduct(cp)
    return vol

class TestCase(unittest.TestCase) :
    def setUp(self):
        pass

    def _test0Cdk2(self):
        fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','DistGeomHelpers',
                                            'test_data','cis_trans_cases.csv')
        
        ofile = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','DistGeomHelpers',
                                            'test_data','embedDistOpti.sdf')
        self.assertTrue(compareWithOld(fileN, ofile))

    def test1Small(self):
        #writer = Chem.SDWriter("test.sdf")
        # single double and tripple atoms cases should not fail
        mol = Chem.MolFromSmiles('O')
        rdDistGeom.EmbedMolecule(mol,10,1)
        conf = mol.GetConformer()
        self.assertTrue(lstEq(conf.GetAtomPosition(0), [0.0, 0.0, 0.0]))
        #writer.write(mol)
        
        mol = Chem.MolFromSmiles('CO')
        rdDistGeom.EmbedMolecule(mol, 10,1)
        conf = mol.GetConformer()
        self.assertTrue(lstEq(conf.GetAtomPosition(0), [0.69192, 0.0, 0.0]))
        self.assertTrue(lstEq(conf.GetAtomPosition(1), [-0.69192, 0.0, 0.0]))
        #writer.write(mol)

        mol = Chem.MolFromSmiles('CCC')
        rdDistGeom.EmbedMolecule(mol,10,1)
        conf = mol.GetConformer()
        self.assertTrue(lstEq(conf.GetAtomPosition(0), [-1.21676, -0.2989, 0.0]))
        self.assertTrue(lstEq(conf.GetAtomPosition(1), [-0.00604, 0.59337, 0.0]))
        self.assertTrue(lstEq(conf.GetAtomPosition(2), [1.22281, -0.29446, 0.0]))
        #writer.write(mol)
        
        mol = Chem.MolFromSmiles('O=C=O')
        rdDistGeom.EmbedMolecule(mol,10,1)
        conf = mol.GetConformer()
        
        #writer.write(mol)
        self.assertTrue(lstEq(conf.GetAtomPosition(0), [-1.2180, -0.06088, 0.0]))
        self.assertTrue(lstEq(conf.GetAtomPosition(1), [-0.00408, 0.12116, 0.0]))
        self.assertTrue(lstEq(conf.GetAtomPosition(2), [1.22207, -0.060276, 0.0]))
        
        mol = Chem.MolFromSmiles('C=C=C=C')
        rdDistGeom.EmbedMolecule(mol,10,1)
        conf = mol.GetConformer()
        
        #writer.write(mol)

        d1 = computeDist(conf.GetAtomPosition(0), conf.GetAtomPosition(1))
        self.assertTrue(feq(d1, 1.31, 0.01))
        d2 = computeDist(conf.GetAtomPosition(0), conf.GetAtomPosition(2))
        self.assertTrue(feq(d2, 2.59, 0.05))
        d3 = computeDist(conf.GetAtomPosition(0), conf.GetAtomPosition(3))
        self.assertTrue(feq(d3, 3.84, 0.1))
        d4 = computeDist(conf.GetAtomPosition(1), conf.GetAtomPosition(2))
        self.assertTrue(feq(d4, 1.29, 0.01))
        d5 = computeDist(conf.GetAtomPosition(1), conf.GetAtomPosition(3))
        self.assertTrue(feq(d5, 2.54, 0.1))
        d6 = computeDist(conf.GetAtomPosition(2), conf.GetAtomPosition(3))
        self.assertTrue(feq(d6, 1.31, 0.01))

    def test2Utils(self):
        mol = Chem.MolFromSmiles('CC')
        bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)
        self.assertTrue(bm[1,0]>0)
        self.assertTrue(bm[0,1]>0)
        self.assertTrue(bm[0,1]>=bm[1,0])
        self.assertTrue(bm[1,0]<1.510)
        self.assertTrue(bm[0,1]>1.510)

    def test3MultiConf(self):
        mol = Chem.MolFromSmiles("CC(C)(C)c(cc12)n[n]2C(=O)/C=C(N1)/COC")
        cids = rdDistGeom.EmbedMultipleConfs(mol,10,maxAttempts=30,randomSeed=100)
        energies = [112.98, 103.57, 110.78, 100.40, 95.37, 101.64,
                    114.72, 112.65, 124.53, 107.50]
        nenergies = []
        for cid in cids:
            ff = ChemicalForceFields.UFFGetMoleculeForceField(mol, 10.0, cid)
            ee = ff.CalcEnergy()
            nenergies.append(ee)
        #print(['%.2f'%x for x in nenergies])
        #print(nenergies)
        self.assertTrue(lstEq(energies, nenergies,tol=1e-2))
            
    def test4OrderDependence(self) :
        self.assertTrue(compareOrder("CC(C)(C)C(=O)NC(C1)CC(N2C)CCC12",
                                     "CN1C2CCC1CC(NC(=O)C(C)(C)C)C2"))
        #issue 230
        self.assertTrue(compareOrder("C#CC(C)(C)N(CN1)C\\N=C/1SC",
                                     "CSC1=NCN(C(C)(C)C#C)CN1"))
        #issue 232
        self.assertTrue(compareOrder("CC(C)(C)C(=O)NC(C1)CC(N2C)CCC12",
                                     "CN1C2CCC1CC(NC(=O)C(C)(C)C)C2"))
        
    def test5Issue285(self):
        m = Chem.MolFromSmiles('CNC=O')
        cs = rdDistGeom.EmbedMultipleConfs(m,10)
        for i,ci in enumerate(cs):
            for j in range(i+1,len(cs)):
                cj = cs[j]
                self.assertTrue(Chem.MolToMolBlock(m,confId=ci)!=Chem.MolToMolBlock(m,confId=cj))

    def test6RmsPruning(self):
        smiles = ['CC(C)CC(NC(C1[N+]CCC1)=O)C([O-])=O',
                  'CC(NC(CO)C(O)c1ccc([N+]([O-])=O)cc1)=O',
                  'CC([N+])C(NC(C)C(N1C(C=O)CCC1)=O)=O',
                  'CC(NC1C(O)C=C(C([O-])=O)OC1C(O)C(O)CO)=O',
                  'CCCC=C(NC(C1CC1(C)C)=O)C([O-])=O',
                  'OCC(O)C(O)C(Cn1c2c(cc(C)c(C)c2)nc-2c(=O)[nH]c(=O)nc12)O']


        nconfs = []
        expected = [5, 6, 6, 6, 6, 3]
        for smi in smiles:
            mol = Chem.MolFromSmiles(smi)
            cids = rdDistGeom.EmbedMultipleConfs(mol, 50, maxAttempts=30,
                                                 randomSeed=100, pruneRmsThresh=1.5)
            nconfs.append(len(cids))
            
        d = [abs(x-y) for x,y in zip(expected,nconfs)]

        self.assertTrue(max(d)<=1)

    def test6Chirality(self):
        # turn on chirality and we should get chiral volume that is pretty consistent and
        # positive
        tgtVol=13.0
        smiles = "Cl[C@](C)(F)Br"
        mol = Chem.MolFromSmiles(smiles)
        cids = rdDistGeom.EmbedMultipleConfs(mol, 30, maxAttempts=30,
                                             randomSeed=100)
        self.assertTrue(len(cids)==30)
        for cid in cids:
            conf = mol.GetConformer(cid)
            vol = computeChiralVol(conf.GetAtomPosition(0),
                                   conf.GetAtomPosition(2),
                                   conf.GetAtomPosition(3),
                                   conf.GetAtomPosition(4))
            self.assertTrue(abs(vol-tgtVol)<1)

        # turn of chirality and now we should see both chiral forms
        smiles = "ClC(C)(F)Br"
        mol = Chem.MolFromSmiles(smiles)
        cids = rdDistGeom.EmbedMultipleConfs(mol, 30, maxAttempts=30,
                                             randomSeed=120)
        self.assertTrue(len(cids)==30)
        nPos=0
        nNeg=0
        for cid in cids:
            conf = mol.GetConformer(cid)
            vol = computeChiralVol(conf.GetAtomPosition(0),
                                   conf.GetAtomPosition(2),
                                   conf.GetAtomPosition(3),
                                   conf.GetAtomPosition(4))
            self.assertTrue(abs(vol-tgtVol)<1 or abs(vol+tgtVol)<1)
            if vol<0: nNeg+=1
            else: nPos+=1
        self.assertTrue(nPos>0)
        self.assertTrue(nNeg>0)

        tgtVol=5.0
        for i in range(10):
            smiles = "Cl[C@H](F)Br"
            mol = Chem.MolFromSmiles(smiles)
            ci = rdDistGeom.EmbedMolecule(mol, 30, (i+1)*10)
            conf = mol.GetConformer(ci)
            vol = computeChiralVol(conf.GetAtomPosition(0),
                                   conf.GetAtomPosition(1),
                                   conf.GetAtomPosition(2),
                                   conf.GetAtomPosition(3))
            self.assertTrue(abs(vol-tgtVol)<1,"%s %s"%(vol,tgtVol))

        tgtVol=3.5
        expected = [-3.62, -3.67, -3.72,  3.91,  3.95,  3.98,  3.90,  3.94,  3.98,  3.91]
        nPos=0
        nNeg=0
        for i in range(30):
            smiles = "ClC(F)Br"
            mol = Chem.MolFromSmiles(smiles)
            ci = rdDistGeom.EmbedMolecule(mol, 30, (i+1)*10)
            conf = mol.GetConformer(ci)
            vol = computeChiralVol(conf.GetAtomPosition(0),
                                   conf.GetAtomPosition(1),
                                   conf.GetAtomPosition(2),
                                   conf.GetAtomPosition(3))
            self.assertTrue(abs(vol-tgtVol)<1 or abs(vol+tgtVol)<1)
            if vol<0: nNeg+=1
            else: nPos+=1

        self.assertTrue(nPos>0)
        self.assertTrue(nNeg>0)

        smiles = "Cl[C@H](F)Br"
        m = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(m)
        cids = rdDistGeom.EmbedMultipleConfs(mol, 10, maxAttempts=30,
                                             randomSeed=100)
        self.assertTrue(len(cids)==10)
        tgtVol=10.5
        for cid in cids:
            conf = mol.GetConformer(cid)
            vol = computeChiralVol(conf.GetAtomPosition(0),
                                   conf.GetAtomPosition(2),
                                   conf.GetAtomPosition(3),
                                   conf.GetAtomPosition(4))
            self.assertTrue(abs(vol-tgtVol)<2.)
        
        # let's try a little more complicated system
        expectedV1 = -2.0
        expectedV2 = -2.9
        
        for i in range(5):
            smi = "C1=CC=C(C=C1)[C@H](OC1=C[NH]N=C1)C(=O)[NH]C[C@H](Cl)C1=CC=NC=C1"
            mol = Chem.MolFromSmiles(smi)
            ci = rdDistGeom.EmbedMolecule(mol, 30, randomSeed=(i+1)*15)
            ff = ChemicalForceFields.UFFGetMoleculeForceField(mol, 10.0, ci)
            ff.Minimize()
            
            conf = mol.GetConformer(ci)
            vol1 = computeChiralVol(conf.GetAtomPosition(6),
                                   conf.GetAtomPosition(3),
                                   conf.GetAtomPosition(7),
                                   conf.GetAtomPosition(13))
            self.assertTrue(abs(vol1-expectedV1)<1 or abs(vol1+expectedV1)<1)
            if vol1<0: nNeg+=1
            else: nPos+=1


            vol2 = computeChiralVol(conf.GetAtomPosition(17),
                                    conf.GetAtomPosition(16),
                                    conf.GetAtomPosition(18),
                                    conf.GetAtomPosition(19))
            self.assertTrue(abs(vol2-expectedV2)<1 or abs(vol2+expectedV2)<1)

        # remove the chiral specification and we should see other chiral
        # forms of the compound
        expectedV1 = 2.0 #[-2.30, -2.31, -2.30,  2.30, -1.77]
        expectedV2 = 2.8 #[2.90,  2.89,  2.69, -2.90, -2.93]
        
        self.assertTrue(nPos>0)
        self.assertTrue(nNeg>0)
        for i in range(5):
            smi = "C1=CC=C(C=C1)C(OC1=C[NH]N=C1)C(=O)[NH]CC(Cl)C1=CC=NC=C1"
            mol = Chem.MolFromSmiles(smi)
            ci = rdDistGeom.EmbedMolecule(mol, 30, (i+1)*10)
            ff = ChemicalForceFields.UFFGetMoleculeForceField(mol, 10.0, ci)
            ff.Minimize()
            
            conf = mol.GetConformer(ci)
            vol1 = computeChiralVol(conf.GetAtomPosition(6),
                                   conf.GetAtomPosition(3),
                                   conf.GetAtomPosition(7),
                                   conf.GetAtomPosition(13))
            vol2 = computeChiralVol(conf.GetAtomPosition(17),
                                    conf.GetAtomPosition(16),
                                    conf.GetAtomPosition(18),
                                    conf.GetAtomPosition(19))
            self.assertTrue(abs(abs(vol1)-expectedV1)<1.0)
            self.assertTrue(abs(abs(vol2)-expectedV2)<1.0)


    def test7ConstrainedEmbedding(self):
        ofile = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','DistGeomHelpers',
                             'test_data','constrain1.sdf')
        suppl = Chem.SDMolSupplier(ofile);
        ref = next(suppl)
        probe = copy.deepcopy(ref)

        cMap={}
        for i in range(5):
            cMap[i]=ref.GetConformer().GetAtomPosition(i)
        ci = rdDistGeom.EmbedMolecule(probe,coordMap=cMap,randomSeed=23)
        self.assertTrue(ci>-1);
        algMap = list(zip(range(5),range(5)))
        ssd = rdMolAlign.AlignMol(probe,ref,atomMap=algMap)
        self.assertTrue(ssd<0.1)
            
if __name__ == '__main__':
  unittest.main()
