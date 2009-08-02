# $Id$
# 
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as rdMD
from rdkit import DataStructs
from rdkit import RDConfig
import unittest

def feq(v1, v2, tol=1.e-4) :
  return abs(v1-v2) < tol

class TestCase(unittest.TestCase) :
  def setUp(self):
    pass
  def testAtomPairTypes(self):
    params = rdMD.AtomPairsParameters
    mol = Chem.MolFromSmiles("C=C");
    self.failUnless(rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(0))==\
                    rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(1)))
    self.failUnless(rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(0))==\
                    1 | (1 | 1<<params.numPiBits)<<params.numBranchBits)
    
    mol = Chem.MolFromSmiles("C#CO");
    self.failUnless(rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(0))!=\
                    rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(1)))
    self.failUnless(rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(0))==\
                    1 | (2 | 1<<params.numPiBits)<<params.numBranchBits)
    self.failUnless(rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(1))==\
                    2 | (2 | 1<<params.numPiBits)<<params.numBranchBits)
    self.failUnless(rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(2))==\
                    1 | (0 | 3<<params.numPiBits)<<params.numBranchBits)
    self.failUnless(rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(1),1)==\
                    1 | (2 | 1<<params.numPiBits)<<params.numBranchBits)
    self.failUnless(rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(1),2)==\
                    0 | (2 | 1<<params.numPiBits)<<params.numBranchBits)

  def testAtomPairs(self):
    m = Chem.MolFromSmiles('CCC')
    fp1 = rdMD.GetAtomPairFingerprint(m)
    fp2 = rdMD.GetAtomPairFingerprint(m,minLength=1,maxLength=2)
    nz1 = fp1.GetNonzeroElements()
    self.failUnlessEqual(len(nz1),2)
    nz2 = fp2.GetNonzeroElements()
    self.failUnlessEqual(len(nz2),2)

    fp2 = rdMD.GetAtomPairFingerprint(m,minLength=1,maxLength=1)
    nz2 = fp2.GetNonzeroElements()
    self.failUnlessEqual(len(nz2),1)

  def testHashedAtomPairs(self):
    m = Chem.MolFromSmiles('c1ccccc1')
    fp1 = rdMD.GetHashedAtomPairFingerprint(m,2048)
    fp2 = rdMD.GetHashedAtomPairFingerprint(m,2048,1,3)
    self.failUnless(fp1==fp2)
    fp2 = rdMD.GetHashedAtomPairFingerprint(m,2048,1,2)
    sim= DataStructs.DiceSimilarity(fp1,fp2)
    self.failUnless(sim>0.0 and sim<1.0)

    m = Chem.MolFromSmiles('c1ccccn1')
    fp2 = rdMD.GetHashedAtomPairFingerprint(m,2048)
    sim= DataStructs.DiceSimilarity(fp1,fp2)
    self.failUnless(sim>0.0 and sim<1.0)
    
  def testRootedAtomPairs(self):
    m = Chem.MolFromSmiles('Oc1ccccc1')
    fp1 = rdMD.GetAtomPairFingerprint(m)
    fp2 = rdMD.GetAtomPairFingerprint(m,fromAtoms=(0,))
    nz1 = fp1.GetNonzeroElements()
    nz2 = fp2.GetNonzeroElements()
    for k,v in nz2.iteritems():
      self.failUnless(v<=nz1[k])

  def testTopologicalTorsions(self):
    mol = Chem.MolFromSmiles("CC");
    fp = rdMD.GetTopologicalTorsionFingerprint(mol)
    self.failUnless(fp.GetTotalVal()==0)
    
    mol = Chem.MolFromSmiles("CCCC");
    fp = rdMD.GetTopologicalTorsionFingerprint(mol)
    self.failUnless(fp.GetTotalVal()==1)
    fp = rdMD.GetTopologicalTorsionFingerprint(mol,3)
    self.failUnless(fp.GetTotalVal()==2)
    
    mol = Chem.MolFromSmiles("CCCO");
    fp = rdMD.GetTopologicalTorsionFingerprint(mol)
    self.failUnless(fp.GetTotalVal()==1)
    fp = rdMD.GetTopologicalTorsionFingerprint(mol,3)
    self.failUnless(fp.GetTotalVal()==2)

    mol = Chem.MolFromSmiles("CCCCCCCCCCC");
    fp = rdMD.GetTopologicalTorsionFingerprint(mol,7)
    self.failUnlessRaises(ValueError,lambda : rdMD.GetTopologicalTorsionFingerprint(mol,8))
    
  def testHashedTopologicalTorsions(self):
    mol = Chem.MolFromSmiles("c1ncccc1");
    fp1 = rdMD.GetHashedTopologicalTorsionFingerprint(mol)
    mol = Chem.MolFromSmiles("n1ccccc1");
    fp2 = rdMD.GetHashedTopologicalTorsionFingerprint(mol)
    self.failUnlessEqual(DataStructs.DiceSimilarity(fp1,fp2),1.0)

  def testRootedTorsions(self):
    m = Chem.MolFromSmiles('Oc1ccccc1')
    fp1 = rdMD.GetTopologicalTorsionFingerprint(m)
    fp2 = rdMD.GetTopologicalTorsionFingerprint(m,fromAtoms=(0,))
    nz1 = fp1.GetNonzeroElements()
    nz2 = fp2.GetNonzeroElements()
    for k,v in nz2.iteritems():
      self.failUnless(v<=nz1[k])

    m = Chem.MolFromSmiles('COCC')
    fp1 = rdMD.GetTopologicalTorsionFingerprint(m)
    self.failUnlessEqual(len(fp1.GetNonzeroElements()),1)
    fp1 = rdMD.GetTopologicalTorsionFingerprint(m,fromAtoms=(0,))
    self.failUnlessEqual(len(fp1.GetNonzeroElements()),1)
    fp1 = rdMD.GetTopologicalTorsionFingerprint(m,fromAtoms=(1,))
    self.failUnlessEqual(len(fp1.GetNonzeroElements()),0)

  def testMorganFingerprints(self):
    mol = Chem.MolFromSmiles('CC(F)(Cl)C(F)(Cl)C')
    fp = rdMD.GetMorganFingerprint(mol,0)
    self.failUnless(len(fp.GetNonzeroElements())==4)
    fp = rdMD.GetMorganFingerprint(mol,1)
    self.failUnless(len(fp.GetNonzeroElements())==8)
    fp = rdMD.GetMorganFingerprint(mol,2)
    self.failUnless(len(fp.GetNonzeroElements())==9)

    mol = Chem.MolFromSmiles('CC(F)(Cl)[C@](F)(Cl)C')
    fp = rdMD.GetMorganFingerprint(mol,0)
    self.failUnless(len(fp.GetNonzeroElements())==4)
    fp = rdMD.GetMorganFingerprint(mol,1)
    self.failUnless(len(fp.GetNonzeroElements())==9)
    fp = rdMD.GetMorganFingerprint(mol,2)
    self.failUnless(len(fp.GetNonzeroElements())==10)

    mol = Chem.MolFromSmiles('CCCCC')
    fp = rdMD.GetMorganFingerprint(mol,0,fromAtoms=(0,))
    self.failUnless(len(fp.GetNonzeroElements())==1)

    
  def testCrippen(self):
    mol = Chem.MolFromSmiles("NCO");
    contribs = rdMD._CalcCrippenContribs(mol)

  def testMolWt(self):
    mol = Chem.MolFromSmiles("C");
    amw = rdMD._CalcMolWt(mol);
    self.failUnless(feq(amw,16.043,.001));
    amw = rdMD._CalcMolWt(mol,True);
    self.failUnless(feq(amw,12.011,.001));
    mol2 = Chem.AddHs(mol);
    amw = rdMD._CalcMolWt(mol2);
    self.failUnless(feq(amw,16.043,.001));
    amw = rdMD._CalcMolWt(mol2,True);
    self.failUnless(feq(amw,12.011,.001));

  def testPairValues(self):
    import base64
    testD=(('CCCO','AQAAAAQAAAAAAIAABgAAACGECAABAAAAIoQIAAEAAABBhAgAAQAAACNEGAABAAAAQUQYAAEAAABC\nRBgAAQAAAA==\n'),
           ('CNc1ccco1','AQAAAAQAAAAAAIAAEAAAACOECgABAAAAJIQKAAIAAABBhQoAAgAAAEKFCgABAAAAIsQKAAEAAABB\nxQoAAQAAAELFCgACAAAAIYQQAAEAAABChRAAAQAAAEOFEAACAAAAYYUQAAEAAAAjhBoAAQAAAEGF\nGgABAAAAQoUaAAIAAABhhRoAAQAAAEKIGgABAAAA\n'),
           )
    for smi,txt in testD:
      pkl = base64.decodestring(txt)
      fp = rdMD.GetAtomPairFingerprint(Chem.MolFromSmiles(smi))
      fp2 = DataStructs.IntSparseIntVect(pkl)
      self.failUnlessEqual(DataStructs.DiceSimilarity(fp,fp2),1.0)
      self.failUnlessEqual(fp,fp2)
           
  def testTorsionValues(self):
    import base64
    testD=(('CCCO','AQAAAAgAAAD/////DwAAAAEAAAAAAAAAIECAAAMAAAABAAAA\n'),
           ('CNc1ccco1','AQAAAAgAAAD/////DwAAAAkAAAAAAAAAIICkSAEAAAABAAAAKVKgSQEAAAABAAAAKVCgUAEAAAAB\nAAAAKVCgUQEAAAABAAAAKVCkCAIAAAABAAAAKdCkCAIAAAABAAAAKVCgSAMAAAABAAAAKVCkSAMA\nAAABAAAAIICkSAMAAAABAAAA\n'),
           )
    for smi,txt in testD:
      pkl = base64.decodestring(txt)
      fp = rdMD.GetTopologicalTorsionFingerprint(Chem.MolFromSmiles(smi))
      fp2 = DataStructs.LongSparseIntVect(pkl)
      self.failUnlessEqual(DataStructs.DiceSimilarity(fp,fp2),1.0)
      self.failUnlessEqual(fp,fp2)
           
  def test24RDKFingerprint(self):
    from rdkit import DataStructs
    m1 = Chem.MolFromSmiles('C1=CC=CC=C1')
    fp1 = Chem.RDKFingerprint(m1)
    self.failUnless(len(fp1)==2048)
    m2 = Chem.MolFromSmiles('C1=CC=CC=C1')
    fp2 = Chem.RDKFingerprint(m2)

    tmp = DataStructs.TanimotoSimilarity(fp1,fp2)
    self.failUnless(tmp==1.0,tmp)

    m2 = Chem.MolFromSmiles('C1=CC=CC=N1')
    fp2 = Chem.RDKFingerprint(m2)
    self.failUnless(len(fp2)==2048)
    tmp = DataStructs.TanimotoSimilarity(fp1,fp2)
    self.failUnless(tmp<1.0,tmp)
    self.failUnless(tmp>0.0,tmp)

    fp3 = Chem.RDKFingerprint(m1,tgtDensity=0.3)
    self.failUnless(len(fp3)<2048)
    
  def test55LayeredFingerprint(self):
    m1 = Chem.MolFromSmiles('CC(C)C')
    fp1 = Chem.LayeredFingerprint(m1)
    self.failUnlessEqual(len(fp1),2048)
    atomCounts=[0]*m1.GetNumAtoms()
    fp2 = Chem.LayeredFingerprint(m1,atomCounts=atomCounts)
    self.failUnlessEqual(fp1,fp2)
    self.failUnlessEqual(atomCounts,[4,7,4,4])

    fp2 = Chem.LayeredFingerprint(m1,atomCounts=atomCounts)
    self.failUnlessEqual(fp1,fp2)
    self.failUnlessEqual(atomCounts,[8,14,8,8])

    pbv=DataStructs.ExplicitBitVect(2048)
    fp3 = Chem.LayeredFingerprint(m1,setOnlyBits=pbv)
    self.failUnlessEqual(fp3.GetNumOnBits(),0)

    fp3 = Chem.LayeredFingerprint(m1,setOnlyBits=fp2)
    self.failUnlessEqual(fp3,fp2)

    m2=Chem.MolFromSmiles('CC')
    fp4 = Chem.LayeredFingerprint(m2)
    atomCounts=[0]*m1.GetNumAtoms()
    fp3 = Chem.LayeredFingerprint(m1,setOnlyBits=fp4,atomCounts=atomCounts)
    self.failUnlessEqual(atomCounts,[1,3,1,1])

    m2=Chem.MolFromSmiles('CCC')
    fp4 = Chem.LayeredFingerprint(m2)
    atomCounts=[0]*m1.GetNumAtoms()
    fp3 = Chem.LayeredFingerprint(m1,setOnlyBits=fp4,atomCounts=atomCounts)
    self.failUnlessEqual(atomCounts,[3,6,3,3])

    

if __name__ == '__main__':
  unittest.main()
