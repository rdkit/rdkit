# $Id$
# 
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as rdMD
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit import RDConfig
from rdkit.Geometry import rdGeometry as rdG
import unittest

def feq(v1, v2, tol=1.e-4) :
  return abs(v1-v2) < tol

class TestCase(unittest.TestCase) :
  def setUp(self):
    pass
  def testAtomPairTypes(self):
    params = rdMD.AtomPairsParameters
    mol = Chem.MolFromSmiles("C=C");
    self.assertTrue(rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(0))==\
                    rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(1)))
    self.assertTrue(rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(0))==\
                    1 | (1 | 1<<params.numPiBits)<<params.numBranchBits)
    
    mol = Chem.MolFromSmiles("C#CO");
    self.assertTrue(rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(0))!=\
                    rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(1)))
    self.assertTrue(rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(0))==\
                    1 | (2 | 1<<params.numPiBits)<<params.numBranchBits)
    self.assertTrue(rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(1))==\
                    2 | (2 | 1<<params.numPiBits)<<params.numBranchBits)
    self.assertTrue(rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(2))==\
                    1 | (0 | 3<<params.numPiBits)<<params.numBranchBits)
    self.assertTrue(rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(1),1)==\
                    1 | (2 | 1<<params.numPiBits)<<params.numBranchBits)
    self.assertTrue(rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(1),2)==\
                    0 | (2 | 1<<params.numPiBits)<<params.numBranchBits)

  def testAtomPairs(self):
    m = Chem.MolFromSmiles('CCC')
    fp1 = rdMD.GetAtomPairFingerprint(m)
    fp2 = rdMD.GetAtomPairFingerprint(m,minLength=1,maxLength=2)
    nz1 = fp1.GetNonzeroElements()
    self.assertEqual(len(nz1),2)
    nz2 = fp2.GetNonzeroElements()
    self.assertEqual(len(nz2),2)

    fp2 = rdMD.GetAtomPairFingerprint(m,minLength=1,maxLength=1)
    nz2 = fp2.GetNonzeroElements()
    self.assertEqual(len(nz2),1)

  def testHashedAtomPairs(self):
    m = Chem.MolFromSmiles('c1ccccc1')
    fp1 = rdMD.GetHashedAtomPairFingerprint(m,2048)
    fp2 = rdMD.GetHashedAtomPairFingerprint(m,2048,1,3)
    self.assertTrue(fp1==fp2)
    fp2 = rdMD.GetHashedAtomPairFingerprint(m,2048,1,2)
    sim= DataStructs.DiceSimilarity(fp1,fp2)
    self.assertTrue(sim>0.0 and sim<1.0)

    m = Chem.MolFromSmiles('c1ccccn1')
    fp2 = rdMD.GetHashedAtomPairFingerprint(m,2048)
    sim= DataStructs.DiceSimilarity(fp1,fp2)
    self.assertTrue(sim>0.0 and sim<1.0)


    m = Chem.MolFromSmiles('c1ccccc1')
    fp1 = rdMD.GetHashedAtomPairFingerprintAsBitVect(m,2048)
    m = Chem.MolFromSmiles('c1ccccn1')
    fp2 = rdMD.GetHashedAtomPairFingerprintAsBitVect(m,2048)
    sim= DataStructs.DiceSimilarity(fp1,fp2)
    self.assertTrue(sim>0.0 and sim<1.0)

    
  def testRootedAtomPairs(self):
    m = Chem.MolFromSmiles('Oc1ccccc1')
    fp1 = rdMD.GetAtomPairFingerprint(m)
    fp2 = rdMD.GetAtomPairFingerprint(m,fromAtoms=(0,))
    nz1 = fp1.GetNonzeroElements()
    nz2 = fp2.GetNonzeroElements()
    for k,v in nz2.items():
      self.assertTrue(v<=nz1[k])

  def testTopologicalTorsions(self):
    mol = Chem.MolFromSmiles("CC");
    fp = rdMD.GetTopologicalTorsionFingerprint(mol)
    self.assertTrue(fp.GetTotalVal()==0)
    
    mol = Chem.MolFromSmiles("CCCC");
    fp = rdMD.GetTopologicalTorsionFingerprint(mol)
    self.assertTrue(fp.GetTotalVal()==1)
    fp = rdMD.GetTopologicalTorsionFingerprint(mol,3)
    self.assertTrue(fp.GetTotalVal()==2)
    
    mol = Chem.MolFromSmiles("CCCO");
    fp = rdMD.GetTopologicalTorsionFingerprint(mol)
    self.assertTrue(fp.GetTotalVal()==1)
    fp = rdMD.GetTopologicalTorsionFingerprint(mol,3)
    self.assertTrue(fp.GetTotalVal()==2)

    mol = Chem.MolFromSmiles("CCCCCCCCCCC");
    fp = rdMD.GetTopologicalTorsionFingerprint(mol,7)
    self.assertRaises(ValueError,lambda : rdMD.GetTopologicalTorsionFingerprint(mol,8))
    
  def testHashedTopologicalTorsions(self):
    mol = Chem.MolFromSmiles("c1ncccc1");
    fp1 = rdMD.GetHashedTopologicalTorsionFingerprint(mol)
    mol = Chem.MolFromSmiles("n1ccccc1");
    fp2 = rdMD.GetHashedTopologicalTorsionFingerprint(mol)
    self.assertEqual(DataStructs.DiceSimilarity(fp1,fp2),1.0)

  def testRootedTorsions(self):
    m = Chem.MolFromSmiles('Oc1ccccc1')
    fp1 = rdMD.GetTopologicalTorsionFingerprint(m)
    fp2 = rdMD.GetTopologicalTorsionFingerprint(m,fromAtoms=(0,))
    nz1 = fp1.GetNonzeroElements()
    nz2 = fp2.GetNonzeroElements()
    for k,v in nz2.items():
      self.assertTrue(v<=nz1[k])

    m = Chem.MolFromSmiles('COCC')
    fp1 = rdMD.GetTopologicalTorsionFingerprint(m)
    self.assertEqual(len(fp1.GetNonzeroElements()),1)
    fp1 = rdMD.GetTopologicalTorsionFingerprint(m,fromAtoms=(0,))
    self.assertEqual(len(fp1.GetNonzeroElements()),1)
    fp1 = rdMD.GetTopologicalTorsionFingerprint(m,fromAtoms=(1,))
    self.assertEqual(len(fp1.GetNonzeroElements()),0)

  def testMorganFingerprints(self):
    mol = Chem.MolFromSmiles('CC(F)(Cl)C(F)(Cl)C')
    fp = rdMD.GetMorganFingerprint(mol,0)
    self.assertTrue(len(fp.GetNonzeroElements())==4)
    
    mol = Chem.MolFromSmiles('CC')
    fp = rdMD.GetMorganFingerprint(mol,0)
    self.assertTrue(len(fp.GetNonzeroElements())==1)
    self.assertTrue(list(fp.GetNonzeroElements().values())[0]==2)
    fp = rdMD.GetMorganFingerprint(mol,0,useCounts=False)
    self.assertTrue(len(fp.GetNonzeroElements())==1)
    self.assertTrue(list(fp.GetNonzeroElements().values())[0]==1)

    mol = Chem.MolFromSmiles('CC(F)(Cl)C(F)(Cl)C')
    fp = rdMD.GetHashedMorganFingerprint(mol,0)
    self.assertTrue(len(fp.GetNonzeroElements())==4)
    fp = rdMD.GetMorganFingerprint(mol,1)
    self.assertTrue(len(fp.GetNonzeroElements())==8)
    fp = rdMD.GetHashedMorganFingerprint(mol,1)
    self.assertTrue(len(fp.GetNonzeroElements())==8)
    fp = rdMD.GetMorganFingerprint(mol,2)
    self.assertTrue(len(fp.GetNonzeroElements())==9)

    mol = Chem.MolFromSmiles('CC(F)(Cl)[C@](F)(Cl)C')
    fp = rdMD.GetMorganFingerprint(mol,0)
    self.assertTrue(len(fp.GetNonzeroElements())==4)
    fp = rdMD.GetMorganFingerprint(mol,1)
    self.assertTrue(len(fp.GetNonzeroElements())==8)
    fp = rdMD.GetMorganFingerprint(mol,2)
    self.assertTrue(len(fp.GetNonzeroElements())==9)
    fp = rdMD.GetMorganFingerprint(mol,0,useChirality=True)
    self.assertTrue(len(fp.GetNonzeroElements())==4)
    fp = rdMD.GetMorganFingerprint(mol,1,useChirality=True)
    self.assertTrue(len(fp.GetNonzeroElements())==9)
    fp = rdMD.GetMorganFingerprint(mol,2,useChirality=True)
    self.assertTrue(len(fp.GetNonzeroElements())==10)

    mol = Chem.MolFromSmiles('CCCCC')
    fp = rdMD.GetMorganFingerprint(mol,0,fromAtoms=(0,))
    self.assertTrue(len(fp.GetNonzeroElements())==1)

    mol = Chem.MolFromSmiles('CC1CC1')
    vs1 = rdMD.GetConnectivityInvariants(mol)
    self.assertEqual(len(vs1),mol.GetNumAtoms())
    fp1 = rdMD.GetMorganFingerprint(mol,2,invariants=vs1)
    fp2 = rdMD.GetMorganFingerprint(mol,2)
    self.assertEqual(fp1,fp2)
    
    vs2 = rdMD.GetConnectivityInvariants(mol,False)
    self.assertEqual(len(vs2),mol.GetNumAtoms())
    self.assertNotEqual(vs1,vs2)
    fp1 = rdMD.GetMorganFingerprint(mol,2,invariants=vs2)
    self.assertNotEqual(fp1,fp2)

    mol = Chem.MolFromSmiles('Cc1ccccc1')
    vs1 = rdMD.GetFeatureInvariants(mol)
    self.assertEqual(len(vs1),mol.GetNumAtoms())
    self.assertEqual(vs1[0],0)
    self.assertNotEqual(vs1[1],0)
    self.assertEqual(vs1[1],vs1[2])
    self.assertEqual(vs1[1],vs1[3])
    self.assertEqual(vs1[1],vs1[4])

    mol = Chem.MolFromSmiles('FCCCl')
    vs1 = rdMD.GetFeatureInvariants(mol)
    self.assertEqual(len(vs1),mol.GetNumAtoms())
    self.assertEqual(vs1[1],0)
    self.assertEqual(vs1[2],0)
    self.assertNotEqual(vs1[0],0)
    self.assertEqual(vs1[0],vs1[3])
    
    fp1 = rdMD.GetMorganFingerprint(mol,0,invariants=vs1)
    fp2 = rdMD.GetMorganFingerprint(mol,0,useFeatures=True)
    self.assertEqual(fp1,fp2)

    
  def testCrippen(self):
    mol = Chem.MolFromSmiles("n1ccccc1CO");
    contribs = rdMD._CalcCrippenContribs(mol)
    self.assertEqual(len(contribs),mol.GetNumAtoms());

    ts = [0]*mol.GetNumAtoms()
    contribs = rdMD._CalcCrippenContribs(mol,force=True,atomTypes=ts)
    self.assertEqual(ts,[59, 25, 25, 25, 25, 28, 17, 69])

    ls = ['']*mol.GetNumAtoms()
    contribs = rdMD._CalcCrippenContribs(mol,force=True,atomTypeLabels=ls)
    self.assertEqual(ls,['N11', 'C18', 'C18', 'C18', 'C18', 'C21', 'C10', 'O2'])

  def testMolWt(self):
    mol = Chem.MolFromSmiles("C");
    amw = rdMD._CalcMolWt(mol);
    self.assertTrue(feq(amw,16.043,.001));
    amw = rdMD._CalcMolWt(mol,True);
    self.assertTrue(feq(amw,12.011,.001));
    mol2 = Chem.AddHs(mol);
    amw = rdMD._CalcMolWt(mol2);
    self.assertTrue(feq(amw,16.043,.001));
    amw = rdMD._CalcMolWt(mol2,True);
    self.assertTrue(feq(amw,12.011,.001));

    mol = Chem.MolFromSmiles("C");
    amw = rdMD.CalcExactMolWt(mol);
    self.assertTrue(feq(amw,16.031,.001));

    
  def testPairValues(self):
    import base64
    testD=(('CCCO',b'AQAAAAQAAAAAAIAABgAAACGECAABAAAAIoQIAAEAAABBhAgAAQAAACNEGAABAAAAQUQYAAEAAABC\nRBgAAQAAAA==\n'),
           ('CNc1ccco1',b'AQAAAAQAAAAAAIAAEAAAACOECgABAAAAJIQKAAIAAABBhQoAAgAAAEKFCgABAAAAIsQKAAEAAABB\nxQoAAQAAAELFCgACAAAAIYQQAAEAAABChRAAAQAAAEOFEAACAAAAYYUQAAEAAAAjhBoAAQAAAEGF\nGgABAAAAQoUaAAIAAABhhRoAAQAAAEKIGgABAAAA\n'),
           )
    for smi,txt in testD:
      pkl = base64.decodestring(txt)
      fp = rdMD.GetAtomPairFingerprint(Chem.MolFromSmiles(smi))
      fp2 = DataStructs.IntSparseIntVect(pkl)
      self.assertEqual(DataStructs.DiceSimilarity(fp,fp2),1.0)
      self.assertEqual(fp,fp2)
           
  def testTorsionValues(self):
    import base64
    testD=(('CCCO',b'AQAAAAgAAAD/////DwAAAAEAAAAAAAAAIECAAAMAAAABAAAA\n'),
           ('CNc1ccco1',b'AQAAAAgAAAD/////DwAAAAkAAAAAAAAAIICkSAEAAAABAAAAKVKgSQEAAAABAAAAKVCgUAEAAAAB\nAAAAKVCgUQEAAAABAAAAKVCkCAIAAAABAAAAKdCkCAIAAAABAAAAKVCgSAMAAAABAAAAKVCkSAMA\nAAABAAAAIICkSAMAAAABAAAA\n'),
           )
    for smi,txt in testD:
      pkl = base64.decodestring(txt)
      fp = rdMD.GetTopologicalTorsionFingerprint(Chem.MolFromSmiles(smi))
      fp2 = DataStructs.LongSparseIntVect(pkl)
      self.assertEqual(DataStructs.DiceSimilarity(fp,fp2),1.0)
      self.assertEqual(fp,fp2)

  def testAtomPairOptions(self):
    m1 = Chem.MolFromSmiles('c1ccccc1')
    m2 = Chem.MolFromSmiles('c1ccccn1')

    fp1 = rdMD.GetAtomPairFingerprint(m1)
    fp2 = rdMD.GetAtomPairFingerprint(m2)
    self.assertNotEqual(fp1,fp2)
    
    fp1 = rdMD.GetAtomPairFingerprint(m1,atomInvariants=[1]*6)
    fp2 = rdMD.GetAtomPairFingerprint(m2,atomInvariants=[1]*6)
    self.assertEqual(fp1,fp2)

    fp1 = rdMD.GetAtomPairFingerprint(m1,atomInvariants=[1]*6)
    fp2 = rdMD.GetAtomPairFingerprint(m2,atomInvariants=[2]*6)
    self.assertNotEqual(fp1,fp2)

    fp1 = rdMD.GetHashedAtomPairFingerprintAsBitVect(m1)
    fp2 = rdMD.GetHashedAtomPairFingerprintAsBitVect(m2)
    self.assertNotEqual(fp1,fp2)
    
    fp1 = rdMD.GetHashedAtomPairFingerprintAsBitVect(m1,atomInvariants=[1]*6)
    fp2 = rdMD.GetHashedAtomPairFingerprintAsBitVect(m2,atomInvariants=[1]*6)
    self.assertEqual(fp1,fp2)

    fp1 = rdMD.GetHashedAtomPairFingerprintAsBitVect(m1,atomInvariants=[1]*6)
    fp2 = rdMD.GetHashedAtomPairFingerprintAsBitVect(m2,atomInvariants=[2]*6)
    self.assertNotEqual(fp1,fp2)

    fp1 = rdMD.GetTopologicalTorsionFingerprint(m1)
    fp2 = rdMD.GetTopologicalTorsionFingerprint(m2)
    self.assertNotEqual(fp1,fp2)
    
    fp1 = rdMD.GetTopologicalTorsionFingerprint(m1,atomInvariants=[1]*6)
    fp2 = rdMD.GetTopologicalTorsionFingerprint(m2,atomInvariants=[1]*6)
    self.assertEqual(fp1,fp2)

    fp1 = rdMD.GetTopologicalTorsionFingerprint(m1,atomInvariants=[1]*6)
    fp2 = rdMD.GetTopologicalTorsionFingerprint(m2,atomInvariants=[2]*6)
    self.assertNotEqual(fp1,fp2)

    fp1 = rdMD.GetHashedTopologicalTorsionFingerprintAsBitVect(m1)
    fp2 = rdMD.GetHashedTopologicalTorsionFingerprintAsBitVect(m2)
    self.assertNotEqual(fp1,fp2)
    
    fp1 = rdMD.GetHashedTopologicalTorsionFingerprintAsBitVect(m1,atomInvariants=[1]*6)
    fp2 = rdMD.GetHashedTopologicalTorsionFingerprintAsBitVect(m2,atomInvariants=[1]*6)
    self.assertEqual(fp1,fp2)

    fp1 = rdMD.GetHashedTopologicalTorsionFingerprintAsBitVect(m1,atomInvariants=[1]*6)
    fp2 = rdMD.GetHashedTopologicalTorsionFingerprintAsBitVect(m2,atomInvariants=[2]*6)
    self.assertNotEqual(fp1,fp2)

  def testMolFormula(self):
    m = Chem.MolFromSmiles("[2H]C([3H])O")
    formula = rdMD.CalcMolFormula(m)
    self.assertEqual(formula,'CH4O')
    formula = rdMD.CalcMolFormula(m,separateIsotopes=True)
    self.assertEqual(formula,'CH2DTO')
    formula = rdMD.CalcMolFormula(m,separateIsotopes=True,abbreviateHIsotopes=False)
    self.assertEqual(formula,'CH2[2H][3H]O')

    m = Chem.MolFromSmiles("[2H][13CH2]CO")
    formula = rdMD.CalcMolFormula(m)
    self.assertEqual(formula,'C2H6O')
    formula = rdMD.CalcMolFormula(m,separateIsotopes=True)
    self.assertEqual(formula,'C[13C]H5DO')




if __name__ == '__main__':
  unittest.main()
