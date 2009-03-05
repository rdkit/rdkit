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

  def testHashedAtomPairs(self):
    m = Chem.MolFromSmiles('c1ccccc1')
    fp1 = rdMD.GetHashedAtomPairFingerprint(m)
    m = Chem.MolFromSmiles('c1ccccn1')
    fp2 = rdMD.GetHashedAtomPairFingerprint(m)
    sim= DataStructs.TanimotoSimilarity(fp1,fp2)
    self.failUnless(sim>0.0 and sim<1.0)
    
  def testRootedAtomPairs(self):
    m = Chem.MolFromSmiles('Oc1ccccc1')
    fp1 = rdMD.GetAtomPairFingerprint(m)
    fp2 = rdMD.GetAtomPairFingerprint(m,includingAtoms=(0,))
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
    
  def testRootedTorsions(self):
    m = Chem.MolFromSmiles('Oc1ccccc1')
    fp1 = rdMD.GetTopologicalTorsionFingerprint(m)
    fp2 = rdMD.GetTopologicalTorsionFingerprint(m,includingAtoms=(0,))
    nz1 = fp1.GetNonzeroElements()
    nz2 = fp2.GetNonzeroElements()
    for k,v in nz2.iteritems():
      self.failUnless(v<=nz1[k])

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


if __name__ == '__main__':
  unittest.main()
