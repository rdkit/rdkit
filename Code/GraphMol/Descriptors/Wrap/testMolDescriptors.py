import Chem
from Chem import rdMolDescriptors as rdMD
import RDConfig
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
    

    
if __name__ == '__main__':
  unittest.main()
