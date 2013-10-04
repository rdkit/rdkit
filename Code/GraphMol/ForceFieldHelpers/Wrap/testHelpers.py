from rdkit import Chem
from rdkit.Chem import ChemicalForceFields
from rdkit import RDConfig
import unittest
import os

def feq(v1,v2,tol2=1e-4):
  return abs(v1-v2)<=tol2

class TestCase(unittest.TestCase):
  def setUp(self) :
    self.dirName=os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
                              'ForceFieldHelpers','UFF','test_data')

  def test1(self) :
    fName = os.path.join(self.dirName,'benzene.mol')
    m = Chem.MolFromMolFile(fName)
    self.failIf(ChemicalForceFields.UFFOptimizeMolecule(m))
    # make sure that keyword arguments work:
    m = Chem.MolFromMolFile(fName)
    self.failUnless(ChemicalForceFields.UFFOptimizeMolecule(m,maxIters=1))

    m = Chem.MolFromMolFile(fName)
    self.failIf(ChemicalForceFields.UFFOptimizeMolecule(m,vdwThresh=2.0))

    m = Chem.MolFromMolFile(fName)
    self.failIf(ChemicalForceFields.UFFOptimizeMolecule(m,confId=-1))

    m = Chem.MolFromMolFile(fName)
    self.failUnlessRaises(ValueError,lambda :ChemicalForceFields.UFFOptimizeMolecule(m,confId=1))


  def test2(self) :
    fName = os.path.join(self.dirName,'benzene.mol')
    m = Chem.MolFromMolFile(fName)
    ff = ChemicalForceFields.UFFGetMoleculeForceField(m)
    self.failUnless(ff)
    e1 = ff.CalcEnergy()
    r = ff.Minimize()
    self.failUnless(r==0)
    e2 = ff.CalcEnergy()
    self.failUnless(e2<e1);
    
    # test keyword args:
    r = ff.Minimize(forceTol=1e-8)
    self.failUnless(r==0)

    # test keyword args:
    r = ff.Minimize(energyTol=1e-3)
    self.failUnless(r==0)

  def test3(self) :
    molB = """


  4  4  0  0  0  0  0  0  0  0999 V2000
   -0.8500    0.4512   -0.6671 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3307   -0.9436   -0.3641 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6796   -0.4074    0.5894 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5011    0.8998   -0.1231 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
  1  4  1  0
M  END"""
    m = Chem.MolFromMolBlock(molB)
    
    ff = ChemicalForceFields.UFFGetMoleculeForceField(m)
    self.failUnless(ff)
    e1 = ff.CalcEnergy()
    r = ff.Minimize()
    self.failUnless(r==0)
    e2 = ff.CalcEnergy()
    self.failUnless(e2<e1);
    
  def test4(self) :
    m = Chem.MolFromSmiles('[Cu](C)(C)(C)(C)C')
    self.failIf(ChemicalForceFields.UFFHasAllMoleculeParams(m))

    m = Chem.MolFromSmiles('C(C)(C)(C)C')
    self.failUnless(ChemicalForceFields.UFFHasAllMoleculeParams(m))


  def test5(self) :
    fName = os.path.join(self.dirName,'benzene.mol')
    m = Chem.MolFromMolFile(fName)
    self.failIf(ChemicalForceFields.MMFFOptimizeMolecule(m))
    # make sure that keyword arguments work:
    m = Chem.MolFromMolFile(fName)
    self.failUnless(ChemicalForceFields.MMFFOptimizeMolecule(m, maxIters = 1))

    m = Chem.MolFromMolFile(fName)
    self.failIf(ChemicalForceFields.MMFFOptimizeMolecule(m, nonBondedThresh = 2.0))

    m = Chem.MolFromMolFile(fName)
    self.failIf(ChemicalForceFields.MMFFOptimizeMolecule(m, confId = -1))

    m = Chem.MolFromMolFile(fName)
    self.failUnlessRaises(ValueError, lambda :ChemicalForceFields.MMFFOptimizeMolecule(m, confId = 1))


  def test6(self) :
    fName = os.path.join(self.dirName, 'benzene.mol')
    m = Chem.MolFromMolFile(fName)
    mp = ChemicalForceFields.MMFFGetMoleculeProperties(m)
    ff = ChemicalForceFields.MMFFGetMoleculeForceField(m, mp)
    self.failUnless(ff)
    e1 = ff.CalcEnergy()
    r = ff.Minimize()
    self.failUnless(r == 0)
    e2 = ff.CalcEnergy()
    self.failUnless(e2 < e1);
    
    # test keyword args:
    r = ff.Minimize(forceTol = 1.0e-8)
    self.failUnless(r == 0)

    # test keyword args:
    r = ff.Minimize(energyTol = 1.0e-3)
    self.failUnless(r == 0)

  def test7(self) :
    molB = """


  4  4  0  0  0  0  0  0  0  0999 V2000
   -0.8500    0.4512   -0.6671 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3307   -0.9436   -0.3641 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6796   -0.4074    0.5894 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.5011    0.8998   -0.1231 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
  1  4  1  0
M  END"""
    m = Chem.MolFromMolBlock(molB)
    
    mp = ChemicalForceFields.MMFFGetMoleculeProperties(m)
    ff = ChemicalForceFields.MMFFGetMoleculeForceField(m, mp)
    self.failUnless(ff)
    e1 = ff.CalcEnergy()
    r = ff.Minimize()
    self.failUnless(r == 0)
    e2 = ff.CalcEnergy()
    self.failUnless(e2 < e1);
    
  def test8(self) :
    m = Chem.MolFromSmiles('[Cu](C)(C)(C)(C)C')
    self.failIf(ChemicalForceFields.MMFFHasAllMoleculeParams(m))

    m = Chem.MolFromSmiles('C(C)(C)(C)C')
    self.failUnless(ChemicalForceFields.MMFFHasAllMoleculeParams(m))





if __name__== '__main__':
    unittest.main()

    
