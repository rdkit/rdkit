# $Id$
#
# Created by Greg Landrum, July 2008
#
from rdkit import RDConfig
import os
import unittest
from rdkit import DataStructs, Chem
from rdkit.Avalon import pyAvalonTools

struchk_conf_path = os.path.join(RDConfig.RDDataDir, 'struchk', '')
struchk_log_path = ''
STRUCHK_INIT = '''-ta %(struchk_conf_path)scheckfgs.trn
-or
-ca %(struchk_conf_path)scheckfgs.chk
-cc
-cl 3
-cs
-cn 999
-l %(struchk_log_path)sstruchk.log'''%locals()


def feq(v1,v2,tol=1e-4):
  return abs(v1-v2)<tol
class TestCase(unittest.TestCase):
  def setUp(self) :
    pass

  def test1(self):
    m1 = Chem.MolFromSmiles('c1cccnc1')
    smi = pyAvalonTools.GetCanonSmiles(m1)
    self.failUnless(smi=='c1ccncc1')
    smi = pyAvalonTools.GetCanonSmiles('c1cccnc1',True)
    self.failUnless(smi=='c1ccncc1')
  
  def test2(self):
    tgts=['CC1=CC(=O)C=CC1=O','c2ccc1SC(=Nc1c2)SSC4=Nc3ccccc3S4','[O-][N+](=O)c1cc(Cl)c(O)c(c1)[N+]([O-])=O',
          'N=C1NC=C(S1)[N+]([O-])=O','Nc3ccc2C(=O)c1ccccc1C(=O)c2c3',
          'OC(=O)c1ccccc1C3=C2C=CC(=O)C(Br)=C2Oc4c3ccc(O)c4Br','CN(C)C2C(=O)c1ccccc1C(=O)C=2Cl',
          'Cc3ccc2C(=O)c1ccccc1C(=O)c2c3[N+]([O-])=O',r'C/C(=N\O)/C(/C)=N/O',
          'c1ccc(cc1)P(c2ccccc2)c3ccccc3']
    d= file(os.path.join(RDConfig.RDDataDir,'NCI','first_200.props.sdf'),'r').read()
    mbs = d.split('$$$$\n')[:10]
    smis = [pyAvalonTools.GetCanonSmiles(mb,False) for mb in mbs]
    self.failUnless(smis==tgts)
    smis = [pyAvalonTools.GetCanonSmiles(smi,True) for smi in smis]
    self.failUnless(smis==tgts)

  def test3(self):
    bv = pyAvalonTools.GetAvalonFP(Chem.MolFromSmiles('c1ccccn1'))
    self.failUnlessEqual(len(bv),512)
    self.failUnlessEqual(bv.GetNumOnBits(),20)
    bv = pyAvalonTools.GetAvalonFP(Chem.MolFromSmiles('c1ccccc1'))
    self.failUnlessEqual(bv.GetNumOnBits(),8)
    bv = pyAvalonTools.GetAvalonFP(Chem.MolFromSmiles('c1nnccc1'))
    self.failUnlessEqual(bv.GetNumOnBits(),30)
    bv = pyAvalonTools.GetAvalonFP(Chem.MolFromSmiles('c1ncncc1'))
    self.failUnlessEqual(bv.GetNumOnBits(),27)

    bv = pyAvalonTools.GetAvalonFP(Chem.MolFromSmiles('c1ncncc1'),nBits=1024)
    self.failUnlessEqual(len(bv),1024)
    self.failUnless(bv.GetNumOnBits()>27)

  def test4(self):
    bv = pyAvalonTools.GetAvalonFP('c1ccccn1',True)
    self.failUnlessEqual(bv.GetNumOnBits(),20)
    bv = pyAvalonTools.GetAvalonFP('c1ccccc1',True)
    self.failUnlessEqual(bv.GetNumOnBits(),8)
    bv = pyAvalonTools.GetAvalonFP('c1nnccc1',True)
    self.failUnlessEqual(bv.GetNumOnBits(),30)
    bv = pyAvalonTools.GetAvalonFP('c1ncncc1',True)
    self.failUnlessEqual(bv.GetNumOnBits(),27)
    bv = pyAvalonTools.GetAvalonFP('c1ncncc1',True,nBits=1024)
    self.failUnlessEqual(len(bv),1024)
    self.failUnless(bv.GetNumOnBits()>27)

    bv = pyAvalonTools.GetAvalonFP(Chem.MolToMolBlock(Chem.MolFromSmiles('c1ccccn1')),False)
    self.failUnlessEqual(len(bv),512)
    self.failUnlessEqual(bv.GetNumOnBits(),20)
    bv = pyAvalonTools.GetAvalonFP(Chem.MolToMolBlock(Chem.MolFromSmiles('c1ccccc1')),False)
    self.failUnlessEqual(bv.GetNumOnBits(),8)

  def test4b(self):
    words = pyAvalonTools.GetAvalonFPAsWords(Chem.MolFromSmiles('c1ccccn1'))
    words2 = pyAvalonTools.GetAvalonFPAsWords(Chem.MolFromSmiles('Cc1ccccn1'))
    self.failUnlessEqual(len(words),len(words2))
    for i,word in enumerate(words):
      self.failUnlessEqual(word&words2[i],word)

  def test5(self):
    m = Chem.MolFromSmiles('c1ccccc1C1(CC1)N')
    pyAvalonTools.Generate2DCoords(m)
    self.failUnlessEqual(m.GetNumConformers(),1)
    self.failUnless(m.GetConformer(0).Is3D()==False)

  def test6(self):
    mb=pyAvalonTools.Generate2DCoords('c1ccccc1C1(CC1)N',True)
    m = Chem.MolFromMolBlock(mb)
    self.failUnlessEqual(m.GetNumConformers(),1)
    self.failUnless(m.GetConformer(0).Is3D()==False)

  def testRDK151(self):
    smi="C[C@H](F)Cl"
    m = Chem.MolFromSmiles(smi)
    temp = pyAvalonTools.GetCanonSmiles(smi,True)
    self.failUnlessEqual(temp,smi)
    temp = pyAvalonTools.GetCanonSmiles(m)
    self.failUnlessEqual(temp,smi)

  def testStruChk(self):
    smi_good='c1ccccc1C1(CC-C(C)C1)C'
    smi_bad='c1c(R)cccc1C1(CC-C(C)C1)C'
    r = pyAvalonTools.InitializeCheckMol(STRUCHK_INIT)
    self.failUnlessEqual(r, 0)
    (err, fixed_mol) = pyAvalonTools.CheckMolecule(smi_good, True)
    self.failUnlessEqual(err, 0)
    mol = Chem.MolFromSmiles(smi_good)   
    (err, fixed_mol)=pyAvalonTools.CheckMolecule(mol)
    self.failUnlessEqual(err, 0)

    (err, fixed_mol)=pyAvalonTools.CheckMoleculeString(smi_good,True)
    self.failUnlessEqual(err, 0)
    self.failIfEqual(fixed_mol,"")
    self.failUnless(fixed_mol.find('M  END')>0)

    (err, fixed_mol)=pyAvalonTools.CheckMolecule(smi_bad, False)
    self.failIfEqual(err, 0)
    self.failIf(fixed_mol)

    (err, fixed_mol)=pyAvalonTools.CheckMoleculeString(smi_bad, False)
    self.failIfEqual(err, 0)
    self.failIf(fixed_mol)
    pyAvalonTools.CloseCheckMolFiles()
    
if __name__ == '__main__':
    unittest.main()
