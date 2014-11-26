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
STRUCHK_INIT = '''-tm
-ta %(struchk_conf_path)scheckfgs.trn
-tm
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
    self.assertTrue(smi=='c1ccncc1')
    smi = pyAvalonTools.GetCanonSmiles('c1cccnc1',True)
    self.assertTrue(smi=='c1ccncc1')
  
  def test2(self):
    tgts=['CC1=CC(=O)C=CC1=O','c2ccc1SC(=Nc1c2)SSC4=Nc3ccccc3S4','[O-][N+](=O)c1cc(Cl)c(O)c(c1)[N+]([O-])=O',
          'N=C1NC=C(S1)[N+]([O-])=O','Nc3ccc2C(=O)c1ccccc1C(=O)c2c3',
          'OC(=O)c1ccccc1C3=C2C=CC(=O)C(Br)=C2Oc4c3ccc(O)c4Br','CN(C)C2C(=O)c1ccccc1C(=O)C=2Cl',
          'Cc3ccc2C(=O)c1ccccc1C(=O)c2c3[N+]([O-])=O',r'C/C(=N\O)/C(/C)=N/O',
          'c1ccc(cc1)P(c2ccccc2)c3ccccc3']
    with open(os.path.join(RDConfig.RDDataDir,'NCI','first_200.props.sdf'),'r') as f:
      d = f.read()
    mbs = d.split('$$$$\n')[:10]
    smis = [pyAvalonTools.GetCanonSmiles(mb,False) for mb in mbs]
    self.assertTrue(smis==tgts)
    smis = [pyAvalonTools.GetCanonSmiles(smi,True) for smi in smis]
    self.assertTrue(smis==tgts)

  def test3(self):
    bv = pyAvalonTools.GetAvalonFP(Chem.MolFromSmiles('c1ccccn1'))
    self.assertEqual(len(bv),512)
    self.assertEqual(bv.GetNumOnBits(),20)
    bv = pyAvalonTools.GetAvalonFP(Chem.MolFromSmiles('c1ccccc1'))
    self.assertEqual(bv.GetNumOnBits(),8)
    bv = pyAvalonTools.GetAvalonFP(Chem.MolFromSmiles('c1nnccc1'))
    self.assertEqual(bv.GetNumOnBits(),30)
    bv = pyAvalonTools.GetAvalonFP(Chem.MolFromSmiles('c1ncncc1'))
    self.assertEqual(bv.GetNumOnBits(),27)

    bv = pyAvalonTools.GetAvalonFP(Chem.MolFromSmiles('c1ncncc1'),nBits=1024)
    self.assertEqual(len(bv),1024)
    self.assertTrue(bv.GetNumOnBits()>27)

  def test4(self):
    bv = pyAvalonTools.GetAvalonFP('c1ccccn1',True)
    self.assertEqual(bv.GetNumOnBits(),20)
    bv = pyAvalonTools.GetAvalonFP('c1ccccc1',True)
    self.assertEqual(bv.GetNumOnBits(),8)
    bv = pyAvalonTools.GetAvalonFP('c1nnccc1',True)
    self.assertEqual(bv.GetNumOnBits(),30)
    bv = pyAvalonTools.GetAvalonFP('c1ncncc1',True)
    self.assertEqual(bv.GetNumOnBits(),27)
    bv = pyAvalonTools.GetAvalonFP('c1ncncc1',True,nBits=1024)
    self.assertEqual(len(bv),1024)
    self.assertTrue(bv.GetNumOnBits()>27)

    bv = pyAvalonTools.GetAvalonFP(Chem.MolToMolBlock(Chem.MolFromSmiles('c1ccccn1')),False)
    self.assertEqual(len(bv),512)
    self.assertEqual(bv.GetNumOnBits(),20)
    bv = pyAvalonTools.GetAvalonFP(Chem.MolToMolBlock(Chem.MolFromSmiles('c1ccccc1')),False)
    self.assertEqual(bv.GetNumOnBits(),8)

  def test4b(self):
    words = pyAvalonTools.GetAvalonFPAsWords(Chem.MolFromSmiles('c1ccccn1'))
    words2 = pyAvalonTools.GetAvalonFPAsWords(Chem.MolFromSmiles('Cc1ccccn1'))
    self.assertEqual(len(words),len(words2))
    for i,word in enumerate(words):
      self.assertEqual(word&words2[i],word)

  def test5(self):
    m = Chem.MolFromSmiles('c1ccccc1C1(CC1)N')
    pyAvalonTools.Generate2DCoords(m)
    self.assertEqual(m.GetNumConformers(),1)
    self.assertTrue(m.GetConformer(0).Is3D()==False)

  def test6(self):
    mb=pyAvalonTools.Generate2DCoords('c1ccccc1C1(CC1)N',True)
    m = Chem.MolFromMolBlock(mb)
    self.assertEqual(m.GetNumConformers(),1)
    self.assertTrue(m.GetConformer(0).Is3D()==False)

  def testRDK151(self):
    smi="C[C@H](F)Cl"
    m = Chem.MolFromSmiles(smi)
    temp = pyAvalonTools.GetCanonSmiles(smi,True)
    self.assertEqual(temp,smi)
    temp = pyAvalonTools.GetCanonSmiles(m)
    self.assertEqual(temp,smi)

  def testStruChk(self):
    smi_good='c1ccccc1C1(CC-C(C)C1)C'
    smi_bad='c1c(R)cccc1C1(CC-C(C)C1)C'
    r = pyAvalonTools.InitializeCheckMol(STRUCHK_INIT)
    self.assertEqual(r, 0)
    (err, fixed_mol) = pyAvalonTools.CheckMolecule(smi_good, True)
    self.assertEqual(err, 0)
    mol = Chem.MolFromSmiles(smi_good)   
    (err, fixed_mol)=pyAvalonTools.CheckMolecule(mol)
    self.assertEqual(err, 0)

    (err, fixed_mol)=pyAvalonTools.CheckMoleculeString(smi_good,True)
    self.assertEqual(err, 0)
    self.assertNotEqual(fixed_mol,"")
    self.assertTrue(fixed_mol.find('M  END')>0)

    (err, fixed_mol)=pyAvalonTools.CheckMolecule(smi_bad, False)
    self.assertNotEqual(err, 0)
    self.assertFalse(fixed_mol)

    (err, fixed_mol)=pyAvalonTools.CheckMoleculeString(smi_bad, False)
    self.assertNotEqual(err, 0)
    self.assertFalse(fixed_mol)
    pyAvalonTools.CloseCheckMolFiles()
    
#   def testIsotopeBug(self):
#     mb="""D isotope problem.mol
#   Mrv0541 08141217122D          

#   4  3  0  0  0  0            999 V2000
#    -3.2705    0.5304    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    -2.5561    0.9429    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    -1.8416    0.5304    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
#    -2.5561    1.7679    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#   1  2  1  0  0  0  0
#   2  3  1  0  0  0  0
#   2  4  1  0  0  0  0
# M  ISO  1   3   2
# M  END
# """
#     csmi = pyAvalonTools.GetCanonSmiles(mb,False)
#     self.assertEqual(csmi,'[2H]C(C)C')
#     mb="""D isotope problem.mol
#   Mrv0541 08141217122D          

#   4  3  0  0  0  0            999 V2000
#    -3.2705    0.5304    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    -2.5561    0.9429    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    -1.8416    0.5304    0.0000 H   2  0  0  0  0  0  0  0  0  0  0  0
#    -2.5561    1.7679    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#   1  2  1  0  0  0  0
#   2  3  1  0  0  0  0
#   2  4  1  0  0  0  0
# M  ISO  1   3   2
# M  END
# """
#     csmi = pyAvalonTools.GetCanonSmiles(mb,False)
#     self.assertEqual(csmi,'[2H]C(C)C')

#     mb="""D isotope problem.mol
#   Mrv0541 08141217122D          

#   4  3  0  0  0  0            999 V2000
#    -3.2705    0.5304    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    -2.5561    0.9429    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    -1.8416    0.5304    0.0000 D   0  0  0  0  0  0  0  0  0  0  0  0
#    -2.5561    1.7679    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#   1  2  1  0  0  0  0
#   2  3  1  0  0  0  0
#   2  4  1  0  0  0  0
# M  END
# """
#     csmi = pyAvalonTools.GetCanonSmiles(mb,False)
#     self.assertEqual(csmi,'[2H]C(C)C')
    
#   def testChiralPBug(self):
#     mb="""Untitled Document-1
#   Mrv0541 08161213182D          

#   5  4  0  0  0  0            999 V2000
#    -1.1196    1.1491    0.0000 P   0  0  2  0  0  0  0  0  0  0  0  0
#    -0.4052    1.5616    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
#    -1.9446    1.1491    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
#    -1.3332    1.9460    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
#    -0.7071    0.4346    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
#   1  2  1  0  0  0  0
#   1  3  1  0  0  0  0
#   1  4  2  0  0  0  0
#   1  5  1  1  0  0  0
# M  END
# """
#     r = pyAvalonTools.InitializeCheckMol(STRUCHK_INIT)
#     self.assertEqual(r, 0)
#     (err, fixed_mol) = pyAvalonTools.CheckMolecule(mb, False)
#     self.assertEqual(err, 0)
#     self.assertTrue(fixed_mol)
#     self.assertNotEqual(fixed_mol.GetAtomWithIdx(0).GetChiralTag(),Chem.rdchem.ChiralType.CHI_UNSPECIFIED)

if __name__ == '__main__':
    unittest.main()
