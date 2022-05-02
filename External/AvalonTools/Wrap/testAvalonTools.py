# $Id$
#
# Created by Greg Landrum, July 2008
#

from rdkit import RDConfig
import os, sys
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
-l %(struchk_log_path)sstruchk.log''' % locals()

STRUCHK_INIT_IN_MEMORY_LOGGING = '''-tm
-ta %(struchk_conf_path)scheckfgs.trn
-tm
-or
-ca %(struchk_conf_path)scheckfgs.chk
-cc
-cl 3
-cs
-cn 999''' % locals()

#
# PubChem test molecule converted from the RCSB
#  mondo errors here.
atom_clash = """2FV9_002_B_2
  RCSB PDB01151502013D
Coordinates from PDB:2FV9:B:2 Model:1 without hydrogens
 32 32  0  0  0  0            999 V2000
   46.8220   28.7360   39.6060   C 0  0  0  0  0  0  0  0  0  0  0  0
   47.3620   28.0340   38.3430   C 0  0  0  0  0  0  0  0  0  0  0  0
   47.5920   29.0540   37.2270   C 0  0  0  0  0  0  0  0  0  0  0  0
   48.4130   28.4900   36.0770   C 0  0  0  0  0  0  0  0  0  0  0  0
   47.1640   31.0380   40.1700   C 0  0  0  0  0  0  0  0  0  0  0  0
   46.6160   27.6990   40.7140   C 0  0  0  0  0  0  0  0  0  0  0  0
   45.0330   26.1340   41.6460   C 0  0  0  0  0  0  0  0  0  0  0  0
   44.1100   26.4530   42.8300   C 0  0  0  0  0  0  0  0  0  0  0  0
   49.2550   34.1450   39.3140   C 0  0  0  0  0  0  0  0  0  0  0  0
   48.2670   33.0140   39.1740   C 0  0  0  0  0  0  0  0  0  0  0  0
   44.3710   25.0810   40.7680   C 0  0  0  0  0  0  0  0  0  0  0  0
   46.3620   26.9550   37.9090   C 0  0  0  0  0  0  0  0  0  0  0  0
   47.7210   29.8260   39.9930   N 0  0  0  0  0  0  0  0  0  0  0  0
   47.5720   27.2580   41.3530   O 0  0  0  0  0  0  0  0  0  0  0  0
   44.7760   23.9020   40.8430   O 0  0  0  0  0  0  0  0  0  0  0  0
   49.7640   36.2780   39.7670   O 0  0  0  0  0  0  0  0  0  0  0  0
   44.8290   27.0700   44.0370   C 0  0  0  0  0  0  0  0  0  0  0  0
   46.0700   26.2600   44.4160   C 0  0  0  0  0  0  0  0  0  0  0  0
   43.8840   27.1450   45.2400   C 0  0  0  0  0  0  0  0  0  0  0  0
   48.7580   35.3740   39.5170   N 0  0  0  0  0  0  0  0  0  0  0  0
   50.4620   33.9280   39.2410   O 0  0  0  0  0  0  0  0  0  0  0  0
   48.1640   32.1600   40.4390   C 0  0  0  0  0  0  0  0  0  0  0  0
   47.7340   32.9890   41.6620   C 0  0  0  0  0  0  0  0  0  0  0  0
   45.9630   31.2790   40.1090   O 0  0  0  0  0  0  0  0  0  0  0  0
   45.3280   27.3570   40.9080   N 0  0  0  0  0  0  0  0  0  0  0  0
   43.4500   25.4420   40.0020   O 0  0  0  0  0  0  0  0  0  0  0  0
   47.8720   32.1800   42.9370   C 0  0  0  0  0  0  0  0  0  0  0  0
   49.0780   32.2360   43.6910   C 0  0  0  0  0  0  0  0  0  0  0  0
   49.1970   31.5090   44.9110   C 0  0  0  0  0  0  0  0  0  0  0  0
   48.1080   30.7170   45.3760   C 0  0  0  0  0  0  0  0  0  0  0  0
   46.9120   30.6330   44.6100   C 0  0  0  0  0  0  0  0  0  0  0  0
   46.7920   31.3650   43.3900   C 0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  6  1  0  0  0  0
  1 13  1  0  0  0  0
  2  3  1  0  0  0  0
  2 12  1  0  0  0  0
  3  4  1  0  0  0  0
  5 13  1  0  0  0  0
  5 22  1  0  0  0  0
  5 24  2  0  0  0  0
  6 14  2  0  0  0  0
  6 25  1  0  0  0  0
  7  8  1  0  0  0  0
  7 11  1  0  0  0  0
  7 25  1  0  0  0  0
  8 17  1  0  0  0  0
  9 10  1  0  0  0  0
  9 20  1  0  0  0  0
  9 21  2  0  0  0  0
 10 22  1  0  0  0  0
 11 15  2  0  0  0  0
 11 26  1  0  0  0  0
 16 20  1  0  0  0  0
 17 18  1  0  0  0  0
 17 19  1  0  0  0  0
 22 23  1  0  0  0  0
 23 27  1  0  0  0  0
 27 28  2  0  0  0  0
 27 32  1  0  0  0  0
 28 29  1  0  0  0  0
 29 30  2  0  0  0  0
 30 31  1  0  0  0  0
 31 32  2  0  0  0  0
M  END
>  <InstanceId>
2FV9_002_B_2

>  <ChemCompId>
002

>  <PdbId>
2FV9

>  <ChainId>
B

>  <ResidueNumber>
2

>  <InsertionCode>


>  <Model>
1

>  <AltIds>


>  <MissingHeavyAtoms>
0

>  <ObservedFormula>
C23 N3 O6

>  <Name>
N-[(2R)-2-BENZYL-4-(HYDROXYAMINO)-4-OXOBUTANOYL]-L-ISOLEUCYL-L-LEUCINE

>  <SystematicName>
(2S)-2-[[(2S,3S)-2-[[(2R)-4-(hydroxyamino)-4-oxo-2-(phenylmethyl)butanoyl]amino]-3-methyl-pentanoyl]amino]-4-methyl-pentanoic acid

>  <Synonyms>


>  <Type>
NON-POLYMER

>  <Formula>
C23 H35 N3 O6

>  <MolecularWeight>
449.541

>  <ModifiedDate>
2011-06-04

>  <Parent>


>  <OneLetterCode>


>  <SubcomponentList>


>  <AmbiguousFlag>


>  <InChI>
InChI=1S/C23H35N3O6/c1-5-15(4)20(22(29)24-18(23(30)31)11-14(2)3)25-21(28)17(13-19(27)26-32)12-16-9-7-6-8-10-16/h6-10,14-15,17-18,20,32H,5,11-13H2,1-4H3,(H,24,29)(H,25,28)(H,26,27)(H,30,31)/t15-,17+,18-,20-/m0/s1

>  <InChIKey>
MWZOULASPWUGJJ-NFBUACBFSA-N

>  <SMILES>
CC[C@H](C)[C@@H](C(=O)N[C@@H](CC(C)C)C(=O)O)NC(=O)[C@H](Cc1ccccc1)CC(=O)NO

$$$$
"""

def feq(v1, v2, tol=1e-4):
  return abs(v1 - v2) < tol


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def test1(self):
    m1 = Chem.MolFromSmiles('c1cccnc1')
    smi = pyAvalonTools.GetCanonSmiles(m1)
    self.assertTrue(smi == 'c1ccncc1')
    smi = pyAvalonTools.GetCanonSmiles('c1cccnc1', True)
    self.assertTrue(smi == 'c1ccncc1')

  def test2(self):
    tgts = ['CC1=CC(=O)C=CC1=O', 'c2ccc1SC(=Nc1c2)SSC4=Nc3ccccc3S4',
            '[O-][N+](=O)c1cc(Cl)c(O)c(c1)[N+]([O-])=O', 'N=C1NC=C(S1)[N+]([O-])=O',
            'Nc3ccc2C(=O)c1ccccc1C(=O)c2c3', 'OC(=O)c1ccccc1C3=C2C=CC(=O)C(Br)=C2Oc4c3ccc(O)c4Br',
            'CN(C)C2C(=O)c1ccccc1C(=O)C=2Cl', 'Cc3ccc2C(=O)c1ccccc1C(=O)c2c3[N+]([O-])=O',
            r'C/C(=N\O)/C(/C)=N/O', 'c1ccc(cc1)P(c2ccccc2)c3ccccc3']
    with open(os.path.join(RDConfig.RDDataDir, 'NCI', 'first_200.props.sdf'), 'r') as f:
      d = f.read()
    mbs = d.split('$$$$\n')[:10]
    smis = [pyAvalonTools.GetCanonSmiles(mb, False) for mb in mbs]
    self.assertTrue(smis == tgts)
    smis = [pyAvalonTools.GetCanonSmiles(smi, True) for smi in smis]
    self.assertTrue(smis == tgts)

  def test3(self):
    bv = pyAvalonTools.GetAvalonFP(Chem.MolFromSmiles('c1ccccn1'))
    self.assertEqual(len(bv), 512)
    self.assertEqual(bv.GetNumOnBits(), 20)
    bv = pyAvalonTools.GetAvalonFP(Chem.MolFromSmiles('c1ccccc1'))
    self.assertEqual(bv.GetNumOnBits(), 8)
    bv = pyAvalonTools.GetAvalonFP(Chem.MolFromSmiles('c1nnccc1'))
    self.assertEqual(bv.GetNumOnBits(), 30)
    bv = pyAvalonTools.GetAvalonFP(Chem.MolFromSmiles('c1ncncc1'))
    self.assertEqual(bv.GetNumOnBits(), 27)

    bv = pyAvalonTools.GetAvalonFP(Chem.MolFromSmiles('c1ncncc1'), nBits=1024)
    self.assertEqual(len(bv), 1024)
    self.assertTrue(bv.GetNumOnBits() > 27)

  def test4(self):
    bv = pyAvalonTools.GetAvalonFP('c1ccccn1', True)
    self.assertEqual(bv.GetNumOnBits(), 20)
    bv = pyAvalonTools.GetAvalonFP('c1ccccc1', True)
    self.assertEqual(bv.GetNumOnBits(), 8)
    bv = pyAvalonTools.GetAvalonFP('c1nnccc1', True)
    self.assertEqual(bv.GetNumOnBits(), 30)
    bv = pyAvalonTools.GetAvalonFP('c1ncncc1', True)
    self.assertEqual(bv.GetNumOnBits(), 27)
    bv = pyAvalonTools.GetAvalonFP('c1ncncc1', True, nBits=1024)
    self.assertEqual(len(bv), 1024)
    self.assertTrue(bv.GetNumOnBits() > 27)

    bv = pyAvalonTools.GetAvalonFP(Chem.MolToMolBlock(Chem.MolFromSmiles('c1ccccn1')), False)
    self.assertEqual(len(bv), 512)
    self.assertEqual(bv.GetNumOnBits(), 20)
    bv = pyAvalonTools.GetAvalonFP(Chem.MolToMolBlock(Chem.MolFromSmiles('c1ccccc1')), False)
    self.assertEqual(bv.GetNumOnBits(), 8)

  def test4b(self):
    words = pyAvalonTools.GetAvalonFPAsWords(Chem.MolFromSmiles('c1ccccn1'))
    words2 = pyAvalonTools.GetAvalonFPAsWords(Chem.MolFromSmiles('Cc1ccccn1'))
    self.assertEqual(len(words), len(words2))
    for i, word in enumerate(words):
      self.assertEqual(word & words2[i], word)

  def test5(self):
    m = Chem.MolFromSmiles('c1ccccc1C1(CC1)N')
    pyAvalonTools.Generate2DCoords(m)
    self.assertEqual(m.GetNumConformers(), 1)
    self.assertFalse(m.GetConformer(0).Is3D())

  def test6(self):
    mb = pyAvalonTools.Generate2DCoords('c1ccccc1C1(CC1)N', True)
    m = Chem.MolFromMolBlock(mb)
    self.assertEqual(m.GetNumConformers(), 1)
    self.assertFalse(m.GetConformer(0).Is3D())

  def testGitHub1062(self):
    s0 = r'C/C=C\C'
    m1 = Chem.MolFromSmiles(s0)
    s1 = Chem.MolToSmiles(m1)
    pyAvalonTools.Generate2DCoords(m1)
    mb = Chem.MolToMolBlock(m1)
    m2 = Chem.MolFromMolBlock(mb)
    s2 = Chem.MolToSmiles(m2)
    self.assertEqual(s1, s2)

    # repeat the test with an input smiles that is not canonical
    # to verify that the implementation is not sensitive to the
    # ordering of atoms
    s0 = r'C/C=C(F)\C'
    m1 = Chem.MolFromSmiles(s0)
    s1 = Chem.MolToSmiles(m1)
    self.assertNotEqual(s1, s0)
    pyAvalonTools.Generate2DCoords(m1)
    mb = Chem.MolToMolBlock(m1)
    m2 = Chem.MolFromMolBlock(mb)
    s2 = Chem.MolToSmiles(m2)
    self.assertEqual(s1, s2)


  def testRDK151(self):
    smi = "C[C@H](F)Cl"
    m = Chem.MolFromSmiles(smi)
    temp = pyAvalonTools.GetCanonSmiles(smi, True)
    self.assertEqual(temp, smi)
    temp = pyAvalonTools.GetCanonSmiles(m)
    self.assertEqual(temp, smi)

  def testStruChk(self):
    smi_good = 'c1ccccc1C1(CC-C(C)C1)C'
    smi_bad = 'c1c(R)cccc1C1(CC-C(C)C1)C'
    r = pyAvalonTools.InitializeCheckMol(STRUCHK_INIT)
    self.assertEqual(r, 0)
    (err, fixed_mol) = pyAvalonTools.CheckMolecule(smi_good, True)
    self.assertEqual(err, 0)
    mol = Chem.MolFromSmiles(smi_good)
    (err, fixed_mol) = pyAvalonTools.CheckMolecule(mol)
    self.assertEqual(err, 0)

    (err, fixed_mol) = pyAvalonTools.CheckMoleculeString(smi_good, True)
    self.assertEqual(err, 0)
    self.assertNotEqual(fixed_mol, "")
    self.assertTrue(fixed_mol.find('M  END') > 0)

    (err, fixed_mol) = pyAvalonTools.CheckMolecule(smi_bad, False)
    self.assertNotEqual(err, 0)
    self.assertFalse(fixed_mol)

    (err, fixed_mol) = pyAvalonTools.CheckMoleculeString(smi_bad, False)
    self.assertNotEqual(err, 0)
    self.assertFalse(fixed_mol)
    pyAvalonTools.CloseCheckMolFiles()

  def testStruChkInMemoryLog(self):
    r = pyAvalonTools.InitializeCheckMol(STRUCHK_INIT_IN_MEMORY_LOGGING)
    try:
      (err, fixed_mol) = pyAvalonTools.CheckMoleculeString(atom_clash, False)
      log =  pyAvalonTools.GetCheckMolLog()
      self.assertTrue("of average bond length from bond" in log)

      # make sure that the log is cleared for the next molecule
      (err, fixed_mol) = pyAvalonTools.CheckMoleculeString("c1ccccc1", True)
      log =  pyAvalonTools.GetCheckMolLog()
      self.assertFalse(log)

    finally:
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

  def testAvalonCountFPs(self):
    # need to go to longer bit counts to avoid collisions:
    cv1 = pyAvalonTools.GetAvalonCountFP('c1ccccc1', True, nBits=6000)
    cv2 = pyAvalonTools.GetAvalonCountFP('c1ccccc1.c1ccccc1', True, nBits=6000)
    for idx, v in cv1.GetNonzeroElements().items():
      self.assertEqual(2 * v, cv2[idx])

    cv1 = pyAvalonTools.GetAvalonCountFP(Chem.MolFromSmiles('c1ccccc1'), nBits=6000)
    cv2 = pyAvalonTools.GetAvalonCountFP(Chem.MolFromSmiles('c1ccccc1.c1ccccc1'), nBits=6000)
    for idx, v in cv1.GetNonzeroElements().items():
      self.assertEqual(2 * v, cv2[idx])


if __name__ == '__main__':
  unittest.main()
