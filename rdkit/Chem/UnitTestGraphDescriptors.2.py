# $Id$
#
#  Copyright (C) 2003-2006  Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""unit testing code for graph-theoretical descriptors

"""
from __future__ import print_function
from rdkit import RDConfig
import unittest,os.path
from rdkit import Chem
from rdkit.Chem import GraphDescriptors,MolSurf,Lipinski,Crippen

def feq(n1,n2,tol=1e-4):
  return abs(n1-n2)<=tol

class TestCase(unittest.TestCase):
  def setUp(self):
    if doLong:
      print('\n%s: '%self.shortDescription(),end='')

    
  def testBertzCTShort(self):
   """ test calculation of Bertz 'C(T)' index

   """
   data = [('C=CC=C',21.01955),
           ('O=CC=O',25.01955),
           ('FCC(=O)CF',46.7548875),
           ('O=C1C=CC(=O)C=C1',148.705216),
           ('C12C(F)=C(O)C(F)C1C(F)=C(O)C(F)2',315.250442),
           ('C12CC=CCC1C(=O)C3CC=CCC3C(=O)2',321.539522)]

   for smi,CT in data:
     m = Chem.MolFromSmiles(smi)
     newCT = GraphDescriptors.BertzCT(m, forceDMat = 1)
     assert feq(newCT,CT,1e-3),'mol %s (CT calc = %f) should have CT = %f'%(smi,newCT,CT)

  def _testBertzCTLong(self):
   """ test calculation of Bertz 'C(T)' index

    NOTE: this is a backwards compatibility test, because of the changes
      w.r.t. the treatment of aromatic atoms in the new version, we need
      to ignore molecules with aromatic rings...

     
   """
   col = 1
   with open(os.path.join(RDConfig.RDCodeDir,'Chem','test_data','PP_descrs_regress.2.csv'),'r') as inF:
     lineNum=0
     for line in inF:
       lineNum+=1
       if line[0] != '#':
         splitL = line.split(',')
         smi = splitL[0]
         try:
           m = Chem.MolFromSmiles(smi)
         except:
           m = None

         assert m,'line %d, smiles: %s'%(lineNum,smi)
         useIt=1
         for atom in m.GetAtoms():
           if atom.GetIsAromatic():
             useIt=0
             break
         if useIt:
           tgtVal = float(splitL[col])
           try:
             val = GraphDescriptors.BertzCT(m)
           except:
             val = 666
           assert feq(val,tgtVal,1e-4),'line %d, mol %s (CT calc = %f) should have CT = %f'%(lineNum,smi,val,tgtVal)

  def __testDesc(self,fileN,col,func):
   with open(os.path.join(RDConfig.RDCodeDir,'Chem','test_data',fileN),'r') as inF:
     lineNum=0
     for line in inF:
       lineNum+=1
       if line[0] != '#':
         splitL = line.split(',')
         smi = splitL[0]
         try:
           m = Chem.MolFromSmiles(smi)
         except:
           m = None
         assert m,'line %d, smiles: %s'%(lineNum,smi)
         useIt=1
         if useIt:
           tgtVal = float(splitL[col])
           if not feq(tgtVal,666.0):
             try:
               val = func(m)
             except:
               val = 666
             assert feq(val,tgtVal,1e-4),'line %d, mol %s (calc = %f) should have val = %f'%(lineNum,smi,val,tgtVal)

  def testChi0Long(self):
   """ test calculation of Chi0
     
   """
   col = 2
   self.__testDesc('PP_descrs_regress.csv',col,GraphDescriptors.Chi0)

  def _testChi0Long2(self):
   """ test calculation of Chi0
     
   """
   col = 2
   self.__testDesc('PP_descrs_regress.2.csv',col,GraphDescriptors.Chi0)

  def testHallKierAlphaLong(self):
   """ test calculation of the Hall-Kier Alpha value
     
   """
   col = 3
   self.__testDesc('PP_descrs_regress.csv',col,GraphDescriptors.HallKierAlpha)

  def _testHallKierAlphaLong2(self):
   """ test calculation of the Hall-Kier Alpha value
     
   """
   col = 3
   self.__testDesc('PP_descrs_regress.2.csv',col,GraphDescriptors.HallKierAlpha)

  def testIpc(self):
   """ test calculation of Ipc.

   """
   data = [('CCCCC',1.40564,11.24511),('CCC(C)C',1.37878, 9.65148),('CC(C)(C)C',0.72193,3.60964),('CN(CC)CCC',1.67982,31.91664),('C1CCCCC1',1.71997,34.39946),('CC1CCCCC1',1.68562,47.19725),('Cc1ccccc1',1.68562,47.19725),('CC(C)=C(C)C',1.36096,13.60964),('C#N',1.00000,2.00000),('OC#N',0.91830,2.75489)]
   for smi,res1,res2 in data:
    m = Chem.MolFromSmiles(smi)
    Ipc = GraphDescriptors.Ipc(m, forceDMat=1)
    Ipc_avg = GraphDescriptors.Ipc(m,avg = 1, forceDMat=1)
    assert feq(Ipc_avg,res1,1e-3),'mol %s (Ipc_avg=%f) should have Ipc_avg=%f'%(smi,Ipc_avg,res1)
    assert feq(Ipc,res2,1e-3),'mol %s (Ipc=%f) should have Ipc=%f'%(smi,Ipc,res2)
    Ipc = GraphDescriptors.Ipc(m)
    Ipc_avg = GraphDescriptors.Ipc(m,avg = 1)
    assert feq(Ipc_avg,res1,1e-3),'2nd pass: mol %s (Ipc_avg=%f) should have Ipc_avg=%f'%(smi,Ipc_avg,res1)
    assert feq(Ipc,res2,1e-3),'2nd pass: mol %s (Ipc=%f) should have Ipc=%f'%(smi,Ipc,res2)

  def _testIpcLong(self):
   """ test calculation of Ipc
     
   """
   col = 4
   self.__testDesc('PP_descrs_regress.csv',col,GraphDescriptors.Ipc)

  def _testIpcLong2(self):
   """ test calculation of Ipc
     
   """
   col = 4
   self.__testDesc('PP_descrs_regress.2.csv',col,GraphDescriptors.Ipc)

  def testKappa1(self):
    """ test calculation of the Hall-Kier kappa1 value

     corrected data from Tables 3 and 6 of Rev. Comp. Chem. vol 2, 367-422, (1991)

    """
    data = [('C12CC2C3CC13',2.344),
            ('C1CCC12CC2',3.061),
            ('C1CCCCC1',4.167),
            ('CCCCCC',6.000),
            ('CCC(C)C1CCC(C)CC1',9.091),
            ('CC(C)CC1CCC(C)CC1',9.091),
            ('CC(C)C1CCC(C)CCC1',9.091)
             ]
    for smi,res in data:
      m = Chem.MolFromSmiles(smi)
      kappa = GraphDescriptors.Kappa1(m)
      assert feq(kappa,res,1e-3),'mol %s (kappa1=%f) should have kappa1=%f'%(smi,kappa,res)
      

  def testKappa2(self):
    """ test calculation of the Hall-Kier kappa2 value

     corrected data from Tables 5 and 6 of Rev. Comp. Chem. vol 2, 367-422, (1991)

    """
    data = [
      ('[C+2](C)(C)(C)(C)(C)C',0.667),
      ('[C+](C)(C)(C)(C)(CC)',1.240),
      ('C(C)(C)(C)(CCC)',2.3444),
      ('CC(C)CCCC',4.167),
      ('CCCCCCC',6.000),
      ('CCCCCC',5.000),
      ('CCCCCCC',6.000),
      ('C1CCCC1',1.440),
      ('C1CCCC1C',1.633),
      ('C1CCCCC1',2.222),
      ('C1CCCCCC1',3.061),
      ('CCCCC',4.00),
      ('CC=CCCC',4.740),
      ('C1=CN=CN1',0.884),
      ('c1ccccc1',1.606),
      ('c1cnccc1',1.552),
      ('n1ccncc1',1.500),
      ('CCCCF',3.930),
      ('CCCCCl',4.290),
      ('CCCCBr',4.480),
      ('CCC(C)C1CCC(C)CC1',4.133),
      ('CC(C)CC1CCC(C)CC1',4.133),
      ('CC(C)C1CCC(C)CCC1',4.133)
      ]
    for smi,res in data:
      #print smi
      m = Chem.MolFromSmiles(smi)
      kappa = GraphDescriptors.Kappa2(m)
      assert feq(kappa,res,1e-3),'mol %s (kappa2=%f) should have kappa2=%f'%(smi,kappa,res)
      
  def testKappa3(self):
    """ test calculation of the Hall-Kier kappa3 value

     corrected data from Tables 3 and 6 of Rev. Comp. Chem. vol 2, 367-422, (1991)

    """
    data = [
      ('C[C+](C)(C)(C)C(C)(C)C',2.000),
      ('CCC(C)C(C)(C)(CC)',2.380),
      ('CCC(C)CC(C)CC',4.500),
      ('CC(C)CCC(C)CC',5.878),
      ('CC(C)CCCC(C)C',8.000),
      ('CCC(C)C1CCC(C)CC1',2.500),
      ('CC(C)CC1CCC(C)CC1',3.265),
      ('CC(C)C1CCC(C)CCC1',2.844)
      ]
    for smi,res in data:
      m = Chem.MolFromSmiles(smi)
      kappa = GraphDescriptors.Kappa3(m)
      assert feq(kappa,res,1e-3),'mol %s (kappa3=%f) should have kappa3=%f'%(smi,kappa,res)
      
  def testKappa3Long(self):
   """ test calculation of kappa3
     
   """
   col = 5
   self.__testDesc('PP_descrs_regress.csv',col,GraphDescriptors.Kappa3)

  def _testKappa3Long2(self):
   """ test calculation of kappa3
     
   """
   col = 5
   self.__testDesc('PP_descrs_regress.2.csv',col,GraphDescriptors.Kappa3)

  def _testLabuteASALong(self):
   """ test calculation of Labute's ASA value
     
   """
   col = 6
   self.__testDesc('PP_descrs_regress.csv',col,lambda x:MolSurf.LabuteASA(x,includeHs=1))

  def _testLabuteASALong2(self):
   """ test calculation of Labute's ASA value
     
   """
   col = 6
   self.__testDesc('PP_descrs_regress.2.csv',col,lambda x:MolSurf.LabuteASA(x,includeHs=1))

  def _testTPSAShortNCI(self):
    " Short TPSA test "
    inName = RDConfig.RDDataDir+'/NCI/first_200.tpsa.csv'
    with open(inName,'r') as inF:
      lines = inF.readlines()
    for line in lines:
      if line[0] != '#':
        line.strip()
        smi,ans = line.split(',')
        ans = float(ans)

        mol = Chem.MolFromSmiles(smi)

        calc = MolSurf.TPSA(mol)
        assert feq(calc,ans),'bad TPSA for SMILES %s (%.2f != %.2f)'%(smi,calc,ans)

  def _testTPSALongNCI(self):
    " Long TPSA test "
    fileN = 'tpsa_regr.csv'
    with open(os.path.join(RDConfig.RDCodeDir,'Chem','test_data',fileN),'r') as inF:
      lines = inF.readlines()
    lineNo = 0
    for line in lines:
      lineNo+=1
      if line[0] != '#':
        line.strip()
        smi,ans = line.split(',')
        ans = float(ans)

        try:
          mol = Chem.MolFromSmiles(smi)
        except:
          import traceback
          traceback.print_exc()
          mol = None
        assert mol,"line %d, failed for smiles: %s"%(lineNo,smi)

      
        calc = MolSurf.TPSA(mol)
        assert feq(calc,ans),'line %d: bad TPSA for SMILES %s (%.2f != %.2f)'%(lineNo,smi,calc,ans)

  def testTPSALong(self):
   """ test calculation of TPSA
     
   """
   col = 28
   self.__testDesc('PP_descrs_regress.csv',col,MolSurf.TPSA)

            
  def _testTPSALong2(self):
   """ test calculation of TPSA
     
   """
   col = 28
   self.__testDesc('PP_descrs_regress.2.csv',col,MolSurf.TPSA)

  def _testLipinskiLong(self):
   """ test calculation of Lipinski params
     
   """
   fName = 'PP_descrs_regress.csv'
   # we can't do H Acceptors for these pyridine-containing molecules
   #  because the values will be wrong for EVERY one.
   #col = 29
   #self.__testDesc(fName,col,Lipinski.NumHAcceptors)

   col = 30
   self.__testDesc(fName,col,Lipinski.NumHDonors)

   col = 31
   self.__testDesc(fName,col,Lipinski.NumHeteroatoms)

   col = 32
   self.__testDesc(fName,col,Lipinski.NumRotatableBonds)

  def _testHAcceptorsLong(self):
   """ test calculation of Lipinski params
     
   """
   fName = 'Block_regress.Lip.csv'
   col = 1
   self.__testDesc(fName,col,Lipinski.NumHAcceptors)

  def _testHDonorsLong(self):
   """ test calculation of Lipinski params
     
   """
   fName = 'Block_regress.Lip.csv'
   col = 2
   self.__testDesc(fName,col,Lipinski.NumHDonors)

  def _testHeterosLong(self):
   """ test calculation of Lipinski params
     
   """
   fName = 'Block_regress.Lip.csv'
   col = 3
   self.__testDesc(fName,col,Lipinski.NumHeteroatoms)

  def _testRotBondsLong(self):
   """ test calculation of Lipinski params
     
   """
   fName = 'Block_regress.Lip.csv'
   col = 4
   self.__testDesc(fName,col,Lipinski.NumRotatableBonds)

  def _testLogPLong(self):
   """ test calculation of Lipinski params
     
   """
   fName = 'PP_descrs_regress.csv'
   col = 33
   self.__testDesc(fName,col,lambda x:Crippen.MolLogP(x,includeHs=1))
            
  def _testLogPLong2(self):
   """ test calculation of Lipinski params
     
   """
   fName = 'PP_descrs_regress.2.csv'
   col = 33
   self.__testDesc(fName,col,lambda x:Crippen.MolLogP(x,includeHs=1))

  def _testMOELong(self):
   """ test calculation of MOE-type descriptors
     
   """
   fName = 'PP_descrs_regress.VSA.csv'
   col = 1
   self.__testDesc(fName,col,MolSurf.SMR_VSA1)
   col = 2
   self.__testDesc(fName,col,MolSurf.SMR_VSA10)
   col = 3
   self.__testDesc(fName,col,MolSurf.SMR_VSA2)
   col = 4
   self.__testDesc(fName,col,MolSurf.SMR_VSA3)
   col = 5
   self.__testDesc(fName,col,MolSurf.SMR_VSA4)
   col = 6
   self.__testDesc(fName,col,MolSurf.SMR_VSA5)
   col = 7
   self.__testDesc(fName,col,MolSurf.SMR_VSA6)
   col = 8
   self.__testDesc(fName,col,MolSurf.SMR_VSA7)
   col = 9
   self.__testDesc(fName,col,MolSurf.SMR_VSA8)
   col = 10
   self.__testDesc(fName,col,MolSurf.SMR_VSA9)
   col = 11
   self.__testDesc(fName,col,MolSurf.SlogP_VSA1)
   col = 12
   self.__testDesc(fName,col,MolSurf.SlogP_VSA10)
   col = 13
   self.__testDesc(fName,col,MolSurf.SlogP_VSA11)
   col = 14
   self.__testDesc(fName,col,MolSurf.SlogP_VSA12)



  def _testMOELong2(self):
   """ test calculation of MOE-type descriptors
     
   """
   fName = 'PP_descrs_regress.VSA.2.csv'
   col = 1
   self.__testDesc(fName,col,MolSurf.SMR_VSA1)
   col = 2
   self.__testDesc(fName,col,MolSurf.SMR_VSA10)
   col = 11
   self.__testDesc(fName,col,MolSurf.SlogP_VSA1)
   col = 12
   self.__testDesc(fName,col,MolSurf.SlogP_VSA10)
   col = 13
   self.__testDesc(fName,col,MolSurf.SlogP_VSA11)
   col = 14
   self.__testDesc(fName,col,MolSurf.SlogP_VSA12)

  def testBalabanJ(self):
    """ test calculation of the Balaban J value 

      J values are from Balaban's paper and have had roundoff
      errors and typos corrected.
    """
    data = [# alkanes
      ('CC',1.0),('CCC',1.6330),
      ('CCCC',1.9747),('CC(C)C',2.3238),
      ('CCCCC',2.1906),('CC(C)CC',2.5396),('CC(C)(C)C',3.0237),
      ('CCCCCC',2.3391),('CC(C)CCC',2.6272),('CCC(C)CC',2.7542),('CC(C)(C)CC',3.1685),
      ('CC(C)C(C)C',2.9935),

      # cycloalkanes
      ('C1CCCCC1',2.0000),
      ('C1C(C)CCCC1',2.1229),
      ('C1C(CC)CCCC1',2.1250),
      ('C1C(C)C(C)CCC1',2.2794),
      ('C1C(C)CC(C)CC1',2.2307),
      ('C1C(C)CCC(C)C1',2.1924),
      ('C1C(CCC)CCCC1',2.0779),
      ('C1C(C(C)C)CCCC1',2.2284),
      ('C1C(CC)C(C)CCC1',2.2973),
      ('C1C(CC)CC(C)CC1',2.2317),
      ('C1C(CC)CCC(C)C1',2.1804),
      ('C1C(C)C(C)C(C)CC1',2.4133),
      ('C1C(C)C(C)CC(C)C1',2.3462),
      ('C1C(C)CC(C)CC1(C)',2.3409),
      # aromatics
      ('c1ccccc1',3.0000),
      ('c1c(C)cccc1',3.0215),
      ('c1c(CC)cccc1',2.8321),
      ('c1c(C)c(C)ccc1',3.1349),
      ('c1c(C)cc(C)cc1',3.0777),
      ('c1c(C)ccc(C)c1',3.0325),
      ('c1c(CCC)cccc1',2.6149),
      ('c1c(C(C)C)cccc1',2.8483),
      ('c1c(CC)c(C)ccc1',3.0065),
      ('c1c(CC)cc(C)cc1',2.9369),
      ('c1c(CC)ccc(C)c1',2.8816),
      ('c1c(C)c(C)c(C)cc1',3.2478),
      ('c1c(C)c(C)cc(C)c1',3.1717),
      ('c1c(C)cc(C)cc1(C)',3.1657)
      ]
    for smi,res in data:
      m = Chem.MolFromSmiles(smi)
      j = GraphDescriptors.BalabanJ(m,forceDMat=1)
      assert feq(j,res),'mol %s (J=%f) should have J=%f'%(smi,j,res)
      j = GraphDescriptors.BalabanJ(m)
      assert feq(j,res),'second pass: mol %s (J=%f) should have J=%f'%(smi,j,res)
    
  def _testBalabanJLong(self):
   """ test calculation of the balaban j value
     
   """
   fName = 'PP_descrs_regress.rest.2.csv'
   col = 1
   self.__testDesc(fName,col,GraphDescriptors.BalabanJ)
    
  def _testKappa1Long(self):
   """ test calculation of kappa1
     
   """
   fName = 'PP_descrs_regress.rest.2.csv'
   col = 31
   self.__testDesc(fName,col,GraphDescriptors.Kappa1)
    
  def _testKappa2Long(self):
   """ test calculation of kappa2
     
   """
   fName = 'PP_descrs_regress.rest.2.csv'
   col = 32
   self.__testDesc(fName,col,GraphDescriptors.Kappa2)
    

  def _testChi0Long(self):
   fName = 'PP_descrs_regress.rest.2.csv'
   col = 5
   self.__testDesc(fName,col,GraphDescriptors.Chi0)

  def _testChi1Long(self):
   fName = 'PP_descrs_regress.rest.2.csv'
   col = 8
   self.__testDesc(fName,col,GraphDescriptors.Chi1)



  def _testChi0v(self):
    """ test calculation of Chi0v

    """
    data = [('CCCCCC',4.828),('CCC(C)CC',4.992),('CC(C)CCC',4.992),
            ('CC(C)C(C)C',5.155),('CC(C)(C)CC',5.207),
            ('CCCCCO',4.276),('CCC(O)CC',4.439),('CC(O)(C)CC',4.654),('c1ccccc1O',3.834),
            ('CCCl',2.841),('CCBr',3.671),('CCI',4.242)]
    for smi,res in data:
      m = Chem.MolFromSmiles(smi)
      chi = GraphDescriptors.Chi0v(m)
      assert feq(chi,res,1e-3),'mol %s (Chi0v=%f) should have Chi0V=%f'%(smi,chi,res)

  def _testChi0vLong(self):
   fName = 'PP_descrs_regress.rest.2.csv'
   col = 7
   self.__testDesc(fName,col,GraphDescriptors.Chi0v)

  def testChi1v(self):
    """ test calculation of Chi1v

    """
    data = [('CCCCCC',2.914),('CCC(C)CC',2.808),('CC(C)CCC',2.770),
            ('CC(C)C(C)C',2.643),('CC(C)(C)CC',2.561),
            ('CCCCCO',2.523),('CCC(O)CC',2.489),('CC(O)(C)CC',2.284),('c1ccccc1O',2.134)]
    for smi,res in data:
      m = Chem.MolFromSmiles(smi)
      chi = GraphDescriptors.Chi1v(m)
      assert feq(chi,res,1e-3),'mol %s (Chi1v=%f) should have Chi1V=%f'%(smi,chi,res)

  def _testChi1vLong(self):
   fName = 'PP_descrs_regress.rest.2.csv'
   col = 10
   self.__testDesc(fName,col,GraphDescriptors.Chi1v)

  def testPathCounts(self):
    """ FIX: this should be in some other file

    """
    data = [('CCCCCC',(6,5,4,3,2,1)),
            ('CCC(C)CC',(6,5,5,4,1,0)),
            ('CC(C)CCC',(6,5,5,3,2,0)),
            ('CC(C)C(C)C',(6,5,6,4,0,0)),
            ('CC(C)(C)CC',(6,5,7,3,0,0)),
            ('CCCCCO',(6,5,4,3,2,1)),
            ('CCC(O)CC',(6,5,5,4,1,0)),
            ('CC(O)(C)CC',(6,5,7,3,0,0)),
            ('c1ccccc1O',(7,7,8,8,8,8)),
            ]
    for smi,res in data:
      m = Chem.MolFromSmiles(smi)
      for i in range(1,6):
        cnt = len(Chem.FindAllPathsOfLengthN(m,i,useBonds=1))
        assert cnt==res[i],(smi,i,cnt,res[i],Chem.FindAllPathsOfLengthN(m,i,useBonds=1))
        cnt = len(Chem.FindAllPathsOfLengthN(m,i+1,useBonds=0))
        assert cnt==res[i],(smi,i,cnt,res[i],Chem.FindAllPathsOfLengthN(m,i+1,useBonds=1))
  def testChi2v(self):
    """ test calculation of Chi2v

    """
    data = [('CCCCCC',1.707),('CCC(C)CC',1.922),('CC(C)CCC',2.183),
            ('CC(C)C(C)C',2.488),('CC(C)(C)CC',2.914),
            ('CCCCCO',1.431),('CCC(O)CC',1.470),('CC(O)(C)CC',2.166),('c1ccccc1O',1.336),
            ]
    for smi,res in data:
      m = Chem.MolFromSmiles(smi)
      chi = GraphDescriptors.Chi2v(m)
      assert feq(chi,res,1e-3),'mol %s (Chi2v=%f) should have Chi2V=%f'%(smi,chi,res)
  def _testChi2vLong(self):
   fName = 'PP_descrs_regress.rest.2.csv'
   col = 12
   self.__testDesc(fName,col,GraphDescriptors.Chi2v)

  def testChi3v(self):
    """ test calculation of Chi3v

    """
    data = [('CCCCCC',0.957),('CCC(C)CC',1.394),('CC(C)CCC',0.866),('CC(C)C(C)C',1.333),('CC(C)(C)CC',1.061),
            ('CCCCCO',0.762),('CCC(O)CC',0.943),('CC(O)(C)CC',0.865),('c1ccccc1O',0.756)]
    for smi,res in data:
      m = Chem.MolFromSmiles(smi)
      chi = GraphDescriptors.Chi3v(m)
      assert feq(chi,res,1e-3),'mol %s (Chi3v=%f) should have Chi3V=%f'%(smi,chi,res)
  def _testChi3vLong(self):
   fName = 'PP_descrs_regress.rest.2.csv'
   col = 14
   self.__testDesc(fName,col,GraphDescriptors.Chi3v)

  def testChi4v(self):
    """ test calculation of Chi4v

    """
    data = [('CCCCCC',0.500),('CCC(C)CC',0.289),('CC(C)CCC',0.577),
            ('CC(C)C(C)C',0.000),('CC(C)(C)CC',0.000),
            ('CCCCCO',0.362),('CCC(O)CC',0.289),('CC(O)(C)CC',0.000),('c1ccccc1O',0.428)]
    for smi,res in data:
      m = Chem.MolFromSmiles(smi)
      chi = GraphDescriptors.Chi4v(m)
      assert feq(chi,res,1e-3),'mol %s (Chi4v=%f) should have Chi4V=%f'%(smi,chi,res)

            
  def testChi5v(self):
    """ test calculation of Chi5v

    """
    data = [('CCCCCC',0.250),('CCC(C)CC',0.000),('CC(C)CCC',0.000),
            ('CC(C)C(C)C',0.000),('CC(C)(C)CC',0.000),
            ('CCCCCO',0.112),('CCC(O)CC',0.000),('CC(O)(C)CC',0.000),('c1ccccc1O',0.242)]
    for smi,res in data:
      m = Chem.MolFromSmiles(smi)
      chi = GraphDescriptors.ChiNv_(m,5)
      assert feq(chi,res,1e-3),'mol %s (Chi5v=%f) should have Chi5V=%f'%(smi,chi,res)

  def testChi0n(self):
    """ test calculation of Chi0n

    """
    data = [('CCCCCC',4.828),('CCC(C)CC',4.992),('CC(C)CCC',4.992),
            ('CC(C)C(C)C',5.155),('CC(C)(C)CC',5.207),
            ('CCCCCO',4.276),('CCC(O)CC',4.439),('CC(O)(C)CC',4.654),('c1ccccc1O',3.834),
            ('CCCl',2.085),('CCBr',2.085),('CCI',2.085),]
    for smi,res in data:
      m = Chem.MolFromSmiles(smi)
      chi = GraphDescriptors.Chi0n(m)
      assert feq(chi,res,1e-3),'mol %s (Chi0n=%f) should have Chi0n=%f'%(smi,chi,res)
  def _testChi0nLong(self):
    fName = 'PP_descrs_regress.rest.2.csv'
    col = 6
    self.__testDesc(fName,col,GraphDescriptors.Chi0n)

  def testChi1n(self):
    """ test calculation of Chi1n

    """
    data = [('CCCCCC',2.914),('CCC(C)CC',2.808),('CC(C)CCC',2.770),
            ('CC(C)C(C)C',2.643),('CC(C)(C)CC',2.561),
            ('CCCCCO',2.523),('CCC(O)CC',2.489),('CC(O)(C)CC',2.284),('c1ccccc1O',2.134)]
    for smi,res in data:
      m = Chem.MolFromSmiles(smi)
      chi = GraphDescriptors.Chi1n(m)
      assert feq(chi,res,1e-3),'mol %s (Chi1n=%f) should have Chi1N=%f'%(smi,chi,res)
  def _testChi1nLong(self):
    fName = 'PP_descrs_regress.rest.2.csv'
    col = 9
    self.__testDesc(fName,col,GraphDescriptors.Chi1n)

  def testChi2n(self):
    """ test calculation of Chi2n

    """
    data = [('CCCCCC',1.707),('CCC(C)CC',1.922),('CC(C)CCC',2.183),
            ('CC(C)C(C)C',2.488),('CC(C)(C)CC',2.914),
            ('CCCCCO',1.431),('CCC(O)CC',1.470),('CC(O)(C)CC',2.166),('c1ccccc1O',1.336)]
    for smi,res in data:
      m = Chem.MolFromSmiles(smi)
      chi = GraphDescriptors.Chi2n(m)
      assert feq(chi,res,1e-3),'mol %s (Chi2n=%f) should have Chi2N=%f'%(smi,chi,res)
  def _testChi2nLong(self):
    fName = 'PP_descrs_regress.rest.2.csv'
    col = 11
    self.__testDesc(fName,col,GraphDescriptors.Chi2n)


  def testChi3n(self):
    """ test calculation of Chi3n

    """
    data = [('CCCCCC',0.957),('CCC(C)CC',1.394),('CC(C)CCC',0.866),('CC(C)C(C)C',1.333),('CC(C)(C)CC',1.061),
            ('CCCCCO',0.762),('CCC(O)CC',0.943),('CC(O)(C)CC',0.865),('c1ccccc1O',0.756)]
    for smi,res in data:
      m = Chem.MolFromSmiles(smi)
      chi = GraphDescriptors.Chi3n(m)
      assert feq(chi,res,1e-3),'mol %s (Chi3n=%f) should have Chi3N=%f'%(smi,chi,res)
  def _testChi3nLong(self):
    fName = 'PP_descrs_regress.rest.2.csv'
    col = 13
    self.__testDesc(fName,col,GraphDescriptors.Chi3n)


  def testChi4n(self):
    """ test calculation of Chi4n

    """
    data = [('CCCCCC',0.500),('CCC(C)CC',0.289),('CC(C)CCC',0.577),
            ('CC(C)C(C)C',0.000),('CC(C)(C)CC',0.000),
            ('CCCCCO',0.362),('CCC(O)CC',0.289),('CC(O)(C)CC',0.000),('c1ccccc1O',0.428)]
    for smi,res in data:
      m = Chem.MolFromSmiles(smi)
      chi = GraphDescriptors.Chi4n(m)
      assert feq(chi,res,1e-3),'mol %s (Chi4n=%f) should have Chi4N=%f'%(smi,chi,res)

  def testIssue125(self):
    """ test an issue with calculating BalabanJ
    """
    smi = 'O=C(OC)C1=C(C)NC(C)=C(C(OC)=O)C1C2=CC=CC=C2[N+]([O-])=O'
    m1 = Chem.MolFromSmiles(smi)
    m2 = Chem.MolFromSmiles(smi)
    Chem.MolToSmiles(m1)
    j1=GraphDescriptors.BalabanJ(m1)
    j2=GraphDescriptors.BalabanJ(m2)
    assert feq(j1,j2)

  def testOrderDepend(self):
    """ test order dependence of some descriptors:
    """
    data = [('C=CC=C',21.01955,2.73205),
            ('O=CC=O',25.01955,2.73205),
            ('FCC(=O)CF',46.7548875,2.98816),
            ('O=C1C=CC(=O)C=C1',148.705216,2.8265),
            ('C12C(F)=C(O)C(F)C1C(F)=C(O)C(F)2',315.250442,2.4509),
            ('C12CC=CCC1C(=O)C3CC=CCC3C(=O)2',321.539522,1.95986)]

    for smi,CT,bal in data:
      m = Chem.MolFromSmiles(smi)
      newBal = GraphDescriptors.BalabanJ(m, forceDMat = 1)
      assert feq(newBal,bal,1e-4),'mol %s %f!=%f'%(smi,newBal,bal)
      m = Chem.MolFromSmiles(smi)
      newCT = GraphDescriptors.BertzCT(m, forceDMat = 1)
      assert feq(newCT,CT,1e-4),'mol %s (CT calc = %f) should have CT = %f'%(smi,newCT,CT)
      m = Chem.MolFromSmiles(smi)
      newCT = GraphDescriptors.BertzCT(m, forceDMat = 1)
      assert feq(newCT,CT,1e-4),'mol %s (CT calc = %f) should have CT = %f'%(smi,newCT,CT)
      newBal = GraphDescriptors.BalabanJ(m, forceDMat = 1)
      assert feq(newBal,bal,1e-4),'mol %s %f!=%f'%(smi,newBal,bal)

      m = Chem.MolFromSmiles(smi)
      newBal = GraphDescriptors.BalabanJ(m, forceDMat = 1)
      assert feq(newBal,bal,1e-4),'mol %s %f!=%f'%(smi,newBal,bal)
      newCT = GraphDescriptors.BertzCT(m, forceDMat = 1)
      assert feq(newCT,CT,1e-4),'mol %s (CT calc = %f) should have CT = %f'%(smi,newCT,CT)

if __name__ == '__main__':
  import sys,getopt,re
  doLong=0
  if len(sys.argv) >1:
    args,extras=getopt.getopt(sys.argv[1:],'l')
    for arg,val in args:
      if arg=='-l':
        doLong=1
      sys.argv.remove('-l')
  if doLong:
    for methName in dir(TestCase):
      if re.match('_test',methName):
        newName = re.sub('_test','test',methName)
        exec('TestCase.%s = TestCase.%s'%(newName,methName))
        
  unittest.main()

