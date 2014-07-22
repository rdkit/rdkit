# $Id$
#
#  Copyright (C) 2002-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" unit testing code for surface calculations

FIX: add tests for LabuteASA

"""
from __future__ import print_function
from rdkit import RDConfig
import unittest,os
from rdkit.six.moves import cPickle
from rdkit import Chem
from rdkit.Chem import MolSurf
import os.path

def feq(n1,n2,tol=1e-4):
  return abs(n1-n2)<=tol

class TestCase(unittest.TestCase):
  def setUp(self):
    if doLong:
      print('\n%s: '%self.shortDescription(),end='')

  def testTPSAShort(self):
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
        
  def _testTPSALong(self):
    " Longer TPSA test "
    #inName = RDConfig.RDDataDir+'/NCI/first_5k.tpsa.csv'
    inName = os.path.join(RDConfig.RDCodeDir,'Chem','test_data','NCI_5K_TPSA.csv')
    with open(inName,'r') as inF:
      lines = inF.readlines()
    lineNo = 0
    for line in lines:
      lineNo += 1
      if line[0] != '#':
        line.strip()
        smi,ans = line.split(',')
        ans = float(ans)
        try:
          mol = Chem.MolFromSmiles(smi)
        except:
          mol = None
        if not mol:
          print('molecule construction failed on line %d'%lineNo)
        else:
          ok = 1
          try:
            calc = MolSurf.TPSA(mol)
          except:
            ok=0
          assert ok,'Line %d: TPSA Calculation failed for SMILES %s'%(lineNo,smi)
          assert feq(calc,ans),'Line %d: bad TPSA for SMILES %s (%.2f != %.2f)'%(lineNo,smi,calc,ans)

  def testHsAndTPSA(self):
    """
     testing the impact of Hs in the graph on PSA calculations
     This was sf.net issue 1969745 
    """
    mol = Chem.MolFromSmiles('c1c[nH]cc1')
    molH = Chem.AddHs(mol)
    psa = MolSurf.TPSA(mol)
    psaH = MolSurf.TPSA(molH)

    if(psa!=psaH):
      psac = MolSurf.rdMolDescriptors._CalcTPSAContribs(mol)
      psaHc = MolSurf.rdMolDescriptors._CalcTPSAContribs(molH)
      for i,v in enumerate(psac):
        print('\t',i,'\t',v,'\t',psaHc[i])
      while i<len(psaHc):
        print('\t\t\t',psaHc[i])
        i+=1
    self.assertEqual(psa,psaH)

    inName = RDConfig.RDDataDir+'/NCI/first_200.tpsa.csv'
    with open(inName,'r') as inF:
      lines = inF.readlines()
    for line in lines:
      if line[0] != '#':
        line.strip()
        smi,ans = line.split(',')
        ans = float(ans)
        mol = Chem.MolFromSmiles(smi)
        mol = Chem.AddHs(mol)
        calc = MolSurf.TPSA(mol)
        self.assertTrue(feq(calc,ans),'bad TPSA for SMILES %s (%.2f != %.2f)'%(smi,calc,ans))
    if doLong:
      inName = os.path.join(RDConfig.RDCodeDir,'Chem','test_data','NCI_5K_TPSA.csv')
      with open(inName,'r') as inF:
        lines = inF.readlines()
      for line in lines:
        if line[0] != '#':
          line.strip()
          smi,ans = line.split(',')
          ans = float(ans)
          mol = Chem.MolFromSmiles(smi)
          mol = Chem.AddHs(mol)
          calc = MolSurf.TPSA(mol)
          self.assertTrue(feq(calc,ans),'bad TPSA for SMILES %s (%.2f != %.2f)'%(smi,calc,ans))
      
    
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

