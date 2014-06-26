# $Id$
#
#  Copyright (C) 2004-2006  Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" unit testing code for MOE-type descriptors with EStates

"""
from __future__ import print_function
from rdkit import RDConfig
import unittest,os
from rdkit import Chem
from rdkit.Chem.EState import EState_VSA
import os.path

def feq(n1,n2,tol=1e-4):
  return abs(n1-n2)<=tol

class TestCase(unittest.TestCase):
  def setUp(self):
    if doLong:
      print('\n%s: '%self.shortDescription(), end='')
  def test1(self):
    inName = os.path.join(RDConfig.RDCodeDir,'Chem','EState','test_data',
                          'EState_VSA.csv')
    with open(inName,'r') as inF:
      inL = inF.readline()
    names = [x.strip() for x in inL.split(',')[1:]]
    suppl = Chem.SmilesMolSupplier(inName,delimiter=',',nameColumn=-1)
    for mol in suppl:
      self.assertTrue(mol)
      smi = Chem.MolToSmiles(mol)
      for name in names:
        prop = float(mol.GetProp(name))
        func = getattr(EState_VSA,name)
        v = func(mol)
        self.assertTrue(feq(v,prop),'%s: %.4f!=%.4f'%(smi,v,prop))

        
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

