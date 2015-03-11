# $Id$
#
#  Copyright (C) 2007-2010 Greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" General descriptor testing code

"""
from __future__ import print_function
from rdkit import RDConfig
import unittest,os.path
import io
from rdkit.six.moves import cPickle
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
import numpy as np

def feq(n1,n2,tol=1e-4):
  return abs(n1-n2)<=tol

class TestCase(unittest.TestCase):
  def testBadAtomHandling(self):
    smis = ('CC[Pu]','CC[*]')
    for smi in smis:
      m = Chem.MolFromSmiles(smi)
      self.assertTrue(m)
      for nm,fn in Descriptors._descList:
        try:
          v = fn(m)
        except:
          import traceback
          traceback.print_exc()
          self.assertTrue(0,'SMILES: %s'%smi)

  def testMolFormula(self):
    for (smiles, expected) in (  ("[NH4+]", "H4N+"),
                                 ("c1ccccc1", "C6H6"),
                                 ("C1CCCCC1", "C6H12"),
                                 ("c1ccccc1O", "C6H6O"),
                                 ("C1CCCCC1O", "C6H12O"),
                                 ("C1CCCCC1=O", "C6H10O"),
                                 ("N[Na]", "H2NNa"),
                                 ("[C-][C-]", "C2-2"),
                                 ("[H]", "H"),
                                 ("[H-1]", "H-"),
                                 ("[H-1]", "H-"),
                                 ("[CH2]", "CH2"),
                                 ("[He-2]", "He-2"),
                                 ("[U+3]", "U+3"),
                                 ):
      mol = Chem.MolFromSmiles(smiles)
      actual = AllChem.CalcMolFormula(mol)
      self.assertEqual(actual,expected)
  def testMQNDetails(self):
    refFile = os.path.join(RDConfig.RDCodeDir,'Chem','test_data','MQNs_regress.pkl')
    with open(refFile,'r') as intf:
      buf = intf.read().replace('\r\n', '\n').encode('utf-8')
      intf.close()
    with io.BytesIO(buf) as inf:
      pkl = inf.read()
    refData  = cPickle.loads(pkl,encoding='bytes')
    fn = os.path.join(RDConfig.RDCodeDir,'Chem','test_data','aromat_regress.txt')
    ms = [x for x in Chem.SmilesMolSupplier(fn,delimiter='\t')]
    for i,m in enumerate(ms):
      mqns = rdMolDescriptors.MQNs_(m) 
      if mqns!=refData[i][1]:
        indices=[(j,x,y) for j,x,y in zip(range(len(mqns)),mqns,refData[i][1]) if x!=y]
        print(Chem.MolToSmiles(m),indices)
      self.assertEqual(mqns,refData[i][1])
  def testMQN(self):
    tgt = np.array([42917,   274,   870,   621,   135,  1582,    29,  3147,  5463,
        6999,   470,    81, 19055,  4424,   309, 24061, 17820,     1,
        8314, 24146, 16076,  5560,  4262,   646,   746, 13725,  5430,
        2629,   362, 24211, 15939,   292,    41,    20,  1852,  5642,
          31,     9,     1,     2,  3060,  1750])
    fn = os.path.join(RDConfig.RDCodeDir,'Chem','test_data','aromat_regress.txt')
    ms = [x for x in Chem.SmilesMolSupplier(fn,delimiter='\t')]
    vs = np.zeros((42,),np.int32)
    for m in ms:
      vs += rdMolDescriptors.MQNs_(m)
    self.assertFalse(False in (vs==tgt))
    
          
      
# - - - - - 
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
