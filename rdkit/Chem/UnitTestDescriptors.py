# $Id$
#
#  Copyright (C) 2007 Greg Landrum
#
#   @@ All Rights Reserved  @@
#
""" General descriptor testing code

"""
from rdkit import RDConfig
import unittest,os.path
from rdkit import Chem
from rdkit.Chem import AvailDescriptors

def feq(n1,n2,tol=1e-4):
  return abs(n1-n2)<=tol

class TestCase(unittest.TestCase):
  def testBadAtomHandling(self):
    smis = ('CC[Pu]','CC[*]')
    for smi in smis:
      m = Chem.MolFromSmiles(smi)
      self.failUnless(m)
      for nm,fn in AvailDescriptors.descList:
        try:
          v = fn(m)
        except:
          import traceback
          traceback.print_exc()
          self.failUnless(0,'SMILES: %s'%smi)
            
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
