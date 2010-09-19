# $Id$
#
#  Copyright (C) 2007-2010 Greg Landrum
#
#   @@ All Rights Reserved  @@
#
""" General descriptor testing code

"""
from rdkit import RDConfig
import unittest,os.path
from rdkit import Chem
from rdkit.Chem import Descriptors,AvailDescriptors

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
      actual = Descriptors.MolecularFormula(mol)
      self.failUnlessEqual(actual,expected)
      
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
