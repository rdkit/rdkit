# Copyright (C) 2001 greg Landrum and Rational Discovery LLC

"""basic unit testing code for the substructure matching

"""
import RDConfig
import unittest,os,sys

class TestCase(unittest.TestCase):
  def setUp(self):
    #print '\n%s: '%self.shortDescription(),
    # decipher the name of the executable
    if(sys.platform == 'win32'):
      exe = 'SubstructTest2___Win32_Debug/SubstructTest2.exe'
    else:  
      exe = 'test2'
    # update to use the full path
    self.basePath = '%s/Code/GraphMol/Substruct'%(RDConfig.RDBaseDir)
    self.exe = '%s/%s'%(self.basePath,exe)
    
    self.passStr = 'Got a match: \n'
    self.failStr = 'No match\n'

  def testAtomListPass(self):
    """ testing atom list matches which should pass """
    smis = ['C1=CC=CC=C1','C1C=CC=CC=1',
            'N1=CC=CC=C1','C1=NC=CC=C1','C1=CN=CC=C1',
            'C1=CC=NC=C1','C1=CC=CN=C1','C1=CC=CC=N1',
            'P1=CC=CC=C1','C1=PC=CC=C1','C1=CP=CC=C1',
            'C1=CC=PC=C1','C1=CC=CP=C1','C1=CC=CC=P1',
            'C1=C(C)C=CC=C1','C1C=C(CC)C=C(C)C=1',
            ]
    cdxFile = '%s/list-query.cdxml'%(self.basePath)
    for smi in smis:
      p = os.popen('%s %s "%s"'%(self.exe,cdxFile,smi),'r');
      l = p.read()
      assert l.find(self.passStr)!=-1,'no match for (%s): %s'%(smi,l)

  def testAtomListFail(self):
    """ testing atom list matches which should fail """
    smis = ['C1CC=CC=C1','c1ccccc1',
            'O1=CC=CC=C1','C1=NC=CN=C1','C1=CP=CO=C1',
            ]
    cdxFile = '%s/list-query.cdxml'%(self.basePath)
    for smi in smis:
      p = os.popen('%s %s "%s"'%(self.exe,cdxFile,smi),'r');
      l = p.read()
      assert l.find(self.failStr)!=-1,'invalid match for (%s): %s'%(smi,l)

  def testBondListPass(self):
    """ testing bond list matches which should pass """
    smis = ['CCCC=C','C=CCCC','CC=CC=C','C=CC=CC',
            'C1CCC=CC1','C1CC=CC=C1'
            'C=CC=CCCCO','CC=C(C=C)COC',
            ]
    cdxFile = '%s/bond-query.cdxml'%(self.basePath)
    for smi in smis:
      p = os.popen('%s %s "%s"'%(self.exe,cdxFile,smi),'r');
      l = p.read()
      assert l.find(self.passStr)!=-1,'no match for (%s): %s'%(smi,l)

  def testBondListFail(self):
    """ testing bond list matches which should fail """
    smis = ['CCCCC','C=COCC',
            ]
    cdxFile = '%s/bond-query.cdxml'%(self.basePath)
    for smi in smis:
      p = os.popen('%s %s "%s"'%(self.exe,cdxFile,smi),'r');
      l = p.read()
      assert l.find(self.failStr)!=-1,'invalid match for (%s): %s'%(smi,l)

    
if __name__ == '__main__':
  unittest.main()

