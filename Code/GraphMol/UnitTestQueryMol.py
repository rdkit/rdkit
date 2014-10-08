# Copyright (C) 2001 greg Landrum and Rational Discovery LLC

"""basic unit testing code for query mols

"""
from __future__ import print_function

from rdkit import RDConfig
import unittest,os,sys

class TestCase(unittest.TestCase):
  def setUp(self):
    print('\n%s: '%self.shortDescription(),end='')
    # decipher the name of the executable
    if(sys.platform == 'win32'):
      exe = 'QueryMolTest___Win32_Debug/QueryMolTest.exe'
    else:  
      exe = 'querytest.exe'
    # update to use the full path  
    self.exe = '%s/Code/GraphMol/%s'%(RDConfig.RDBaseDir,exe)

  def test1(self):
    """ the basic test """
    res = os.system(self.exe)
    assert res == 0, 'test failed'

if __name__ == '__main__':
  unittest.main()

