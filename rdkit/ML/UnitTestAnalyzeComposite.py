# $Id$
#
#  Copyright (C) 2004-2006  Greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""unit testing code for the AnalyzeComposite functionality

"""
from rdkit import RDConfig
import unittest,os
from rdkit.ML import AnalyzeComposite
import cPickle as pickle

def feq(a,b,tol=1e-4):
  if abs(a-b)>tol: return 0
  else: return 1

class TestCase(unittest.TestCase):
  def setUp(self):
    #print '\n%s: '%self.shortDescription(),
    self.baseDir = os.path.join(RDConfig.RDCodeDir,'ML','test_data')
  def test1_Issue163(self):
    name1 = os.path.join(self.baseDir,'humanoral.1.pkl')
    try:
      c1 = pickle.load(open(name1,'rb'))
    except:
      c1 = None
    self.failUnless(c1)
    name2 = os.path.join(self.baseDir,'humanoral.2.pkl')
    try:
      c2 = pickle.load(open(name2,'rb'))
    except:
      c2 = None
    self.failUnless(c2)

    try:
      res = AnalyzeComposite.ProcessIt([c1,c2],verbose=-1)
    except:
      import traceback
      traceback.print_exc()
      ok=0
    else:
      ok=1
    self.failUnless(ok)

    self.failUnless(res[0][0]=='BALABANJ')
    self.failUnless(res[1][0]=='BERTZCT')
    self.failUnless(res[-1][0]=='FR_ALLYLIC_OXID')
    for entry in res:
      self.failUnless(len(entry)==5)
    
if __name__ == '__main__':
  unittest.main()

