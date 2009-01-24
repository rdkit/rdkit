# $Id$
#
#  Copyright (C) 2003-2006  Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
"""basic unit testing code for the wrapper of the SMARTS matcher

"""
from pyRDKit import RDConfig
import unittest
import os.path
from pyRDKit import Chem


class TestCase(unittest.TestCase):
  def setUp(self):
    #print '\n%s: '%self.shortDescription(),
    fName = os.path.join(RDConfig.RDCodeDir,'Chem','tests','quinone.mol')
    self.m = Chem.MolFromMolFile(fName)
    assert self.m.GetNumAtoms()==8,'bad nAtoms'


  def testMatch(self):
    " testing smarts match "
    p =  Chem.MolFromSmarts('CC(=O)C')
    matches = self.m.GetSubstructMatches(p)
    assert len(matches)==2,'bad UMapList: %s'%(str(res))
    for match in matches:
      assert len(match)==4,'bad match: %s'%(str(match))


  def testOrder(self):
    " testing atom order in smarts match "
    p =  Chem.MolFromSmarts('CC(=O)C')
    matches = self.m.GetSubstructMatches(p)
    m = matches[0]
    atomList = map(lambda x,y=self.m:y.GetAtomWithIdx(x).GetSymbol(),m)
    assert atomList==['C','C','O','C'],'bad atom ordering: %s'%str(atomList)

if __name__ == '__main__':
  unittest.main()
