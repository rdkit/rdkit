# $Id$
#
#  Copyright (C) 2007  greg Landrum
#
#   @@ All Rights Reserved  @@
#
import unittest,cPickle,os,gzip
import Chem
import RDConfig
from Chem.AtomPairs import Pairs,Torsions

class TestCase(unittest.TestCase):
  def setUp(self):
    self.testDataPath=os.path.join(RDConfig.RDCodeDir,'Chem','AtomPairs','test_data')
    inF = gzip.open(os.path.join(self.testDataPath,'mols1000.pkl.gz'),'rb')
    self.mols=cPickle.load(inF)

  def testPairsRegression(self):
    inF = gzip.open(os.path.join(self.testDataPath,'mols1000.aps.pkl.gz'),'rb')
    atomPairs = cPickle.load(inF)
    for i,m in enumerate(self.mols):
      ap = Pairs.GetAtomPairFingerprintAsIntVect(m)
      self.failUnless(ap==atomPairs[i][1])

  def testTorsionsRegression(self):
    inF = gzip.open(os.path.join(self.testDataPath,'mols1000.tts.pkl.gz'),'rb')
    torsions = cPickle.load(inF)
    for i,m in enumerate(self.mols):
      tt = Torsions.GetTopologicalTorsionFingerprintAsIntVect(m)
      self.failUnless(tt==torsions[i][1])

if __name__ == '__main__':
  unittest.main()

