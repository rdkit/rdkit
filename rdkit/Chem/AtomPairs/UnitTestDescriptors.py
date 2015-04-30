# $Id$
#
#  Copyright (C) 2007  greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
from __future__ import print_function
import unittest,os,gzip
from rdkit.six.moves import cPickle  #@UnresolvedImport #pylint: disable=F0401
from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem.AtomPairs import Pairs,Torsions,Utils

class TestCase(unittest.TestCase):
  def setUp(self):
    self.testDataPath=os.path.join(RDConfig.RDCodeDir,'Chem','AtomPairs','test_data')
    inF = gzip.open(os.path.join(self.testDataPath,'mols1000.pkl.gz'),'rb')
    self.mols=cPickle.load(inF, encoding='bytes')

  def testPairsRegression(self):
    inF = gzip.open(os.path.join(self.testDataPath,'mols1000.aps.pkl.gz'),'rb')
    atomPairs = cPickle.load(inF, encoding='bytes')
    for i,m in enumerate(self.mols):
      ap = Pairs.GetAtomPairFingerprint(m)
      #if ap!=atomPairs[i]:
      #  print Chem.MolToSmiles(m)
      #  pd=ap.GetNonzeroElements()
      #  rd=atomPairs[i].GetNonzeroElements()
      #  for k,v in pd.iteritems():
      #    if rd.has_key(k):
      #      if rd[k]!=v: print '>>>1',k,v,rd[k]
      #    else:
      #      print '>>>2',k,v
      #  for k,v in rd.iteritems():
      #    if pd.has_key(k):
      #      if pd[k]!=v: print '>>>3',k,v,pd[k]
      #    else:
      #      print '>>>4',k,v
      self.assertTrue(ap==atomPairs[i])
      self.assertTrue(ap!=atomPairs[i-1])

  def testTorsionsRegression(self):
    inF = gzip.open(os.path.join(self.testDataPath,'mols1000.tts.pkl.gz'),'rb')
    torsions = cPickle.load(inF, encoding='bytes')
    for i,m in enumerate(self.mols):
      tt = Torsions.GetTopologicalTorsionFingerprintAsIntVect(m)
      if tt!=torsions[i]:
        print(Chem.MolToSmiles(m))
        pd=tt.GetNonzeroElements()
        rd=torsions[i].GetNonzeroElements()
        for k,v in pd.iteritems():
          if rd.has_key(k):
            if rd[k]!=v: print('>>>1',k,v,rd[k])
          else:
            print('>>>2',k,v)
        for k,v in rd.iteritems():
          if pd.has_key(k):
            if pd[k]!=v: print('>>>3',k,v,pd[k])
          else:
            print('>>>4',k,v)
       
      self.assertTrue(tt==torsions[i])
      self.assertTrue(tt!=torsions[i-1])

  def testGithub334(self):
    m1 = Chem.MolFromSmiles('N#C')
    self.assertEqual(Utils.NumPiElectrons(m1.GetAtomWithIdx(0)),2)
    self.assertEqual(Utils.NumPiElectrons(m1.GetAtomWithIdx(1)),2)

    m1 = Chem.MolFromSmiles('N#[CH]')
    self.assertEqual(Utils.NumPiElectrons(m1.GetAtomWithIdx(0)),2)
    self.assertEqual(Utils.NumPiElectrons(m1.GetAtomWithIdx(1)),2)

if __name__ == '__main__':
  unittest.main()

