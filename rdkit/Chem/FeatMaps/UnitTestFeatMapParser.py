# $Id$
#
#  Copyright (C) 2006  greg Landrum 
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
from rdkit import RDConfig
import unittest,sys,os,math
from rdkit import Chem
from rdkit.Chem.FeatMaps import FeatMaps,FeatMapParser
from rdkit.Geometry import Point3D

def feq(n1,n2,tol=1e-5):
  return abs(n1-n2)<=tol

class TestCase(unittest.TestCase):
  def setUp(self):
    pass
  
  def test1Basics(self):
    txt="""

ScoreMode=Best
DirScoreMode=DotFullRange

BeginParams
  family=Aromatic radius=2.5 width=1.0 profile=Triangle
  family=Acceptor radius=1.5
EndParams

# optional
BeginPoints
  family=Acceptor pos=(1.0, 0.0, 5.0) weight=1.25 dir=(1, 1, 0)
  family=Aromatic pos=(0.0,1.0,0.0) weight=2.0 dir=(0,0,1) dir=(0,0,-1)
  family=Acceptor pos=(1.0,1.0,2.0) weight=1.25
EndPoints

"""
    p = FeatMapParser.FeatMapParser()
    p.SetData(txt)
    fm = p.Parse()
    self.assertTrue(fm.scoreMode==FeatMaps.FeatMapScoreMode.Best)
    self.assertTrue(fm.dirScoreMode==FeatMaps.FeatDirScoreMode.DotFullRange)
    self.assertTrue(fm.GetNumFeatures()==3)


    feats = fm.GetFeatures()
    self.assertTrue(feq(feats[0].weight,1.25))
    self.assertTrue(feq(feats[1].weight,2.0))
    self.assertTrue(feq(feats[2].weight,1.25))
    
    self.assertTrue(len(feats[0].featDirs)==1)
    self.assertTrue(len(feats[1].featDirs)==2)
    self.assertTrue(len(feats[2].featDirs)==0)
    
    fams = [x.GetFamily() for x in feats]
    self.assertTrue(fams==['Acceptor','Aromatic','Acceptor'])
    
if __name__ == '__main__':
  unittest.main()

