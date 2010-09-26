# $Id$
#
#  Copyright (C) 2004-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
from rdkit import RDConfig
import unittest,sys,os,cPickle
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
import EmbedLib
import gzip

from rdkit import DistanceGeometry as DG
from rdkit import Geometry
import Pharmacophore

def feq(n1,n2,tol=1e-5):
  return abs(n1-n2)<=tol

class TestCase(unittest.TestCase):
  def setUp(self):
    self.fdefBlock = \
                   """DefineFeature HAcceptor1 [N,O;H0]
                      Family HBondAcceptor
                      Weights 1.0
                   EndFeature
                   DefineFeature HDonor1 [N,O;!H0]
                      Family HBondDonor
                      Weights 1.0
                   EndFeature
                   DefineFeature Aromatic1 c1ccccc1
                      Family Aromatic
                      Weights 1.0,1.0,1.0,1.0,1.0,1.0
                   EndFeature\n"""

    self.featFactory = ChemicalFeatures.BuildFeatureFactoryFromString(self.fdefBlock)
    self.feats = [ChemicalFeatures.FreeChemicalFeature('HBondAcceptor', 'HAcceptor1',
                                           Geometry.Point3D(0.0, 0.0, 0.0)),
                  ChemicalFeatures.FreeChemicalFeature('HBondDonor', 'HDonor1',
                                           Geometry.Point3D(2.65, 0.0, 0.0)),
                  ChemicalFeatures.FreeChemicalFeature('Aromatic', 'Aromatic1',
                                           Geometry.Point3D(5.12, 0.908, 0.0)),
                  ]
    self.pcophore=Pharmacophore.Pharmacophore(self.feats)

  def test1Basics(self):
    pcophore = self.pcophore
    self.failUnless(len(pcophore.getFeatures())==3)
    self.failUnless(pcophore.getFeature(0))
    self.failUnless(pcophore.getFeature(1))
    self.failUnless(pcophore.getFeature(2))
    self.failUnlessRaises(IndexError,pcophore.getFeature,3)
    
  def test2BoundSetting(self):
    pcophore = self.pcophore

    pcophore.setUpperBound(0,1,3.0)
    self.failUnless(feq(pcophore.getUpperBound(0,1),3.0))
    self.failUnless(feq(pcophore.getUpperBound(1,0),3.0))
    pcophore.setUpperBound(1,0,5.0)
    self.failUnless(feq(pcophore.getUpperBound(0,1),5.0))
    self.failUnless(feq(pcophore.getUpperBound(1,0),5.0))
    self.failUnlessRaises(IndexError,pcophore.setUpperBound,0,3,2.0)
    self.failUnlessRaises(ValueError,pcophore.setUpperBound,0,3,2.0,checkBounds=True)
    self.failUnlessRaises(IndexError,pcophore.setUpperBound,3,0,2.0)
    self.failUnlessRaises(ValueError,pcophore.setUpperBound,3,0,2.0,checkBounds=True)


    pcophore.setLowerBound(0,1,2.0)
    self.failUnless(feq(pcophore.getLowerBound(0,1),2.0))
    self.failUnless(feq(pcophore.getLowerBound(1,0),2.0))
    pcophore.setLowerBound(1,0,3.0)
    self.failUnless(feq(pcophore.getLowerBound(0,1),3.0))
    self.failUnless(feq(pcophore.getLowerBound(1,0),3.0))
    self.failUnlessRaises(IndexError,pcophore.setLowerBound,0,3,2.0)
    self.failUnlessRaises(ValueError,pcophore.setLowerBound,0,3,2.0,checkBounds=True)
    self.failUnlessRaises(IndexError,pcophore.setLowerBound,3,0,2.0)
    self.failUnlessRaises(ValueError,pcophore.setLowerBound,3,0,2.0,checkBounds=True)

  def test3Bound2DSetting(self):
    pcophore = self.pcophore

    pcophore.setUpperBound2D(0,1,3)
    self.failUnless(pcophore.getUpperBound2D(0,1)==3)
    self.failUnless(pcophore.getUpperBound2D(1,0)==3)
    pcophore.setUpperBound2D(1,0,5)
    self.failUnless(pcophore.getUpperBound2D(0,1)==5)
    self.failUnless(pcophore.getUpperBound2D(1,0)==5)
    self.failUnlessRaises(IndexError,pcophore.setUpperBound2D,0,3,2)
    self.failUnlessRaises(ValueError,pcophore.setUpperBound2D,0,3,2,checkBounds=True)
    self.failUnlessRaises(IndexError,pcophore.setUpperBound2D,3,0,2)
    self.failUnlessRaises(ValueError,pcophore.setUpperBound2D,3,0,2,checkBounds=True)

    pcophore.setLowerBound2D(0,1,3)
    self.failUnless(pcophore.getLowerBound2D(0,1)==3)
    self.failUnless(pcophore.getLowerBound2D(1,0)==3)
    pcophore.setLowerBound2D(1,0,5)
    self.failUnless(pcophore.getLowerBound2D(0,1)==5)
    self.failUnless(pcophore.getLowerBound2D(1,0)==5)
    self.failUnlessRaises(IndexError,pcophore.setLowerBound2D,0,3,2)
    self.failUnlessRaises(ValueError,pcophore.setLowerBound2D,0,3,2,checkBounds=True)
    self.failUnlessRaises(IndexError,pcophore.setLowerBound2D,3,0,2)
    self.failUnlessRaises(ValueError,pcophore.setLowerBound2D,3,0,2,checkBounds=True)


    
if __name__ == '__main__':
  unittest.main()

