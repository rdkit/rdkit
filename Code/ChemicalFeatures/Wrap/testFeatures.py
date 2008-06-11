# $Id$
#
# Copyright (C) 2006-2008 Greg Landrum
#
#  @@ All Rights Reserved @@
#
import RDConfig
import os
import Chem
from Chem import ChemicalFeatures
import unittest
import cPickle
from Geometry import rdGeometry as geom

def feq(v1,v2,tol2=1e-4):
  return abs(v1-v2)<=tol2

def lstFeq(l1, l2, tol=1.e-4):
  if (len(l1) != len(l2)):
    return 0
  for i in range(len(l1)):
    if not feq(l1[i], l2[i], tol):
      return 0
  return 1

def ptFeq(pt1, pt2, tol=0.0001):
  dist = pt1.Distance(pt2)
  return feq(dist, 0.0, tol)


class TestCase(unittest.TestCase):
    def setUp(self):
      pass

    def testBasic(self):
      ffeat = ChemicalFeatures.FreeChemicalFeature()
      ffeat.SetFamily("HBondDonor")
      self.failUnless(ffeat.GetFamily() == "HBondDonor")
      ffeat.SetPos(geom.Point3D(1.0, 2.0, 3.0))
      pos = ffeat.GetPos()
      self.failUnless(ptFeq(pos, geom.Point3D(1.0, 2.0, 3.0)))
      
      ffeat.SetType("HBondDonor1")
      self.failUnless(ffeat.GetType() == "HBondDonor1")

      ffeat = ChemicalFeatures.FreeChemicalFeature("HBondDonor", "HBondDonor1", geom.Point3D(1.0, 2.0, 3.0))
      self.failUnless(ffeat.GetFamily() == "HBondDonor")
      self.failUnless(ffeat.GetType() == "HBondDonor1")
      pos = ffeat.GetPos()
      self.failUnless(ptFeq(pos, geom.Point3D(1.0, 2.0, 3.0)))

      ffeat = ChemicalFeatures.FreeChemicalFeature(type="HBondDonor1", family="HBondDonor",
                                       loc=geom.Point3D(1.0, 2.0, 3.0))
      self.failUnless(ffeat.GetFamily() == "HBondDonor")
      self.failUnless(ffeat.GetType() == "HBondDonor1")
      pos = ffeat.GetPos()
      self.failUnless(ptFeq(pos, geom.Point3D(1.0, 2.0, 3.0)))
      
    def testPickle(self):
      ffeat = ChemicalFeatures.FreeChemicalFeature("HBondDonor", "HBondDonor1", geom.Point3D(1.0, 2.0, 3.0))
      pkl = cPickle.dumps(ffeat)
      ffeat2 = cPickle.loads(pkl)
      self.failUnless(ffeat2.GetFamily()==ffeat.GetFamily());
      self.failUnless(ffeat2.GetType()==ffeat.GetType());
      self.failUnless(ptFeq(ffeat2.GetPos(),ffeat.GetPos()))

      #cPickle.dump(ffeat,file('feat.pkl','wb+'))
      inF = file(os.path.join(RDConfig.RDBaseDir,
                              'Code/ChemicalFeatures/Wrap/testData/feat.pkl'),'rb')
      ffeat2=cPickle.load(inF)
      self.failUnless(ffeat2.GetFamily()==ffeat.GetFamily());
      self.failUnless(ffeat2.GetType()==ffeat.GetType());
      self.failUnless(ptFeq(ffeat2.GetPos(),ffeat.GetPos()))

      
if __name__ == '__main__':
    print "Testing ChemicalFeatures Wrapper code:"
    unittest.main()
    
