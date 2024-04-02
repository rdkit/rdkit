# $Id$
#
# Copyright (C) 2006-2008 Greg Landrum
#
#  @@ All Rights Reserved @@
#

import io
import os
import pickle
import sys
import unittest

from rdkit import Chem, RDConfig
from rdkit.Chem import ChemicalFeatures
from rdkit.Geometry import rdGeometry as geom


def feq(v1, v2, tol2=1e-4):
  return abs(v1 - v2) <= tol2


def lstFeq(l1, l2, tol=1.e-4):
  if len(l1) != len(l2):
    return 0
  for ll1, ll2 in zip(l1, l2):
    if not feq(ll1, ll2, tol):
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
    ffeat.SetId(123)
    pos = ffeat.GetId()
    self.assertTrue(pos == 123)
    ffeat.SetFamily("HBondDonor")
    self.assertTrue(ffeat.GetFamily() == "HBondDonor")
    ffeat.SetPos(geom.Point3D(1.0, 2.0, 3.0))
    pos = ffeat.GetPos()
    self.assertTrue(ptFeq(pos, geom.Point3D(1.0, 2.0, 3.0)))
    ffeat.SetType("HBondDonor1")
    self.assertTrue(ffeat.GetType() == "HBondDonor1")

    ffeat = ChemicalFeatures.FreeChemicalFeature("HBondDonor", "HBondDonor1",
                                                 geom.Point3D(1.0, 2.0, 3.0))
    self.assertTrue(ffeat.GetId() == -1)
    self.assertTrue(ffeat.GetFamily() == "HBondDonor")
    self.assertTrue(ffeat.GetType() == "HBondDonor1")

    ffeat = ChemicalFeatures.FreeChemicalFeature("HBondDonor", "HBondDonor1",
                                                 geom.Point3D(1.0, 2.0, 3.0), id=123)
    self.assertTrue(ffeat.GetId() == 123)
    self.assertTrue(ffeat.GetFamily() == "HBondDonor")
    self.assertTrue(ffeat.GetType() == "HBondDonor1")

    pos = ffeat.GetPos()
    self.assertTrue(ptFeq(pos, geom.Point3D(1.0, 2.0, 3.0)))

    ffeat = ChemicalFeatures.FreeChemicalFeature(id=123, type="HBondDonor1", family="HBondDonor",
                                                 loc=geom.Point3D(1.0, 2.0, 3.0))
    self.assertTrue(ffeat.GetId() == 123)
    self.assertTrue(ffeat.GetFamily() == "HBondDonor")
    self.assertTrue(ffeat.GetType() == "HBondDonor1")
    pos = ffeat.GetPos()
    self.assertTrue(ptFeq(pos, geom.Point3D(1.0, 2.0, 3.0)))

  def testPickle(self):
    ffeat = ChemicalFeatures.FreeChemicalFeature("HBondDonor", "HBondDonor1",
                                                 geom.Point3D(1.0, 2.0, 3.0), 123)
    pkl = pickle.dumps(ffeat)
    ffeat2 = pickle.loads(pkl, encoding='bytes')
    self.assertTrue(ffeat2.GetId() == ffeat.GetId())
    self.assertTrue(ffeat2.GetFamily() == ffeat.GetFamily())
    self.assertTrue(ffeat2.GetType() == ffeat.GetType())
    self.assertTrue(ptFeq(ffeat2.GetPos(), ffeat.GetPos()))

    # Check that the old pickled versions have not been broken
    inTF = open(os.path.join(RDConfig.RDBaseDir, 'Code/ChemicalFeatures/Wrap/testData/feat.pkl'),
                'r')
    buf = inTF.read().replace('\r\n', '\n').encode('utf-8')
    inTF.close()
    inF = io.BytesIO(buf)
    ffeat2 = pickle.load(inF, encoding='bytes')
    # this version (1.0) does not have an id in the byte stream
    self.assertTrue(ffeat2.GetFamily() == ffeat.GetFamily())
    self.assertTrue(ffeat2.GetType() == ffeat.GetType())
    self.assertTrue(ptFeq(ffeat2.GetPos(), ffeat.GetPos()))

    # Test the new version also has the id and works as expected

    # uncomment the following to generate (overwrite) new version of pickled
    # data file
    #pickle.dump(ffeat,file(os.path.join(RDConfig.RDBaseDir, 'Code/ChemicalFeatures/Wrap/testData/featv2.pkl'),'wb+'))
    inTF = open(os.path.join(RDConfig.RDBaseDir, 'Code/ChemicalFeatures/Wrap/testData/featv2.pkl'),
                'r')
    buf = inTF.read().replace('\r\n', '\n').encode('utf-8')
    inTF.close()
    inF = io.BytesIO(buf)
    ffeat2 = pickle.load(inF, encoding='bytes')
    self.assertTrue(ffeat2.GetId() == ffeat.GetId())
    self.assertTrue(ffeat2.GetFamily() == ffeat.GetFamily())
    self.assertTrue(ffeat2.GetType() == ffeat.GetType())
    self.assertTrue(ptFeq(ffeat2.GetPos(), ffeat.GetPos()))


if __name__ == '__main__':
  print("Testing ChemicalFeatures Wrapper code:")
  unittest.main()
