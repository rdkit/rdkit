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
import unittest

from rdkit.Chem.FeatMaps import FeatMapParser, FeatMapUtils
from rdkit.Geometry import Point3D
from rdkit.Chem.FeatMaps.FeatMapUtils import MergeMetric


def feq(n1, n2, tol=1e-4):
  return abs(n1 - n2) <= tol


def pteq(p1, p2, tol=1e-4):
  return feq((p1 - p2).LengthSq(), 0.0, tol)


class TestCase(unittest.TestCase):

  def setUp(self):
    self.paramTxt = """
BeginParams
  family=Acceptor radius=0.5 profile=Box
EndParams
"""
    self.p = FeatMapParser.FeatMapParser()

  def test1Basics(self):
    txt = self.paramTxt + """
BeginPoints
  family=Acceptor pos=(1.0, 0.0, 0.0) weight=1.0
  family=Acceptor pos=(1.1, 0.0, 0.0) weight=1.0
  family=Acceptor pos=(3.0, 0.0, 0.0) weight=1.0
EndPoints
    """
    self.p.SetData(txt)
    fm1 = self.p.Parse()

    self.assertTrue(fm1.GetNumFeatures() == 3)
    self.assertFalse(FeatMapUtils.MergeFeatPoints(fm1))
    self.assertTrue(FeatMapUtils.MergeFeatPoints(fm1, FeatMapUtils.MergeMetric.Distance))
    self.assertTrue(fm1.GetNumFeatures() == 2)
    self.assertTrue(pteq(fm1.GetFeature(0).GetPos(), Point3D(1.05, 0, 0)))
    self.assertTrue(pteq(fm1.GetFeature(1).GetPos(), Point3D(3.0, 0, 0)))

    txt = self.paramTxt + """
BeginPoints
  family=Acceptor pos=(1.0, 0.0, 0.0) weight=1.0
  family=Acceptor pos=(1.1, 0.0, 0.0) weight=1.0
  family=Acceptor pos=(3.0, 0.0, 0.0) weight=1.0
  family=Acceptor pos=(4.0, 0.0, 0.0) weight=1.0
EndPoints
    """
    self.p.SetData(txt)
    fm1 = self.p.Parse()

    self.assertTrue(fm1.GetNumFeatures() == 4)
    self.assertTrue(FeatMapUtils.MergeFeatPoints(fm1, FeatMapUtils.MergeMetric.Distance))
    self.assertTrue(fm1.GetNumFeatures() == 2)
    self.assertTrue(pteq(fm1.GetFeature(0).GetPos(), Point3D(1.05, 0, 0)))
    self.assertTrue(pteq(fm1.GetFeature(1).GetPos(), Point3D(3.5, 0, 0)))

    txt = self.paramTxt + """
BeginPoints
  family=Acceptor pos=(1.0, 0.0, 0.0) weight=1.0
  family=Acceptor pos=(1.2, 0.0, 0.0) weight=1.0
  family=Acceptor pos=(1.3, 0.0, 0.0) weight=1.0
  family=Acceptor pos=(4.0, 0.0, 0.0) weight=1.0
EndPoints
    """
    self.p.SetData(txt)
    fm1 = self.p.Parse()

    self.assertTrue(fm1.GetNumFeatures() == 4)
    self.assertTrue(FeatMapUtils.MergeFeatPoints(fm1, FeatMapUtils.MergeMetric.Distance))
    self.assertTrue(fm1.GetNumFeatures() == 3)
    self.assertTrue(pteq(fm1.GetFeature(0).GetPos(), Point3D(1.00, 0, 0)))
    self.assertTrue(pteq(fm1.GetFeature(1).GetPos(), Point3D(1.25, 0, 0)))
    self.assertTrue(pteq(fm1.GetFeature(2).GetPos(), Point3D(4.0, 0, 0)))

    txt = self.paramTxt + """
BeginPoints
  family=Acceptor pos=(1.0, 0.0, 0.0) weight=1.0
  family=Acceptor pos=(1.2, 0.0, 0.0) weight=3.0
  family=Acceptor pos=(1.3, 0.0, 0.0) weight=1.0
  family=Acceptor pos=(4.0, 0.0, 0.0) weight=1.0
EndPoints
    """
    self.p.SetData(txt)
    fm1 = self.p.Parse()

    self.assertTrue(fm1.GetNumFeatures() == 4)
    self.assertTrue(
      FeatMapUtils.MergeFeatPoints(fm1, FeatMapUtils.MergeMetric.Distance,
                                   mergeMethod=FeatMapUtils.MergeMethod.Average))
    self.assertTrue(fm1.GetNumFeatures() == 3)
    self.assertTrue(pteq(fm1.GetFeature(0).GetPos(), Point3D(1.00, 0, 0)))
    self.assertTrue(pteq(fm1.GetFeature(1).GetPos(), Point3D(1.25, 0, 0)))
    self.assertTrue(pteq(fm1.GetFeature(2).GetPos(), Point3D(4.0, 0, 0)))

    self.p.SetData(txt)
    fm1 = self.p.Parse()
    self.assertTrue(fm1.GetNumFeatures() == 4)
    self.assertTrue(
      FeatMapUtils.MergeFeatPoints(fm1, FeatMapUtils.MergeMetric.Distance,
                                   mergeMethod=FeatMapUtils.MergeMethod.WeightedAverage))
    self.assertTrue(fm1.GetNumFeatures() == 3)
    self.assertTrue(pteq(fm1.GetFeature(0).GetPos(), Point3D(1.00, 0, 0)))
    self.assertTrue(pteq(fm1.GetFeature(1).GetPos(), Point3D(1.225, 0, 0)))
    self.assertTrue(pteq(fm1.GetFeature(2).GetPos(), Point3D(4.0, 0, 0)))

    self.p.SetData(txt)
    fm1 = self.p.Parse()
    self.assertTrue(fm1.GetNumFeatures() == 4)
    self.assertTrue(
      FeatMapUtils.MergeFeatPoints(fm1, FeatMapUtils.MergeMetric.Distance,
                                   mergeMethod=FeatMapUtils.MergeMethod.UseLarger))
    self.assertTrue(fm1.GetNumFeatures() == 3)
    self.assertTrue(pteq(fm1.GetFeature(0).GetPos(), Point3D(1.00, 0, 0)))
    self.assertTrue(pteq(fm1.GetFeature(1).GetPos(), Point3D(1.2, 0, 0)))
    self.assertTrue(pteq(fm1.GetFeature(2).GetPos(), Point3D(4.0, 0, 0)))

    # Order of points doesn't matter
    txt = self.paramTxt + """
BeginPoints
  family=Acceptor pos=(1.0, 0.0, 0.0) weight=1.0
  family=Acceptor pos=(1.3, 0.0, 0.0) weight=1.0
  family=Acceptor pos=(1.2, 0.0, 0.0) weight=3.0
  family=Acceptor pos=(4.0, 0.0, 0.0) weight=1.0
EndPoints
    """
    self.p.SetData(txt)
    fm1 = self.p.Parse()
    self.assertTrue(fm1.GetNumFeatures() == 4)
    self.assertTrue(
      FeatMapUtils.MergeFeatPoints(fm1, FeatMapUtils.MergeMetric.Distance,
                                   mergeMethod=FeatMapUtils.MergeMethod.UseLarger))
    self.assertTrue(fm1.GetNumFeatures() == 3)
    self.assertTrue(pteq(fm1.GetFeature(0).GetPos(), Point3D(1.00, 0, 0)))
    self.assertTrue(pteq(fm1.GetFeature(1).GetPos(), Point3D(1.2, 0, 0)))
    self.assertTrue(pteq(fm1.GetFeature(2).GetPos(), Point3D(4.0, 0, 0)))

    self.assertRaises(ValueError, FeatMapUtils.MergeFeatPoints, fm1,
                      FeatMapUtils.MergeMetric.Distance, mergeMethod='typo')

  def test1BasicsRepeated(self):
    txt = self.paramTxt + """
BeginPoints
  family=Acceptor pos=(0.7, 0.0, 0.0) weight=1.0
  family=Acceptor pos=(1.0, 0.0, 0.0) weight=1.0
  family=Acceptor pos=(1.2, 0.0, 0.0) weight=1.0
  family=Acceptor pos=(1.3, 0.0, 0.0) weight=1.0
  family=Acceptor pos=(4.0, 0.0, 0.0) weight=1.0
EndPoints
    """
    self.p.SetData(txt)
    fm1 = self.p.Parse()

    self.assertTrue(fm1.GetNumFeatures() == 5)
    self.assertTrue(FeatMapUtils.MergeFeatPoints(fm1, FeatMapUtils.MergeMetric.Distance))
    self.assertTrue(fm1.GetNumFeatures() == 4)
    self.assertTrue(pteq(fm1.GetFeature(0).GetPos(), Point3D(0.7, 0, 0)))
    self.assertTrue(pteq(fm1.GetFeature(1).GetPos(), Point3D(1.0, 0, 0)))
    self.assertTrue(pteq(fm1.GetFeature(2).GetPos(), Point3D(1.25, 0, 0)))
    self.assertTrue(pteq(fm1.GetFeature(3).GetPos(), Point3D(4.0, 0, 0)))

    self.assertTrue(FeatMapUtils.MergeFeatPoints(fm1, FeatMapUtils.MergeMetric.Distance))
    self.assertTrue(fm1.GetNumFeatures() == 3)
    self.assertTrue(pteq(fm1.GetFeature(0).GetPos(), Point3D(0.7, 0, 0)))
    self.assertTrue(pteq(fm1.GetFeature(1).GetPos(), Point3D(1.125, 0, 0)))
    self.assertTrue(pteq(fm1.GetFeature(2).GetPos(), Point3D(4.0, 0, 0)))

    self.assertTrue(FeatMapUtils.MergeFeatPoints(fm1, FeatMapUtils.MergeMetric.Distance))
    self.assertTrue(fm1.GetNumFeatures() == 2)
    self.assertTrue(pteq(fm1.GetFeature(0).GetPos(), Point3D(0.9125, 0, 0)))
    self.assertTrue(pteq(fm1.GetFeature(1).GetPos(), Point3D(4.0, 0, 0)))

  def test2ScoreBasics(self):
    txt = self.paramTxt + """
BeginPoints
  family=Acceptor pos=(1.0, 0.0, 0.0) weight=1.0
  family=Acceptor pos=(1.2, 0.0, 0.0) weight=3.0
  family=Acceptor pos=(4.0, 0.0, 0.0) weight=1.0
EndPoints
    """
    self.p.SetData(txt)
    fm1 = self.p.Parse()

    self.assertTrue(fm1.GetNumFeatures() == 3)
    self.assertTrue(
      FeatMapUtils.MergeFeatPoints(fm1, FeatMapUtils.MergeMetric.Overlap,
                                   mergeMethod=FeatMapUtils.MergeMethod.Average))
    self.assertTrue(fm1.GetNumFeatures() == 2)
    self.assertTrue(pteq(fm1.GetFeature(0).GetPos(), Point3D(1.1, 0, 0)))
    self.assertTrue(pteq(fm1.GetFeature(1).GetPos(), Point3D(4.0, 0, 0)))

    txt = self.paramTxt + """
BeginPoints
  family=Acceptor pos=(1.0, 0.0, 0.0) weight=1.0
  family=Acceptor pos=(1.1, 0.0, 0.0) weight=1.0
  family=Acceptor pos=(1.3, 0.0, 0.0) weight=3.0
  family=Acceptor pos=(4.0, 0.0, 0.0) weight=1.0
EndPoints
    """
    self.p.SetData(txt)
    fm1 = self.p.Parse()

    self.assertTrue(fm1.GetNumFeatures() == 4)
    self.assertTrue(
      FeatMapUtils.MergeFeatPoints(fm1, FeatMapUtils.MergeMetric.Overlap,
                                   mergeMethod=FeatMapUtils.MergeMethod.Average))
    self.assertTrue(fm1.GetNumFeatures() == 3)
    self.assertTrue(pteq(fm1.GetFeature(0).GetPos(), Point3D(1.15, 0, 0)))
    self.assertTrue(pteq(fm1.GetFeature(1).GetPos(), Point3D(1.1, 0, 0)))
    self.assertTrue(pteq(fm1.GetFeature(2).GetPos(), Point3D(4.0, 0, 0)))

    txt = self.paramTxt + """
BeginPoints
  family=Acceptor pos=(1.0, 0.0, 0.0) weight=1.0
  family=Acceptor pos=(1.2, 0.0, 0.0) weight=1.0
  family=Acceptor pos=(1.6, 0.0, 0.0) weight=3.0
  family=Acceptor pos=(4.0, 0.0, 0.0) weight=1.0
EndPoints
    """
    self.p.SetData(txt)
    fm1 = self.p.Parse()

    self.assertTrue(fm1.GetNumFeatures() == 4)
    self.assertTrue(
      FeatMapUtils.MergeFeatPoints(fm1, FeatMapUtils.MergeMetric.Overlap,
                                   mergeMethod=FeatMapUtils.MergeMethod.Average))
    self.assertTrue(fm1.GetNumFeatures() == 3)
    self.assertTrue(pteq(fm1.GetFeature(0).GetPos(), Point3D(1.0, 0, 0)))
    self.assertTrue(pteq(fm1.GetFeature(1).GetPos(), Point3D(1.4, 0, 0)))
    self.assertTrue(pteq(fm1.GetFeature(2).GetPos(), Point3D(4.0, 0, 0)))

  def test_CombineFeatMaps(self):
    data1 = self.paramTxt + "BeginPoints\nfamily=Acceptor pos=(1.0, 0.0, 0.0) weight=1.0\nEndPoints"
    data2 = self.paramTxt + "BeginPoints\nfamily=Acceptor pos=(0.0, 0.0, 0.0) weight=1.0\nEndPoints"

    self.p.SetData(data1)
    fm1 = self.p.Parse()
    self.assertEqual(fm1.GetNumFeatures(), 1)
    self.p.SetData(data2)
    fm2 = self.p.Parse()
    self.assertEqual(fm2.GetNumFeatures(), 1)

    # Combine without merging
    fm12 = FeatMapUtils.CombineFeatMaps(fm1, fm2)
    self.assertEqual(fm12.GetNumFeatures(), 2)

    # Combine with distance based merging
    fm12 = FeatMapUtils.CombineFeatMaps(fm1, fm2, mergeMetric=MergeMetric.Distance, mergeTol=1.0)
    self.assertEqual(fm12.GetNumFeatures(), 2)
    fm12 = FeatMapUtils.CombineFeatMaps(fm1, fm2, mergeMetric=MergeMetric.Distance, mergeTol=1.1)
    self.assertEqual(fm12.GetNumFeatures(), 1)

    self.assertRaises(ValueError, FeatMapUtils.CombineFeatMaps, fm1, fm2, mergeMetric='typo')


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
