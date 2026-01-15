import copy
import math
import os
import pickle
import sys
import unittest

from rdkit import RDConfig
from rdkit.Geometry import rdGeometry as geom


def feq(v1, v2, tol=1.0e-4):
  return abs(v1 - v2) < tol


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def test1aPoint3D(self):
    pt = geom.Point3D()
    self.assertTrue(feq(pt.x, 0.0))
    self.assertTrue(feq(pt.y, 0.0))
    self.assertTrue(feq(pt.z, 0.0))

    pt = geom.Point3D(3., 4., 5.)
    self.assertTrue(feq(pt.x, 3.0))
    self.assertTrue(feq(pt.y, 4.0))
    self.assertTrue(feq(pt.z, 5.0))
    self.assertTrue(feq(pt[0], 3.0))
    self.assertTrue(feq(pt[1], 4.0))
    self.assertTrue(feq(pt[2], 5.0))
    self.assertTrue(feq(pt[-3], 3.0))
    self.assertTrue(feq(pt[-2], 4.0))
    self.assertTrue(feq(pt[-1], 5.0))
    lst = list(pt)
    self.assertTrue(feq(lst[0], 3.0))
    self.assertTrue(feq(lst[1], 4.0))
    self.assertTrue(feq(lst[2], 5.0))

    pt2 = geom.Point3D(1., 1., 1.)

    pt3 = pt + pt2
    self.assertTrue(feq(pt3.x, 4.0))
    self.assertTrue(feq(pt3.y, 5.0))
    self.assertTrue(feq(pt3.z, 6.0))

    pt += pt2
    self.assertTrue(feq(pt.x, 4.0))
    self.assertTrue(feq(pt.y, 5.0))
    self.assertTrue(feq(pt.z, 6.0))

    pt3 = pt - pt2
    self.assertTrue(feq(pt3.x, 3.0))
    self.assertTrue(feq(pt3.y, 4.0))
    self.assertTrue(feq(pt3.z, 5.0))

    pt -= pt2
    self.assertTrue(feq(pt.x, 3.0))
    self.assertTrue(feq(pt.y, 4.0))
    self.assertTrue(feq(pt.z, 5.0))

    pt *= 2.0
    self.assertTrue(feq(pt.x, 6.0))
    self.assertTrue(feq(pt.y, 8.0))
    self.assertTrue(feq(pt.z, 10.0))

    pt /= 2
    self.assertTrue(feq(pt.x, 3.0))
    self.assertTrue(feq(pt.y, 4.0))
    self.assertTrue(feq(pt.z, 5.0))
    self.assertTrue(feq(pt.Length(), 7.0711))
    self.assertTrue(feq(pt.LengthSq(), 50.0))
    pt.Normalize()
    self.assertTrue(feq(pt.Length(), 1.0))

    pt1 = geom.Point3D(1.0, 0.0, 0.0)
    pt2 = geom.Point3D(2.0 * math.cos(math.pi / 6), 2.0 * math.sin(math.pi / 6), 0.0)
    ang = pt1.AngleTo(pt2)
    self.assertTrue(feq(ang, math.pi / 6))

    prod = pt1.DotProduct(pt2)
    self.assertTrue(feq(prod, 2.0 * math.cos(math.pi / 6)))

    pt3 = pt1.CrossProduct(pt2)
    self.assertTrue(feq(pt3.x, 0.0))
    self.assertTrue(feq(pt3.y, 0.0))
    self.assertTrue(feq(pt3.z, 1.0))

  def test1bPoint2D(self):
    pt = geom.Point2D()
    self.assertTrue(feq(pt.x, 0.0))
    self.assertTrue(feq(pt.y, 0.0))

    pt = geom.Point2D(3., 4.)
    self.assertTrue(feq(pt.x, 3.0))
    self.assertTrue(feq(pt.y, 4.0))

    self.assertTrue(feq(pt.x, 3.0))
    self.assertTrue(feq(pt.y, 4.0))
    self.assertTrue(feq(pt[0], 3.0))
    self.assertTrue(feq(pt[1], 4.0))
    self.assertTrue(feq(pt[-2], 3.0))
    self.assertTrue(feq(pt[-1], 4.0))
    lst = list(pt)
    self.assertTrue(feq(lst[0], 3.0))
    self.assertTrue(feq(lst[1], 4.0))

    pt2 = geom.Point2D(1., 1.)

    pt3 = pt + pt2
    self.assertTrue(feq(pt3.x, 4.0))
    self.assertTrue(feq(pt3.y, 5.0))

    pt += pt2
    self.assertTrue(feq(pt.x, 4.0))
    self.assertTrue(feq(pt.y, 5.0))

    pt3 = pt - pt2
    self.assertTrue(feq(pt3.x, 3.0))
    self.assertTrue(feq(pt3.y, 4.0))

    pt -= pt2
    self.assertTrue(feq(pt.x, 3.0))
    self.assertTrue(feq(pt.y, 4.0))

    pt *= 2.0
    self.assertTrue(feq(pt.x, 6.0))
    self.assertTrue(feq(pt.y, 8.0))

    pt /= 2
    self.assertTrue(feq(pt.x, 3.0))
    self.assertTrue(feq(pt.y, 4.0))
    self.assertTrue(feq(pt.Length(), 5.0))
    self.assertTrue(feq(pt.LengthSq(), 25.0))
    pt.Normalize()
    self.assertTrue(feq(pt.Length(), 1.0))

    pt1 = geom.Point2D(1.0, 0.0)
    pt2 = geom.Point2D(2.0 * math.cos(math.pi / 6), 2.0 * math.sin(math.pi / 6))
    ang = pt1.AngleTo(pt2)
    self.assertTrue(feq(ang, math.pi / 6))

    prod = pt1.DotProduct(pt2)
    self.assertTrue(feq(prod, 2.0 * math.cos(math.pi / 6)))

  def test1cPointND(self):
    dim = 4
    pt = geom.PointND(4)
    for i in range(dim):
      self.assertTrue(feq(pt[i], 0.0))

    pt[0] = 3
    pt[3] = 4
    self.assertTrue(feq(pt[0], 3.0))
    self.assertTrue(feq(pt[3], 4.0))
    self.assertTrue(feq(pt[-4], 3.0))
    self.assertTrue(feq(pt[-1], 4.0))
    lst = list(pt)
    self.assertTrue(feq(lst[0], 3.0))
    self.assertTrue(feq(lst[3], 4.0))

    pt2 = geom.PointND(4)
    pt2[0] = 1.
    pt2[2] = 1.

    pt3 = pt + pt2
    self.assertTrue(feq(pt3[0], 4.0))
    self.assertTrue(feq(pt3[2], 1.0))
    self.assertTrue(feq(pt3[3], 4.0))

    pt += pt2
    self.assertTrue(feq(pt[0], 4.0))
    self.assertTrue(feq(pt[2], 1.0))
    self.assertTrue(feq(pt[3], 4.0))

    pt3 = pt - pt2
    self.assertTrue(feq(pt3[0], 3.0))
    self.assertTrue(feq(pt3[2], 0.0))
    self.assertTrue(feq(pt3[3], 4.0))

    pt -= pt2
    self.assertTrue(feq(pt[0], 3.0))
    self.assertTrue(feq(pt[2], 0.0))
    self.assertTrue(feq(pt[3], 4.0))

    pt *= 2.0
    self.assertTrue(feq(pt[0], 6.0))
    self.assertTrue(feq(pt[1], 0.0))
    self.assertTrue(feq(pt[2], 0.0))
    self.assertTrue(feq(pt[3], 8.0))

    pt /= 2
    self.assertTrue(feq(pt[0], 3.0))
    self.assertTrue(feq(pt[1], 0.0))
    self.assertTrue(feq(pt[2], 0.0))
    self.assertTrue(feq(pt[3], 4.0))

    self.assertTrue(feq(pt.Length(), 5.0))
    self.assertTrue(feq(pt.LengthSq(), 25.0))
    pt.Normalize()
    self.assertTrue(feq(pt.Length(), 1.0))

    pkl = pickle.dumps(pt)
    pt2 = pickle.loads(pkl)
    self.assertTrue(len(pt) == len(pt2))
    for i in range(len(pt)):
      self.assertTrue(feq(pt2[i], pt[i]))

  def testPointPickles(self):
    pt = geom.Point3D(2.0, -3.0, 1.0)
    pt2 = pickle.loads(pickle.dumps(pt))
    self.assertTrue(feq(pt.x, pt2.x, 1e-6))
    self.assertTrue(feq(pt.y, pt2.y, 1e-6))
    self.assertTrue(feq(pt.z, pt2.z, 1e-6))

    pt = geom.Point2D(2.0, -4.0)
    pt2 = pickle.loads(pickle.dumps(pt))
    self.assertTrue(feq(pt.x, pt2.x, 1e-6))
    self.assertTrue(feq(pt.y, pt2.y, 1e-6))

    pt = geom.PointND(5)
    for i in range(5):
      pt[i] = float(i) * 1.5 - 2.0
    pt2 = pickle.loads(pickle.dumps(pt))
    self.assertTrue(len(pt) == len(pt2))
    for i in range(len(pt)):
      self.assertTrue(feq(pt2[i], pt[i], 1e-6))

  def test8InitPoint2DFromPoint3D(self):
    p3 = geom.Point3D(1., 2., 3.)
    p2 = geom.Point2D(p3)
    self.assertEqual(p2.x, p3.x)
    self.assertEqual(p2.y, p3.y)

  def test1bPoint2D(self):
    pt = geom.Point2D()
    self.assertTrue(feq(pt.x, 0.0))
    self.assertTrue(feq(pt.y, 0.0))

    pt = geom.Point2D(3., 4.)
    self.assertTrue(feq(pt.x, 3.0))
    self.assertTrue(feq(pt.y, 4.0))

    self.assertTrue(feq(pt.x, 3.0))
    self.assertTrue(feq(pt.y, 4.0))
    self.assertTrue(feq(pt[0], 3.0))
    self.assertTrue(feq(pt[1], 4.0))
    self.assertTrue(feq(pt[-2], 3.0))
    self.assertTrue(feq(pt[-1], 4.0))
    lst = list(pt)
    self.assertTrue(feq(lst[0], 3.0))
    self.assertTrue(feq(lst[1], 4.0))

    pt2 = geom.Point2D(1., 1.)

    pt3 = pt + pt2
    self.assertTrue(feq(pt3.x, 4.0))
    self.assertTrue(feq(pt3.y, 5.0))

    pt += pt2
    self.assertTrue(feq(pt.x, 4.0))
    self.assertTrue(feq(pt.y, 5.0))

    pt3 = pt - pt2
    self.assertTrue(feq(pt3.x, 3.0))
    self.assertTrue(feq(pt3.y, 4.0))

    pt -= pt2
    self.assertTrue(feq(pt.x, 3.0))
    self.assertTrue(feq(pt.y, 4.0))

    pt *= 2.0
    self.assertTrue(feq(pt.x, 6.0))
    self.assertTrue(feq(pt.y, 8.0))

    pt /= 2
    self.assertTrue(feq(pt.x, 3.0))
    self.assertTrue(feq(pt.y, 4.0))
    self.assertTrue(feq(pt.Length(), 5.0))
    self.assertTrue(feq(pt.LengthSq(), 25.0))
    pt.Normalize()
    self.assertTrue(feq(pt.Length(), 1.0))

    pt1 = geom.Point2D(1.0, 0.0)
    pt2 = geom.Point2D(2.0 * math.cos(math.pi / 6), 2.0 * math.sin(math.pi / 6))
    ang = pt1.AngleTo(pt2)
    self.assertTrue(feq(ang, math.pi / 6))

    prod = pt1.DotProduct(pt2)
    self.assertTrue(feq(prod, 2.0 * math.cos(math.pi / 6)))

  def test6Dihedrals(self):
    p1 = geom.Point3D(1, 0, 0)
    p2 = geom.Point3D(0, 0, 0)
    p3 = geom.Point3D(0, 1, 0)

    p4 = geom.Point3D(.5, 1, .5)
    ang = geom.ComputeDihedralAngle(p1, p2, p3, p4)
    self.assertAlmostEqual(ang, math.pi / 4, 4)
    ang = geom.ComputeSignedDihedralAngle(p1, p2, p3, p4)
    self.assertAlmostEqual(ang, -math.pi / 4, 4)

    p4 = geom.Point3D(-.5, 1, .5)
    ang = geom.ComputeDihedralAngle(p1, p2, p3, p4)
    self.assertAlmostEqual(ang, 3 * math.pi / 4, 4)
    ang = geom.ComputeSignedDihedralAngle(p1, p2, p3, p4)
    self.assertAlmostEqual(ang, -3 * math.pi / 4, 4)

    p4 = geom.Point3D(.5, 1, -.5)
    ang = geom.ComputeDihedralAngle(p1, p2, p3, p4)
    self.assertAlmostEqual(ang, math.pi / 4, 4)
    ang = geom.ComputeSignedDihedralAngle(p1, p2, p3, p4)
    self.assertAlmostEqual(ang, math.pi / 4, 4)

    p4 = geom.Point3D(-.5, 1, -.5)
    ang = geom.ComputeDihedralAngle(p1, p2, p3, p4)
    self.assertAlmostEqual(ang, 3 * math.pi / 4, 4)
    ang = geom.ComputeSignedDihedralAngle(p1, p2, p3, p4)
    self.assertAlmostEqual(ang, 3 * math.pi / 4, 4)

    p4 = geom.Point3D(0, 1, 1)
    ang = geom.ComputeDihedralAngle(p1, p2, p3, p4)
    self.assertAlmostEqual(ang, math.pi / 2, 4)
    ang = geom.ComputeSignedDihedralAngle(p1, p2, p3, p4)
    self.assertAlmostEqual(ang, -math.pi / 2, 4)

    p4 = geom.Point3D(0, 1, -1)
    ang = geom.ComputeDihedralAngle(p1, p2, p3, p4)
    self.assertAlmostEqual(ang, math.pi / 2, 4)
    ang = geom.ComputeSignedDihedralAngle(p1, p2, p3, p4)
    self.assertAlmostEqual(ang, math.pi / 2, 4)

    p4 = geom.Point3D(1, 1, 0)
    ang = geom.ComputeDihedralAngle(p1, p2, p3, p4)
    self.assertAlmostEqual(ang, 0, 4)
    ang = geom.ComputeSignedDihedralAngle(p1, p2, p3, p4)
    self.assertAlmostEqual(ang, 0, 4)

    p4 = geom.Point3D(-1, 1, 0)
    ang = geom.ComputeDihedralAngle(p1, p2, p3, p4)
    self.assertAlmostEqual(ang, math.pi, 4)
    ang = geom.ComputeSignedDihedralAngle(p1, p2, p3, p4)
    self.assertAlmostEqual(ang, math.pi, 4)


if __name__ == '__main__':
  print("Testing Geometry wrapper")
  unittest.main()
