from __future__ import print_function
import os,sys
import unittest
import copy
import math

from rdkit.six.moves import cPickle

from rdkit import RDConfig
from rdkit import DataStructs
from rdkit.Geometry import rdGeometry as geom

def feq(v1, v2, tol=1.0e-4):
    return abs(v1-v2) < tol

class TestCase(unittest.TestCase):
    def setUp(self) :
        pass

    def test1aPoint3D(self):
        pt = geom.Point3D();
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

        pt3 = pt+pt2
        self.assertTrue(feq(pt3.x, 4.0))
        self.assertTrue(feq(pt3.y, 5.0))
        self.assertTrue(feq(pt3.z, 6.0))
        
        
        pt += pt2
        self.assertTrue(feq(pt.x, 4.0))
        self.assertTrue(feq(pt.y, 5.0))
        self.assertTrue(feq(pt.z, 6.0))

        pt3 = pt-pt2
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
        pt2 = geom.Point3D(2.0*math.cos(math.pi/6), 2.0*math.sin(math.pi/6), 0.0)
        ang = pt1.AngleTo(pt2)
        self.assertTrue(feq(ang, math.pi/6))

        prod = pt1.DotProduct(pt2)
        self.assertTrue(feq(prod, 2.0*math.cos(math.pi/6)))

        pt3 = pt1.CrossProduct(pt2)
        self.assertTrue(feq(pt3.x, 0.0))
        self.assertTrue(feq(pt3.y, 0.0))
        self.assertTrue(feq(pt3.z, 1.0))



    def test1bPoint2D(self):
        pt = geom.Point2D();
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

        pt3 = pt+pt2
        self.assertTrue(feq(pt3.x, 4.0))
        self.assertTrue(feq(pt3.y, 5.0))
        
        pt += pt2
        self.assertTrue(feq(pt.x, 4.0))
        self.assertTrue(feq(pt.y, 5.0))
        
        pt3 = pt-pt2
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
        pt2 = geom.Point2D(2.0*math.cos(math.pi/6), 2.0*math.sin(math.pi/6))
        ang = pt1.AngleTo(pt2)
        self.assertTrue(feq(ang, math.pi/6))

        prod = pt1.DotProduct(pt2)
        self.assertTrue(feq(prod, 2.0*math.cos(math.pi/6)))

    def test1cPointND(self):
        dim=4
        pt = geom.PointND(4);
        for i in range(dim):
            self.assertTrue(feq(pt[i], 0.0))
        
        pt[0]=3
        pt[3]=4
        self.assertTrue(feq(pt[0], 3.0))
        self.assertTrue(feq(pt[3], 4.0))
        self.assertTrue(feq(pt[-4], 3.0))
        self.assertTrue(feq(pt[-1], 4.0))
        lst = list(pt)
        self.assertTrue(feq(lst[0], 3.0))
        self.assertTrue(feq(lst[3], 4.0))


        pt2 = geom.PointND(4)
        pt2[0]=1.
        pt2[2]=1.

        pt3 = pt+pt2
        self.assertTrue(feq(pt3[0], 4.0))
        self.assertTrue(feq(pt3[2], 1.0))
        self.assertTrue(feq(pt3[3], 4.0))
        
        pt += pt2
        self.assertTrue(feq(pt[0], 4.0))
        self.assertTrue(feq(pt[2], 1.0))
        self.assertTrue(feq(pt[3], 4.0))

        pt3 = pt-pt2
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

        pkl = cPickle.dumps(pt)
        pt2 = cPickle.loads(pkl)
        self.assertTrue(len(pt)==len(pt2))
        for i in range(len(pt)):
            self.assertTrue(feq(pt2[i],pt[i]))

    def test3UniformGrid(self):
        ugrid = geom.UniformGrid3D(20, 18, 15)
        self.assertTrue(ugrid.GetNumX() == 40)
        self.assertTrue(ugrid.GetNumY() == 36)
        self.assertTrue(ugrid.GetNumZ() == 30)
        dvect = ugrid.GetOccupancyVect()
        ugrid = geom.UniformGrid3D(20, 18, 15, 0.5, DataStructs.DiscreteValueType.TWOBITVALUE)
        dvect = ugrid.GetOccupancyVect()
        self.assertTrue(dvect.GetValueType() == DataStructs.DiscreteValueType.TWOBITVALUE)

        grd = geom.UniformGrid3D(10.0, 10.0, 10.0, 0.5)
        grd.SetSphereOccupancy(geom.Point3D(-2.0, -2.0, 0.0), 1.5, 0.25)
        grd.SetSphereOccupancy(geom.Point3D(-2.0, 2.0, 0.0), 1.5, 0.25)
        grd.SetSphereOccupancy(geom.Point3D(2.0, -2.0, 0.0), 1.5, 0.25)
        grd.SetSphereOccupancy(geom.Point3D(2.0, 2.0, 0.0), 1.5, 0.25)

        geom.WriteGridToFile(grd, "junk.grd")
        grd2 = geom.UniformGrid3D(10.0, 10.0, 10.0, 0.5)
        grd2.SetSphereOccupancy(geom.Point3D(-2.0, -2.0, 0.0), 1.5, 0.25)
        grd2.SetSphereOccupancy(geom.Point3D(-2.0, 2.0, 0.0), 1.5, 0.25)
        grd2.SetSphereOccupancy(geom.Point3D(2.0, -2.0, 0.0), 1.5, 0.25)

        dist = geom.TanimotoDistance(grd, grd2)
        self.assertTrue(dist == 0.25)
        dist = geom.ProtrudeDistance(grd, grd2)
        self.assertTrue(dist == 0.25)
        dist = geom.ProtrudeDistance(grd2, grd)
        self.assertTrue(dist==0.0)

        grd2 = geom.UniformGrid3D(10.0, 10.0, 10.0, 0.5, DataStructs.DiscreteValueType.FOURBITVALUE)
        grd2.SetSphereOccupancy(geom.Point3D(-2.0, -2.0, 0.0), 1.5, 0.25, 3)
        grd2.SetSphereOccupancy(geom.Point3D(-2.0, 2.0, 0.0), 1.5, 0.25, 3)
        grd2.SetSphereOccupancy(geom.Point3D(2.0, -2.0, 0.0), 1.5, 0.25, 3)
        self.assertRaises(ValueError, lambda : geom.TanimotoDistance(grd, grd2))

        grd2 = geom.UniformGrid3D(10.0, 10.0, 10.0, 1.0)
        self.assertRaises(ValueError, lambda : geom.TanimotoDistance(grd, grd2))

        grd2 = geom.UniformGrid3D(11.0, 10.0, 10.0, 1.0)
        self.assertRaises(ValueError, lambda : geom.TanimotoDistance(grd, grd2))

    def testSymmetry(self):
        grd = geom.UniformGrid3D(10.0, 10.0, 10.0, 0.5)
        grd.SetSphereOccupancy(geom.Point3D(-2.2, -2.0, 0.0), 1.65, 0.25)
        grd.SetSphereOccupancy(geom.Point3D(2.2, -2.0, 0.0), 1.65, 0.25)
        
        bPt1 = geom.Point3D(-4.0, -2.0, -2.0)
        bPt2 = geom.Point3D(4.0, -2.0, -2.0)
        for k in range(8) :
            bPt1 += geom.Point3D(0.0, 0.0, 0.5)
            bPt2 += geom.Point3D(0.0, 0.0, 0.5)
            for j in range(8):
                 bPt1 += geom.Point3D(0.0, 0.5, 0.0)
                 bPt2 += geom.Point3D(0.0, 0.5, 0.0)
                 for i in range(8) :
                     bPt1 += geom.Point3D(0.5, 0.0, 0.0)
                     bPt2 -= geom.Point3D(0.5, 0.0, 0.0)
                     self.assertTrue(grd.GetValPoint(bPt1) == grd.GetValPoint(bPt2))
                 
                 bPt1.x = -4.0
                 bPt2.x = 4.0
            bPt1.y = -2.0
            bPt2.y = -2.0

    def testPointPickles(self):
        pt = geom.Point3D(2.0,-3.0,1.0)
        pt2 = cPickle.loads(cPickle.dumps(pt))
        self.assertTrue(feq(pt.x,pt2.x,1e-6))
        self.assertTrue(feq(pt.y,pt2.y,1e-6))
        self.assertTrue(feq(pt.z,pt2.z,1e-6))

        pt = geom.Point2D(2.0,-4.0)
        pt2 = cPickle.loads(cPickle.dumps(pt))
        self.assertTrue(feq(pt.x,pt2.x,1e-6))
        self.assertTrue(feq(pt.y,pt2.y,1e-6))

    def test4GridPickles(self):
        grd = geom.UniformGrid3D(10.0, 9.0, 8.0, 0.5)
        self.assertTrue(grd.GetNumX() == 20)
        self.assertTrue(grd.GetNumY() == 18)
        self.assertTrue(grd.GetNumZ() == 16)
        grd.SetSphereOccupancy(geom.Point3D(-2.0, -2.0, 0.0), 1.5, 0.25)
        grd.SetSphereOccupancy(geom.Point3D(-2.0, 2.0, 0.0), 1.5, 0.25)
        grd.SetSphereOccupancy(geom.Point3D(2.0, -2.0, 0.0), 1.5, 0.25)
        grd.SetSphereOccupancy(geom.Point3D(2.0, 2.0, 0.0), 1.5, 0.25)

        self.assertTrue(geom.TanimotoDistance(grd,grd)==0.0)

        grd2 = cPickle.loads(cPickle.dumps(grd))
        self.assertTrue(grd2.GetNumX() == 20)
        self.assertTrue(grd2.GetNumY() == 18)
        self.assertTrue(grd2.GetNumZ() == 16)
        self.assertTrue(geom.TanimotoDistance(grd,grd2)==0.0)
        
    def test5GridOps(self):
        grd = geom.UniformGrid3D(10, 10, 10)
        grd.SetSphereOccupancy(geom.Point3D(-2.0, -2.0, 0.0), 1.0, 0.25)
        grd.SetSphereOccupancy(geom.Point3D(-2.0, 2.0, 0.0), 1.0, 0.25)

        grd2 = geom.UniformGrid3D(10, 10, 10)
        grd2.SetSphereOccupancy(geom.Point3D(2.0, -2.0, 0.0), 1.0, 0.25)
        grd2.SetSphereOccupancy(geom.Point3D(2.0, 2.0, 0.0), 1.0, 0.25)

        self.assertTrue(geom.TanimotoDistance(grd,grd)==0.0)
        self.assertTrue(geom.TanimotoDistance(grd,grd2)==1.0)

        grd3 = copy.deepcopy(grd)
        grd3 |= grd2
        self.assertTrue(geom.TanimotoDistance(grd3,grd)==.5)
        self.assertTrue(geom.TanimotoDistance(grd3,grd2)==.5)

        grd3 = copy.deepcopy(grd)
        grd3 += grd2
        self.assertTrue(geom.TanimotoDistance(grd3,grd)==.5)
        self.assertTrue(geom.TanimotoDistance(grd3,grd2)==.5)

        grd3 -= grd
        self.assertTrue(geom.TanimotoDistance(grd3,grd)==1.0)
        self.assertTrue(geom.TanimotoDistance(grd3,grd2)==0)

        grd4 = geom.UniformGrid3D(10, 10, 10)
        grd4.SetSphereOccupancy(geom.Point3D(-2.0, -2.0, 0.0), 1.0, 0.25)
        grd4.SetSphereOccupancy(geom.Point3D(-2.0, 2.0, 0.0), 1.0, 0.25)
        grd4.SetSphereOccupancy(geom.Point3D(2.0, -2.0, 0.0), 1.0, 0.25)
        self.assertTrue(feq(geom.TanimotoDistance(grd4,grd),.3333))
        self.assertTrue(feq(geom.TanimotoDistance(grd4,grd2),.75))

        grd4&=grd2
        self.assertTrue(feq(geom.TanimotoDistance(grd4,grd),1.0))
        self.assertTrue(feq(geom.TanimotoDistance(grd4,grd2),.5))
        

    def test6Dihedrals(self):
        p1 = geom.Point3D(1,0,0)
        p2 = geom.Point3D(0,0,0)
        p3 = geom.Point3D(0,1,0)

        p4 = geom.Point3D(.5,1,.5)
        ang = geom.ComputeDihedralAngle(p1,p2,p3,p4)
        self.assertAlmostEqual(ang,math.pi/4,4)
        ang = geom.ComputeSignedDihedralAngle(p1,p2,p3,p4)
        self.assertAlmostEqual(ang,-math.pi/4,4)

        p4 = geom.Point3D(-.5,1,.5)
        ang = geom.ComputeDihedralAngle(p1,p2,p3,p4)
        self.assertAlmostEqual(ang,3*math.pi/4,4)
        ang = geom.ComputeSignedDihedralAngle(p1,p2,p3,p4)
        self.assertAlmostEqual(ang,-3*math.pi/4,4)

        p4 = geom.Point3D(.5,1,-.5)
        ang = geom.ComputeDihedralAngle(p1,p2,p3,p4)
        self.assertAlmostEqual(ang,math.pi/4,4)
        ang = geom.ComputeSignedDihedralAngle(p1,p2,p3,p4)
        self.assertAlmostEqual(ang,math.pi/4,4)
        
        p4 = geom.Point3D(-.5,1,-.5)
        ang = geom.ComputeDihedralAngle(p1,p2,p3,p4)
        self.assertAlmostEqual(ang,3*math.pi/4,4)
        ang = geom.ComputeSignedDihedralAngle(p1,p2,p3,p4)
        self.assertAlmostEqual(ang,3*math.pi/4,4)

        p4 = geom.Point3D(0,1,1)
        ang = geom.ComputeDihedralAngle(p1,p2,p3,p4)
        self.assertAlmostEqual(ang,math.pi/2,4)
        ang = geom.ComputeSignedDihedralAngle(p1,p2,p3,p4)
        self.assertAlmostEqual(ang,-math.pi/2,4)

        p4 = geom.Point3D(0,1,-1)
        ang = geom.ComputeDihedralAngle(p1,p2,p3,p4)
        self.assertAlmostEqual(ang,math.pi/2,4)
        ang = geom.ComputeSignedDihedralAngle(p1,p2,p3,p4)
        self.assertAlmostEqual(ang,math.pi/2,4)

        p4 = geom.Point3D(1,1,0)
        ang = geom.ComputeDihedralAngle(p1,p2,p3,p4)
        self.assertAlmostEqual(ang,0,4)
        ang = geom.ComputeSignedDihedralAngle(p1,p2,p3,p4)
        self.assertAlmostEqual(ang,0,4)

        p4 = geom.Point3D(-1,1,0)
        ang = geom.ComputeDihedralAngle(p1,p2,p3,p4)
        self.assertAlmostEqual(ang,math.pi,4)
        ang = geom.ComputeSignedDihedralAngle(p1,p2,p3,p4)
        self.assertAlmostEqual(ang,math.pi,4)

    def test7UniformGridIndices(self):
        ugrid = geom.UniformGrid3D(20, 18, 15)
        idx = ugrid.GetGridIndex(3,2,1)
        xi,yi,zi=ugrid.GetGridIndices(idx)
        self.assertEqual(xi,3)
        self.assertEqual(yi,2)
        self.assertEqual(zi,1)
        
if __name__=='__main__':
    print("Testing Geometry wrapper")
    unittest.main()
