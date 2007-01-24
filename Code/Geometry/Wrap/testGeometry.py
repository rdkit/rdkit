import RDConfig
import os,sys
import unittest
import DataStructs
from Geometry import rdGeometry as geom
import cPickle
import math

def feq(v1, v2, tol=1.0e-4):
    return abs(v1-v2) < tol

class TestCase(unittest.TestCase):
    def setUp(self) :
        pass

    def test1Point3D(self):
        pt = geom.Point3D();
        self.failUnless(feq(pt.x, 0.0))
        self.failUnless(feq(pt.y, 0.0))
        self.failUnless(feq(pt.z, 0.0))

        pt = geom.Point3D(3., 4., 5.)
        self.failUnless(feq(pt.x, 3.0))
        self.failUnless(feq(pt.y, 4.0))
        self.failUnless(feq(pt.z, 5.0))
        self.failUnless(feq(pt[0], 3.0))
        self.failUnless(feq(pt[1], 4.0))
        self.failUnless(feq(pt[2], 5.0))
        self.failUnless(feq(pt[-3], 3.0))
        self.failUnless(feq(pt[-2], 4.0))
        self.failUnless(feq(pt[-1], 5.0))
        lst = list(pt)
        self.failUnless(feq(lst[0], 3.0))
        self.failUnless(feq(lst[1], 4.0))
        self.failUnless(feq(lst[2], 5.0))
        
        pt2 = geom.Point3D(1., 1., 1.)

        pt3 = pt+pt2
        self.failUnless(feq(pt3.x, 4.0))
        self.failUnless(feq(pt3.y, 5.0))
        self.failUnless(feq(pt3.z, 6.0))
        
        
        pt += pt2
        self.failUnless(feq(pt.x, 4.0))
        self.failUnless(feq(pt.y, 5.0))
        self.failUnless(feq(pt.z, 6.0))

        pt3 = pt-pt2
        self.failUnless(feq(pt3.x, 3.0))
        self.failUnless(feq(pt3.y, 4.0))
        self.failUnless(feq(pt3.z, 5.0))
        
        pt -= pt2
        self.failUnless(feq(pt.x, 3.0))
        self.failUnless(feq(pt.y, 4.0))
        self.failUnless(feq(pt.z, 5.0))
        
        pt *= 2.0
        self.failUnless(feq(pt.x, 6.0))
        self.failUnless(feq(pt.y, 8.0))
        self.failUnless(feq(pt.z, 10.0))

        pt /= 2
        self.failUnless(feq(pt.x, 3.0))
        self.failUnless(feq(pt.y, 4.0))
        self.failUnless(feq(pt.z, 5.0))
        self.failUnless(feq(pt.Length(), 7.0711))
        self.failUnless(feq(pt.LengthSq(), 50.0))
        pt.Normalize()
        self.failUnless(feq(pt.Length(), 1.0))

        pt1 = geom.Point3D(1.0, 0.0, 0.0)
        pt2 = geom.Point3D(2.0*math.cos(math.pi/6), 2.0*math.sin(math.pi/6), 0.0)
        ang = pt1.AngleTo(pt2)
        self.failUnless(feq(ang, math.pi/6))

        prod = pt1.DotProduct(pt2)
        self.failUnless(feq(prod, 2.0*math.cos(math.pi/6)))

        pt3 = pt1.CrossProduct(pt2)
        self.failUnless(feq(pt3.x, 0.0))
        self.failUnless(feq(pt3.y, 0.0))
        self.failUnless(feq(pt3.z, 1.0))



    def test1Point2D(self):
        pt = geom.Point2D();
        self.failUnless(feq(pt.x, 0.0))
        self.failUnless(feq(pt.y, 0.0))
        
        pt = geom.Point2D(3., 4.)
        self.failUnless(feq(pt.x, 3.0))
        self.failUnless(feq(pt.y, 4.0))
        
        self.failUnless(feq(pt.x, 3.0))
        self.failUnless(feq(pt.y, 4.0))
        self.failUnless(feq(pt[0], 3.0))
        self.failUnless(feq(pt[1], 4.0))
        self.failUnless(feq(pt[-2], 3.0))
        self.failUnless(feq(pt[-1], 4.0))
        lst = list(pt)
        self.failUnless(feq(lst[0], 3.0))
        self.failUnless(feq(lst[1], 4.0))


        pt2 = geom.Point2D(1., 1.)

        pt3 = pt+pt2
        self.failUnless(feq(pt3.x, 4.0))
        self.failUnless(feq(pt3.y, 5.0))
        
        pt += pt2
        self.failUnless(feq(pt.x, 4.0))
        self.failUnless(feq(pt.y, 5.0))
        
        pt3 = pt-pt2
        self.failUnless(feq(pt3.x, 3.0))
        self.failUnless(feq(pt3.y, 4.0))
        
        pt -= pt2
        self.failUnless(feq(pt.x, 3.0))
        self.failUnless(feq(pt.y, 4.0))
                
        pt *= 2.0
        self.failUnless(feq(pt.x, 6.0))
        self.failUnless(feq(pt.y, 8.0))
        
        pt /= 2
        self.failUnless(feq(pt.x, 3.0))
        self.failUnless(feq(pt.y, 4.0))
        self.failUnless(feq(pt.Length(), 5.0))
        self.failUnless(feq(pt.LengthSq(), 25.0))
        pt.Normalize()
        self.failUnless(feq(pt.Length(), 1.0))

        pt1 = geom.Point2D(1.0, 0.0)
        pt2 = geom.Point2D(2.0*math.cos(math.pi/6), 2.0*math.sin(math.pi/6))
        ang = pt1.AngleTo(pt2)
        self.failUnless(feq(ang, math.pi/6))

        prod = pt1.DotProduct(pt2)
        self.failUnless(feq(prod, 2.0*math.cos(math.pi/6)))

    def test3UniformGrid(self):
        ugrid = geom.UniformGrid3D(20, 18, 15)
        self.failUnless(ugrid.GetNumX() == 40)
        self.failUnless(ugrid.GetNumY() == 36)
        self.failUnless(ugrid.GetNumZ() == 30)
        dvect = ugrid.GetOccupancyVect()
        ugrid = geom.UniformGrid3D(20, 18, 15, 0.5, DataStructs.DiscreteValueType.TWOBITVALUE)
        dvect = ugrid.GetOccupancyVect()
        self.failUnless(dvect.GetValueType() == DataStructs.DiscreteValueType.TWOBITVALUE)

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
        self.failUnless(dist == 0.25)

        grd2 = geom.UniformGrid3D(10.0, 10.0, 10.0, 0.5, DataStructs.DiscreteValueType.FOURBITVALUE)
        grd2.SetSphereOccupancy(geom.Point3D(-2.0, -2.0, 0.0), 1.5, 0.25, 3)
        grd2.SetSphereOccupancy(geom.Point3D(-2.0, 2.0, 0.0), 1.5, 0.25, 3)
        grd2.SetSphereOccupancy(geom.Point3D(2.0, -2.0, 0.0), 1.5, 0.25, 3)
        self.failUnlessRaises(ValueError, lambda : geom.TanimotoDistance(grd, grd2))

        grd2 = geom.UniformGrid3D(10.0, 10.0, 10.0, 1.0)
        self.failUnlessRaises(ValueError, lambda : geom.TanimotoDistance(grd, grd2))

        grd2 = geom.UniformGrid3D(11.0, 10.0, 10.0, 1.0)
        self.failUnlessRaises(ValueError, lambda : geom.TanimotoDistance(grd, grd2))

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
                     self.failUnless(grd.GetValPoint(bPt1) == grd.GetValPoint(bPt2))
                 
                 bPt1.x = -4.0
                 bPt2.x = 4.0
            bPt1.y = -2.0
            bPt2.y = -2.0

    def testPointPickles(self):
        pt = geom.Point3D(2.0,-3.0,1.0)
        pt2 = cPickle.loads(cPickle.dumps(pt))
        self.failUnless(feq(pt.x,pt2.x,1e-6))
        self.failUnless(feq(pt.y,pt2.y,1e-6))
        self.failUnless(feq(pt.z,pt2.z,1e-6))

        pt = geom.Point2D(2.0,-4.0)
        pt2 = cPickle.loads(cPickle.dumps(pt))
        self.failUnless(feq(pt.x,pt2.x,1e-6))
        self.failUnless(feq(pt.y,pt2.y,1e-6))

if __name__=='__main__':
    print "Testing Geometry wrapper"
    unittest.main()
