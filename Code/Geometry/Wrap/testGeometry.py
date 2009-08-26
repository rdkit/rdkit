from rdkit import RDConfig
import os,sys
import unittest
from rdkit import DataStructs
from rdkit.Geometry import rdGeometry as geom
import cPickle,copy
import math

def feq(v1, v2, tol=1.0e-4):
    return abs(v1-v2) < tol

class TestCase(unittest.TestCase):
    def setUp(self) :
        pass

    def test1aPoint3D(self):
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



    def test1bPoint2D(self):
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

    def test1cPointND(self):
        dim=4
        pt = geom.PointND(4);
        for i in range(dim):
            self.failUnless(feq(pt[i], 0.0))
        
        pt[0]=3
        pt[3]=4
        self.failUnless(feq(pt[0], 3.0))
        self.failUnless(feq(pt[3], 4.0))
        self.failUnless(feq(pt[-4], 3.0))
        self.failUnless(feq(pt[-1], 4.0))
        lst = list(pt)
        self.failUnless(feq(lst[0], 3.0))
        self.failUnless(feq(lst[3], 4.0))


        pt2 = geom.PointND(4)
        pt2[0]=1.
        pt2[2]=1.

        pt3 = pt+pt2
        self.failUnless(feq(pt3[0], 4.0))
        self.failUnless(feq(pt3[2], 1.0))
        self.failUnless(feq(pt3[3], 4.0))
        
        pt += pt2
        self.failUnless(feq(pt[0], 4.0))
        self.failUnless(feq(pt[2], 1.0))
        self.failUnless(feq(pt[3], 4.0))

        pt3 = pt-pt2
        self.failUnless(feq(pt3[0], 3.0))
        self.failUnless(feq(pt3[2], 0.0))
        self.failUnless(feq(pt3[3], 4.0))
        
        pt -= pt2
        self.failUnless(feq(pt[0], 3.0))
        self.failUnless(feq(pt[2], 0.0))
        self.failUnless(feq(pt[3], 4.0))

        pt *= 2.0
        self.failUnless(feq(pt[0], 6.0))
        self.failUnless(feq(pt[1], 0.0))
        self.failUnless(feq(pt[2], 0.0))
        self.failUnless(feq(pt[3], 8.0))

        
        pt /= 2
        self.failUnless(feq(pt[0], 3.0))
        self.failUnless(feq(pt[1], 0.0))
        self.failUnless(feq(pt[2], 0.0))
        self.failUnless(feq(pt[3], 4.0))

        self.failUnless(feq(pt.Length(), 5.0))
        self.failUnless(feq(pt.LengthSq(), 25.0))
        pt.Normalize()
        self.failUnless(feq(pt.Length(), 1.0))

        pkl = cPickle.dumps(pt)
        pt2 = cPickle.loads(pkl)
        self.failUnless(len(pt)==len(pt2))
        for i in range(len(pt)):
            self.failUnless(feq(pt2[i],pt[i]))

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
