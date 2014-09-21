## Automatically adapted for numpy.oldnumeric Jun 27, 2008 by -c
from __future__ import print_function

import rdkit.Numerics.rdAlignment as rdAlg
from rdkit import Geometry
from rdkit import RDConfig
import os,sys
import unittest
import numpy as np
import math

def lstFeq(l1, l2, tol=1.e-4):
  if (len(list(l1)) != len(list(l2))):
    return 0
  for i in range(len(list(l1))):
    if not feq(l1[i], l2[i], tol):
      return 0
  return 1

def feq(v1,v2,tol2=1e-4):
  return abs(v1-v2)<=tol2

def transformPoint(trans, pt):
  pt2 = list(pt)
  pt2.append(1.0)
  pt2 = np.array(pt2)
  res = np.dot(trans, pt2)
  return res[:3]

class TestCase(unittest.TestCase):
    def setUp(self):
        pass

    def test1Basic(self):
        # passing two numeric arrays
        refPts = np.zeros((2,3), np.float)
        prbPts = np.zeros((2,3), np.float)

        refPts[1,0] = 1.0
                
        prbPts[0,0] = 2.0
        prbPts[0,1] = 2.0
        prbPts[1,0] = 2.0
        prbPts[1,1] = 3.0
                  
        res = rdAlg.GetAlignmentTransform(refPts, prbPts)
        self.assertTrue(feq(res[0], 0.0))
        refLst = list(refPts)
        cnt = 0
        for item in list(prbPts):
          self.assertTrue(lstFeq(transformPoint(res[1], item), refLst[cnt]))
          cnt+= 1
          
        # repeat with with lists or tuples
        refPts = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0))
        prbPts = ((2.0, 2.0, 0.0), (2.0, 3.0, 0.0))
        res = rdAlg.GetAlignmentTransform(refPts, prbPts)
        self.assertTrue(feq(res[0], 0.0))

        refPts = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
        prbPts = [[2.0, 2.0, 0.0], [2.0, 3.0, 0.0]]
        res = rdAlg.GetAlignmentTransform(refPts, prbPts)
        self.assertTrue(feq(res[0], 0.0))

        # mix it up
        refPts = np.zeros((2,3), np.float)
        refPts[1,0] = 1.0
        res = rdAlg.GetAlignmentTransform(refPts, prbPts)
        self.assertTrue(feq(res[0], 0.0))
                
    def test2Weights(self) :
        refPts = np.array([[-math.cos(math.pi/6), -math.sin(math.pi/6), 0.0],
                           [math.cos(math.pi/6), -math.sin(math.pi/6), 0.0],
                           [0.0, 1.0, 0.0]], np.float)
        prbPts = np.array([[-2*math.sin(math.pi/6) + 3.0, 2*math.cos(math.pi/6), 4.0],
                           [-2*math.sin(math.pi/6) + 3.0, -2*math.cos(math.pi/6), 4.0],
                           [5.0, 0.0, 4.0]], np.float)
        res = rdAlg.GetAlignmentTransform(refPts, prbPts)
        self.assertTrue(feq(res[0], 3.0))
        target = [[-1.732, -1., 0.],
                  [1.732, -1., 0.],
                  [0., 2., 0.]]
        cnt = 0
        for item in list(prbPts):
          self.assertTrue(lstFeq(transformPoint(res[1], item), target[cnt]))
          cnt += 1
          
        weights = np.array([1.0, 1.0, 2.0], np.float)
        res = rdAlg.GetAlignmentTransform(refPts, prbPts, weights)
        self.assertTrue(feq(res[0], 3.75))
        cnt = 0
        target = [[-1.732, -1.25, 0.],
                  [1.732, -1.25, 0.],
                  [0., 1.75, 0.]]
        for item in list(prbPts):
          self.assertTrue(lstFeq(transformPoint(res[1], item), target[cnt]))
          cnt += 1
        weights = [1.0, 1.0, 2.0]
        res = rdAlg.GetAlignmentTransform(refPts, prbPts, weights)
        self.assertTrue(feq(res[0], 3.75))

        weights = [1.0, 2.0, 2.0]
        res = rdAlg.GetAlignmentTransform(refPts, prbPts, weights)
        self.assertTrue(feq(res[0], 4.8))

    def test3tetra(self) :
        refPts = np.array([[0.0, 0.0, 0.0],
                           [1.0, 0.0, 0.0],
                           [0.0, 1.0, 0.0],
                           [0.0, 0.0, 1.0]], np.float)
        prbPts = np.array([[2.0, 2.0, 3.0],
                           [3.0, 2.0, 3.0],
                           [2.0, 3.0, 3.0]], np.float)
        self.assertRaises(ValueError,lambda : rdAlg.GetAlignmentTransform(refPts, prbPts))

        prbPts = np.array([[2.0, 2.0, 3.0],
                           [3.0, 2.0, 3.0],
                           [2.0, 3.0, 3.0], 
                           [2.0, 2.0, 4.0]], np.float)
        res = rdAlg.GetAlignmentTransform(refPts, prbPts)
        self.assertTrue(feq(res[0], 0.0))

        wts = [1.0, 1.0, 1.0]
        self.assertRaises(ValueError,lambda : rdAlg.GetAlignmentTransform(refPts, prbPts, wts))

        wts = [1.0, 1.0, 1.0, 1.0]
        res = rdAlg.GetAlignmentTransform(refPts, prbPts, wts)
        self.assertTrue(feq(res[0], 0.0))
        
        # test reflection
        prbPts = np.array([[2.0, 2.0, 3.0],
                           [3.0, 2.0, 3.0],
                           [2.0, 2.0, 4.0], 
                           [2.0, 3.0, 3.0]], np.float)
        res = rdAlg.GetAlignmentTransform(refPts, prbPts, wts)
        self.assertTrue(feq(res[0], 1.0))
        
        res = rdAlg.GetAlignmentTransform(refPts, prbPts, wts, 1)
        self.assertTrue(feq(res[0], 0.0))
        cnt = 0
        refLst = list(refPts)
        for item in list(prbPts):
          self.assertTrue(lstFeq(transformPoint(res[1], item), refLst[cnt]))
          cnt += 1

    def test4points(self) :
        refPts = (Geometry.Point3D(0.0, 0.0, 0.0),
                  Geometry.Point3D(1.0, 0.0, 0.0),
                  Geometry.Point3D(0.0, 1.0, 0.0),
                  Geometry.Point3D(0.0, 0.0, 1.0),
                  )
        prbPts = (Geometry.Point3D(2.0, 2.0, 3.0),
                  Geometry.Point3D(3.0, 2.0, 3.0),
                  Geometry.Point3D(2.0, 3.0, 3.0),
                  Geometry.Point3D(2.0, 2.0, 4.0),
                  )

        res = rdAlg.GetAlignmentTransform(refPts, prbPts)
        self.assertTrue(feq(res[0], 0.0))

    def test5errorHandling(self) :
        refPts = (Geometry.Point3D(0.0, 0.0, 0.0),
                  Geometry.Point3D(1.0, 0.0, 0.0),
                  Geometry.Point3D(0.0, 1.0, 0.0),
                  Geometry.Point3D(0.0, 0.0, 1.0),
                  )
        prbPts = (1,2,3,4,)
        self.assertRaises(ValueError,
                          lambda : rdAlg.GetAlignmentTransform(refPts, prbPts))
        prbPts = ()
        self.assertRaises(ValueError,
                          lambda : rdAlg.GetAlignmentTransform(refPts, prbPts))

        prbPts = 1
        self.assertRaises(ValueError,
                          lambda : rdAlg.GetAlignmentTransform(refPts, prbPts))

        prbPts = (Geometry.Point3D(2.0, 2.0, 3.0),
                  Geometry.Point3D(3.0, 2.0, 3.0),
                  Geometry.Point3D(2.0, 3.0, 3.0),
                  (2.0, 2.0, 5.0),
                  )
        self.assertRaises(ValueError,
                          lambda : rdAlg.GetAlignmentTransform(refPts, prbPts))


        
if __name__ == '__main__':
  print("Testing Alignment Wrapper code:")
  unittest.main()
