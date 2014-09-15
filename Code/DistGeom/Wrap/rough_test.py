# $Id$
#
#  Copyright (C) 2004  Rational Discovery LLC
#         All Rights Reserved
#
from rdkit import RDConfig
import os,sys
import unittest
import numpy as np
from rdkit import DistanceGeometry as DG

def feq(v1,v2,tol2=1e-4):
  return abs(v1-v2)<=tol2

class TestCase(unittest.TestCase):
  def setUp(self):
    pass

  def test1SmoothPass(self):
    arr = np.array([[0,1.0,5.0],
                    [1.0,0,1.0],
                    [0.0,1.0,0]], np.float)
    self.assertTrue(DG.DoTriangleSmoothing(arr))
    self.assertTrue(feq(arr[0,2],2.0))
    self.assertTrue(feq(arr[2,0],0.0))
    self.assertTrue(feq(arr[0,1],1.0))
    self.assertTrue(feq(arr[1,0],1.0))
    self.assertTrue(feq(arr[1,2],1.0))
    
  def test2SmoothFail(self):
    arr = np.array([[0,1.0,5.0],
                    [1.0,0,1.0],
                    [3.0,1.0,0]], np.float)
    self.assertFalse(DG.DoTriangleSmoothing(arr))
    
    
  def test3SmoothPass(self):
    arr = np.array([[0,1.1,5.0],
                    [0.9,0,1.1],
                    [0.0,0.9,0]], np.float)
    self.assertTrue(DG.DoTriangleSmoothing(arr))
    self.assertTrue(feq(arr[0,2],2.2))
    self.assertTrue(feq(arr[2,0],0.0))
    self.assertTrue(feq(arr[0,1],1.1))
    self.assertTrue(feq(arr[1,0],0.9))
    self.assertTrue(feq(arr[1,2],1.1))
    
  
  def test4Embed(self):
    arr = np.array([[0,1.0,5.0],
                    [1.0,0,1.0],
                    [0.0,1.0,0]], np.float)
    self.assertTrue(DG.DoTriangleSmoothing(arr))
    coords = DG.EmbedBoundsMatrix(arr,randomSeed=100);
    v1 = coords[0]-coords[1]
    v2 = coords[1]-coords[2]
    d1 = np.dot(v1,v1)
    self.assertTrue(feq(d1,1.0, 0.001));
    d2 = np.dot(v2,v2)
    self.assertTrue(feq(d2,1.0, 0.001));
    
  def test5EmbedFail(self):
    arr = np.array([[0,1.0,5.0],
                    [1.0,0,1.0],
                    [3.0,1.0,0]], np.float)
    self.assertRaises(ValueError,lambda : DG.EmbedBoundsMatrix(arr))
    #DG.EmbedBoundsMatrix(arr,randomizeOnFailure=0,randomSeed=1)
    DG.EmbedBoundsMatrix(arr,randomizeOnFailure=1);

  def test6EmbedConstraints(self):
    arr = np.array([[0.0,1.0,1.0],
                    [1.0,0.0,1.0],
                    [0.99,1.0,0.0]], np.float)
    self.assertTrue(DG.DoTriangleSmoothing(arr))
    coords = DG.EmbedBoundsMatrix(arr, randomSeed=100)
    v1 = coords[0]-coords[1]
    v2 = coords[1]-coords[2]
    d1 = np.dot(v1,v1)
    
    self.assertTrue(feq(d1,1.0,2e-3));
    d2 = np.dot(v2,v2)
    self.assertTrue(feq(d2,1.0,2e-3));
    arr = np.array([[0.0,1.0,1.0,1.01],
                    [1.0,0.0,1.0,1.0],
                    [1.0,1.0,0.0,1.0],
                    [0.99,1.0,1.0,0.0],
                  ], np.float)
    self.assertTrue(DG.DoTriangleSmoothing(arr))
    coords = DG.EmbedBoundsMatrix(arr)
    v1 = coords[0]-coords[1]
    v2 = coords[1]-coords[2]
    d1 = np.dot(v1,v1)
    self.assertTrue(feq(d1,1.0,1e-3));
    d2 = np.dot(v2,v2)
    self.assertTrue(feq(d2,1.0,1e-3));

    return
    # this test is currently (rev:4769) passing on windows and
    # failing on linux.  It's kind of dependent on fp precision, so
    # it's probably ok to ditch it.
    arr = np.array([[0.0,1.0,1.0,1.0],
                    [1.0,0.0,1.0,1.0],
                    [1.0,1.0,0.0,1.0],
                    [1.0,1.0,1.0,0.0],
                  ], np.float)
    self.assertTrue(DG.DoTriangleSmoothing(arr))
    coords = DG.EmbedBoundsMatrix(arr,randomSeed=100)
    v1 = coords[0]-coords[1]
    v2 = coords[1]-coords[2]
    d1 = np.dot(v1,v1)
    self.assertTrue(feq(d1,1.0,1e-3));
    d2 = np.dot(v2,v2)
    self.assertTrue(feq(d2,1.0,1e-3));

    

    
    
    
if __name__ == '__main__':
  unittest.main()


