# $Id$
#
#  Copyright (C) 2002-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""unit testing code for the signature utils

"""
import unittest
from rdkit.Chem.Pharm2D import Utils

class TestCase(unittest.TestCase):
  def testCounts(self):
    vals = [
      ((0,1,2),4),
      ((0,0,0),0),
      ((2,2,2),9),
      ((1,1,2),7),
      ((1,2,2),8),
      ]
    for combo,tgt in vals:
      res = Utils.CountUpTo(3,3,combo)
      assert res==tgt,'Bad res (%d) for combo %s'%(res,str((combo,tgt)))

  def testGetTriangles(self):
    vals = [
      (2,0,[]),
      (3,1,((0,1,2),)),
      (4,2,((0,1,3),(1,2,4))),
      (5,3,((0,1,4),(1,2,5),(2,3,6))),
      ]
    for tpl in vals:
      nPts,cnt,tris = tpl
      r = Utils.GetTriangles(nPts)
      assert len(r)==cnt,'bad triangle length %d for probe %s'%(len(r),str(tpl))
      assert r==tris,'bad triangle list %s for probe %s'%(str(r),str(tpl))
      
  def testDistTriangleInequality(self):
    bins = [(1,2),(2,3),(5,6)]
    vals = [
      ((0,0,0),1),
      ((0,0,1),1),
      ((0,0,2),0),
      ((0,1,2),1),
      ((1,1,2),1),
      ]
    for tpl in vals:
      ds,tgt = tpl
      distBins = [bins[x] for x in ds]
      r = Utils.BinsTriangleInequality(distBins[0],distBins[1],distBins[2])
      assert r==tgt,'bad result %d for probe %s'%(r,str(tpl))

  def testLimitPharmacophores(self):
    bins = [(1,2),(2,3),(5,6)]
    vals = [
      ((0,0,0),1),
      ((0,0,1),1),
      ((0,0,2),0),
      ((0,1,2),1),
      ((1,1,2),1),
      ]
    for tpl in vals:
      ds,tgt = tpl
      r = Utils.ScaffoldPasses(ds,bins)
      assert r==tgt,'bad result %d for probe %s'%(r,str(tpl))

  def testGetPossiblePharmacophores(self):
    bins = [(1,2),(2,3),(5,6)]
    vals = [
      (2,3),
      (3,24),
      ]
    for tpl in vals:
      num,tgt = tpl
      pphores = Utils.GetPossibleScaffolds(num,bins)
      cnt = len(pphores)
      assert cnt==tgt,'bad pharmacophore count %d for probe %s'%(cnt,str(tpl))

    
if __name__ == '__main__':
  unittest.main()

