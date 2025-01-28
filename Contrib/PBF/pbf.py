#
# Copyright (C) 2015 Greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.

import numpy as np
from numpy import linalg

from rdkit import Chem
from rdkit.Chem import AllChem


def GetBestFitPlane(pts, weights=None):
  if weights is None:
    wSum = len(pts)
    origin = np.sum(pts, 0)
  origin /= wSum

  sumXX = 0
  sumXY = 0
  sumXZ = 0
  sumYY = 0
  sumYZ = 0
  sumZZ = 0
  sums = np.zeros((3, 3), np.double)
  for pt in pts:
    dp = pt - origin
    for i in range(3):
      sums[i, i] += dp[i] * dp[i]
      for j in range(i + 1, 3):
        sums[i, j] += dp[i] * dp[j]
        sums[j, i] += dp[i] * dp[j]
  sums /= wSum
  vals, vects = linalg.eigh(sums)
  order = np.argsort(vals)
  normal = vects[:, order[0]]
  plane = np.zeros((4, ), np.double)
  plane[:3] = normal
  plane[3] = -1 * normal.dot(origin)
  return plane


def PBFRD(mol, confId=-1):
  conf = mol.GetConformer(confId)
  if not conf.Is3D():
    return 0

  pts = np.array([list(conf.GetAtomPosition(x)) for x in range(mol.GetNumAtoms())])
  plane = GetBestFitPlane(pts)
  denom = np.dot(plane[:3], plane[:3])
  denom = denom**0.5
  # add up the distance from the plane for each point:
  res = 0.0
  for pt in pts:
    res += np.abs(pt.dot(plane[:3]) + plane[3])
  res /= denom
  res /= len(pts)
  return res


if __name__ == '__main__':
  suppl = Chem.SDMolSupplier('./testData/egfr.sdf', removeHs=False)
  expected = open('./testData/egfr.out', 'r')
  for m in suppl:
    res = PBFRD(m)
    inl = next(expected).strip().split()
    expect = float(inl[1])
    assert abs(res - expect) < 1e-4
