#
# Copyright (C) 2016 Greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.

# generates reference data for the PMI descriptors
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from numpy import linalg

def GetMoments(mol, includeWeights):
  conf = mol.GetConformer()
  if includeWeights:
    weights = [x.GetMass() for x in mol.GetAtoms()]
  else:
    weights = [1.0]*mol.GetNumAtoms()

  pts = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
  wSum = sum(weights)
  origin = np.sum([pts[i]*weights[i] for i in range(mol.GetNumAtoms())], 0)
  origin /= wSum
  sumXX = 0
  sumXY = 0
  sumXZ = 0
  sumYY = 0
  sumYZ = 0
  sumZZ = 0
  sums = np.zeros((3, 3), np.double)
  for j,pt in enumerate(pts):
    dp = weights[j]*(pt - origin)
    for i in range(3):
      sums[i, i] += dp[i] * dp[i]
      for j in range(i + 1, 3):
        sums[i, j] += dp[i] * dp[j]
        sums[j, i] += dp[i] * dp[j]
  sums /= wSum
  vals, vects = linalg.eigh(sums)
  vals = sorted(vals)
  return vals

if __name__ == '__main__':
  suppl = Chem.SDMolSupplier('./PBF_egfr.sdf', removeHs=False)
  output = open('./PMI_egfr.out', 'w+')
  for m in suppl:
    i1,i2,i3 = GetMoments(m,True)
    mi1,mi2,mi3 = GetMoments(m,False)
    print("%s %.4f %.4f %.4f %.4f %.4f %.4f"%(m.GetProp("_Name"),i1,i2,i3,mi1,mi2,mi3),file=output)
