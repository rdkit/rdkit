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


def GetMoments(mol, includeWeights):
  conf = mol.GetConformer()
  if includeWeights:
    masses = np.array([x.GetMass() for x in mol.GetAtoms()])
  else:
    masses = [1.0] * mol.GetNumAtoms()

  ps = conf.GetPositions()
  mps = [x*y for x,y in zip(ps,masses)]
  centroid = np.sum(mps,axis=0)/sum(masses)
  cps = ps - centroid
  xx = xy = xz = yy = yz = zz = 0.0
  for m,p in zip(masses,cps):
    xx += m*(p[1]*p[1] + p[2]*p[2])
    yy += m*(p[0]*p[0] + p[2]*p[2])
    zz += m*(p[1]*p[1] + p[0]*p[0])
    xy -= m*p[0]*p[1]
    xz -= m*p[0]*p[2]
    yz -= m*p[1]*p[2]
  covm = np.array([[xx,xy,xz],[xy,yy,yz],[xz,yz,zz]])
  res = np.linalg.eigvals(covm)
  return sorted(res)


if __name__ == '__main__':
  suppl = Chem.SDMolSupplier('./PBF_egfr.sdf', removeHs=False)
  output = open('./PMI_egfr.out', 'w+')
  for m in suppl:
    i1, i2, i3 = GetMoments(m, True)
    mi1, mi2, mi3 = GetMoments(m, False)
    print("%s %.4f %.4f %.4f %.4f %.4f %.4f" % (m.GetProp("_Name"), i1, i2, i3, mi1, mi2, mi3),
          file=output)
