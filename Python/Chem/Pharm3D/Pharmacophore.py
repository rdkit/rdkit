# $Id: Pharmacophore.py 5083 2006-03-11 17:56:19Z NightlyBuild $
#
# Copyright (C) 2004-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
import Geometry
from Numeric import *
import Chem
from Chem import ChemicalFeatures

class Pharmacophore:
  def __init__(self, feats,initMats=True):
    self._initializeFeats(feats)
    nf = len(feats)
    self._boundsMat = zeros((nf, nf), Float)
    self._boundsMat2D = zeros((nf, nf), Int)
    if initMats:
      self._initializeMatrices()

  def _initializeFeats(self,feats):
    self._feats = []
    for feat in feats:
      if isinstance(feat,ChemicalFeatures.MolChemicalFeature):
        pos = feat.GetPos()
        newFeat = ChemicalFeatures.FreeChemicalFeature(feat.GetFamily(),feat.GetType(),
                                           Geometry.Point3D(pos[0],pos[1],pos[2]))
        self._feats.append(newFeat)
      else:
        self._feats.append(feat)

  def _initializeMatrices(self):
    # initialize the bounds matrix with distances to start with
    nf = len(self._feats)
    for i in range(1, nf):
      loci = self._feats[i].GetPos()
      for j in range(i):
        locj = self._feats[j].GetPos()
        dist = loci.Distance(locj)
        self._boundsMat[i,j] = dist
        self._boundsMat[j,i] = dist
    for i in range(nf):
      for j in range(i+1,nf):
        self._boundsMat2D[i,j] = 1000


  def getFeatures(self):
    return self._feats

  def getFeature(self, i):
    return self._feats[i]

  def getUpperBound(self, i, j):
    if (i > j):
      j,i = i,j
    return self._boundsMat[i, j]

  def getLowerBound(self, i, j):
    if (j > i):
      j,i = i,j
    return self._boundsMat[i,j]

  def _checkBounds(self,i,j):
    " raises ValueError on failure "
    nf = len(self._feats)
    if (i < 0) or (i >= nf):
      raise ValueError, "Index out of bound"
    if (j < 0) or (j >= nf):
      raise ValueError, "Index out of bound"
    return True

  def setUpperBound(self, i, j, val, checkBounds=False):
    if (checkBounds): self._checkBounds(i,j)
    if (i > j):
      j,i = i,j
    self._boundsMat[i,j] = val

  def setLowerBound(self, i, j, val, checkBounds=False):
    if (checkBounds): self._checkBounds(i,j)
    if (j > i):
      j,i = i,j
    self._boundsMat[i,j] = val

  def getUpperBound2D(self, i, j):
    if (i > j):
      j,i = i,j
    return self._boundsMat2D[i, j]

  def getLowerBound2D(self, i, j):
    if (j > i):
      j,i = i,j
    return self._boundsMat2D[i,j]


  def setUpperBound2D(self, i, j, val, checkBounds=False):
    if (checkBounds): self._checkBounds(i,j)
    if (i > j):
      j,i = i,j
    self._boundsMat2D[i,j] = val

  def setLowerBound2D(self, i, j, val, checkBounds=False):
    if (checkBounds): self._checkBounds(i,j)
    if (j > i):
      j,i = i,j
    self._boundsMat2D[i,j] = val

  def __str__(self):
    res = ' '*13
    for i,iFeat in enumerate(self._feats):
      res += '% 12s '%iFeat.GetFamily()
    res += '\n'
    for i,iFeat in enumerate(self._feats):
      res += '% 12s '%iFeat.GetFamily()
      for j,jFeat in enumerate(self._feats):
        if j<i:
          res += '% 12.3f '%self.getLowerBound(i,j)
        elif j>i:
          res += '% 12.3f '%self.getUpperBound(i,j)
        else:
          res += '% 12.3f '%0.0
      res += '\n'
    return res
      
