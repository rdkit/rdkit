#
# Copyright (C) 2004-2008 Greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import numpy

from rdkit import Geometry
from rdkit.Chem import ChemicalFeatures
from rdkit.RDLogger import logger

logger = logger()


class Pharmacophore:

  def __init__(self, feats, initMats=True):
    self._initializeFeats(feats)
    nf = len(feats)
    self._boundsMat = numpy.zeros((nf, nf), numpy.float)
    self._boundsMat2D = numpy.zeros((nf, nf), numpy.int)
    if initMats:
      self._initializeMatrices()

  def _initializeFeats(self, feats):
    self._feats = []
    for feat in feats:
      if isinstance(feat, ChemicalFeatures.MolChemicalFeature):
        newFeat = ChemicalFeatures.FreeChemicalFeature(
          feat.GetFamily(), feat.GetType(), feat.GetPos())
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
        self._boundsMat[i, j] = dist
        self._boundsMat[j, i] = dist
    for i in range(nf):
      for j in range(i + 1, nf):
        self._boundsMat2D[i, j] = 1000

  def getFeatures(self):
    return self._feats

  def getFeature(self, i):
    return self._feats[i]

  def getUpperBound(self, i, j):
    if i > j:
      j, i = i, j
    return self._boundsMat[i, j]

  def getLowerBound(self, i, j):
    if j > i:
      j, i = i, j
    return self._boundsMat[i, j]

  def _checkBounds(self, i, j):
    " raises ValueError on failure "
    nf = len(self._feats)
    if (0 <= i < nf) and (0 <= j < nf):
      return True
    raise ValueError("Index out of bound")

  def setUpperBound(self, i, j, val, checkBounds=False):
    if checkBounds:
      self._checkBounds(i, j)
    if i > j:
      j, i = i, j
    self._boundsMat[i, j] = val

  def setLowerBound(self, i, j, val, checkBounds=False):
    if checkBounds:
      self._checkBounds(i, j)
    if j > i:
      j, i = i, j
    self._boundsMat[i, j] = val

  def getUpperBound2D(self, i, j):
    if i > j:
      j, i = i, j
    return self._boundsMat2D[i, j]

  def getLowerBound2D(self, i, j):
    if j > i:
      j, i = i, j
    return self._boundsMat2D[i, j]

  def setUpperBound2D(self, i, j, val, checkBounds=False):
    if checkBounds:
      self._checkBounds(i, j)
    if i > j:
      j, i = i, j
    self._boundsMat2D[i, j] = val

  def setLowerBound2D(self, i, j, val, checkBounds=False):
    if checkBounds:
      self._checkBounds(i, j)
    if j > i:
      j, i = i, j
    self._boundsMat2D[i, j] = val

  def __str__(self):
    res = ' ' * 14
    for i, iFeat in enumerate(self._feats):
      res += '%13s ' % iFeat.GetFamily()
    res += '\n'
    for i, iFeat in enumerate(self._feats):
      res += '%13s ' % iFeat.GetFamily()
      for j, _ in enumerate(self._feats):
        if j < i:
          res += '%13.3f ' % self.getLowerBound(i, j)
        elif j > i:
          res += '%13.3f ' % self.getUpperBound(i, j)
        else:
          res += '% 13.3f ' % 0.0
      res += '\n'
    return res


class ExplicitPharmacophore(Pharmacophore):
  """ this is a pharmacophore with explicit point locations and radii
  """

  def __init__(self, feats=None, radii=None, inputData=None, initMats=True):
    if feats is not None and radii is not None:
      self._initializeFeats(feats, radii)
    elif isinstance(inputData, file):
      self.initFromFile(inputData)
    elif isinstance(inputData, str):
      self.initFromString(inputData)
    else:
      raise TypeError("At least one of (feats, radii) and inputData must be specified.")

    nf = len(self._feats)
    self._boundsMat = numpy.zeros((nf, nf), numpy.float)
    self._boundsMat2D = numpy.zeros((nf, nf), numpy.int)
    if initMats:
      self._initializeMatrices()
      self.applyRadiiToBounds()

  def _initializeFeats(self, feats, radii):
    if len(feats) != len(radii):
      raise ValueError('len(feats)!=len(radii)')
    self._feats = []
    self._radii = []
    for feat, rad in zip(feats, radii):
      if isinstance(feat, ChemicalFeatures.MolChemicalFeature):
        newFeat = ChemicalFeatures.FreeChemicalFeature(
            feat.GetFamily(), feat.GetType(), feat.GetPos())
      else:
        newFeat = feat
      self._feats.append(newFeat)
      self._radii.append(rad)

  def getRadii(self):
    return self._radii

  def getRadius(self, i):
    return self._radii[i]

  def setRadius(self, i, rad):
    self._radii[i] = rad

  def applyRadiiToBounds(self):
    """
    From Nikolaus Stiefl
    https://github.com/rdkit/UGM_2016/blob/master/Notebooks/Stiefl_
      RDKitPh4FullPublication.ipynb
    """
    self._initializeMatrices()
    radii = self.getRadii()
    for i in range(len(radii)):
      for j in range(i+1, len(radii)):
        sumRadii = radii[i] + radii[j]
        self.setLowerBound(i, j, max(self.getLowerBound(i,j) - sumRadii, 0))
        self.setUpperBound(i, j, self.getUpperBound(i, j) + sumRadii)

  def initFromString(self, text):
    lines = text.split(r'\n')
    self.initFromLines(lines)

  def initFromFile(self, inF):
    self.initFromLines(inF.readlines())

  def initFromLines(self, lines):
    import re
    spaces = re.compile('[\ \t]+')

    feats = []
    rads = []
    for lineNum, line in enumerate(lines):
      txt = line.split('#')[0].strip()
      if txt:
        splitL = spaces.split(txt)
        if len(splitL) < 5:
          logger.error('Input line %d only contains %d fields, 5 are required. Read failed.' %
                       (lineNum, len(splitL)))
          return
        fName = splitL[0]
        try:
          xP = float(splitL[1])
          yP = float(splitL[2])
          zP = float(splitL[3])
          rad = float(splitL[4])
        except ValueError:
          logger.error('Error parsing a number of line %d. Read failed.' % (lineNum))
          return
        feats.append(
          ChemicalFeatures.FreeChemicalFeature(fName, fName, Geometry.Point3D(xP, yP, zP)))
        rads.append(rad)
    self._initializeFeats(feats, rads)

  def getFeatureStr(self):
    """
    Returns a representation of the features that constitute the Pharmacophore.
    The __str__ method, on the other hand, returns a representation of the
    bounds matrix.
    """
    res = ' ' * 14
    res += '%13s %13s %13s ' % ('x', 'y', 'z')
    res += '%13s' % "Radius"
    res += '\n'
    for feat, rad in zip(self._feats, self._radii):
      res += '% 13s ' % feat.GetFamily()
      p = feat.GetPos()
      res += '% 13.3f % 13.3f % 13.3f ' % (p.x, p.y, p.z)
      res += '% 13.3f' % rad
      res += '\n'
    return res