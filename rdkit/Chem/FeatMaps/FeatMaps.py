# $Id$
#
# Copyright (C) 2006 Greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import math

from rdkit.Chem.FeatMaps.FeatMapPoint import FeatMapPoint


class FeatMapScoreMode(object):
  All = 0
  """ score each feature in the probe against every matching
      feature in the FeatMap.
  """

  Closest = 1
  """ score each feature in the probe against the closest
      matching feature in the FeatMap.
  """

  Best = 2
  """ score each feature in the probe against the matching
      feature in the FeatMap that leads to the highest score
  """


class FeatDirScoreMode(object):
  Ignore = 0
  """ ignore feature directions
  """

  DotFullRange = 1
  """ Use the dot product and allow negative contributions when
      directions are anti-parallel.
      e.g. score = dot(f1Dir,f2Dir)
  """

  DotPosRange = 2
  """ Use the dot product and scale contributions to lie between
      zero and one.
      e.g. score = ( dot(f1Dir,f2Dir) + 1 ) / 2
  """


class FeatMapParams(object):
  """ one of these should be instantiated for each
  feature type in the feature map
  """
  radius = 2.5
  " cutoff radius "

  width = 1.0
  " width parameter (e.g. the gaussian sigma) "

  class FeatProfile(object):
    " scoring profile of the feature "
    Gaussian = 0
    Triangle = 1
    Box = 2

  featProfile = FeatProfile.Gaussian


class FeatMap(object):
  dirScoreMode = FeatDirScoreMode.Ignore
  scoreMode = FeatMapScoreMode.All
  params = {}

  def __init__(self, params=None, feats=None, weights=None):
    if params:
      self.params = params
    self._initializeFeats(feats, weights)

  def _initializeFeats(self, feats, weights):
    self._feats = []
    if feats:
      if len(feats) != len(weights):
        raise ValueError('feats and weights lists must be the same length')
      for feat, weight in zip(feats, weights):
        self.AddFeature(feat, weight)

  def AddFeature(self, feat, weight=None):
    if self.params and not feat.GetFamily() in self.params:
      raise ValueError('feature family %s not found in params' % feat.GetFamily())

    newFeat = FeatMapPoint()
    newFeat.initFromFeat(feat)
    newFeat.weight = weight

    self.AddFeatPoint(newFeat)

  def AddFeatPoint(self, featPt):
    if not isinstance(featPt, FeatMapPoint):
      raise ValueError('addFeatPoint() must be called with a FeatMapPoint instance')
    if self.params and not featPt.GetFamily() in self.params:
      raise ValueError('feature family %s not found in params' % featPt.GetFamily())
    self._feats.append(featPt)

  def GetFeatures(self):
    return self._feats

  def GetNumFeatures(self):
    return len(self._feats)

  def GetFeature(self, i):
    return self._feats[i]

  def DropFeature(self, i):
    del self._feats[i]

  def _loopOverMatchingFeats(self, oFeat):
    for sIdx, sFeat in enumerate(self._feats):
      if sFeat.GetFamily() == oFeat.GetFamily():
        yield sIdx, sFeat

  def GetFeatFeatScore(self, feat1, feat2, typeMatch=True):
    """ feat1 is one of our feats
        feat2 is any Feature

    """
    if typeMatch and feat1.GetFamily() != feat2.GetFamily():
      return 0.0
    d2 = feat1.GetDist2(feat2)
    params = self.params[feat1.GetFamily()]
    if d2 > params.radius * params.radius:
      return 0.0

    if params.featProfile == FeatMapParams.FeatProfile.Gaussian:
      score = math.exp(-d2 / params.width)
    elif params.featProfile == FeatMapParams.FeatProfile.Triangle:
      d = math.sqrt(d2)
      if d < params.width:
        score = 1. - d / params.width
      else:
        score = 0.0
    elif params.featProfile == FeatMapParams.FeatProfile.Box:
      score = 1.0
    score *= feat1.weight

    if self.dirScoreMode != FeatDirScoreMode.Ignore:
      dirScore = feat1.GetDirMatch(feat2)
      if self.dirScoreMode == FeatDirScoreMode.DotPosRange:
        dirScore = (dirScore + 1.0) / 2.0
      elif self.dirScoreMode != FeatDirScoreMode.DotFullRange:
        raise NotImplementedError('bad feature dir score mode')
      score *= dirScore

    return score

  def ScoreFeats(self, featsToScore, mapScoreVect=None, featsScoreVect=None,
                 featsToFeatMapIdx=None):
    nFeats = len(self._feats)
    if mapScoreVect is not None:
      if len(mapScoreVect) != nFeats:
        raise ValueError('if provided, len(mapScoreVect) should equal numFeats')
      for i in range(nFeats):
        mapScoreVect[i] = 0.0
    else:
      mapScoreVect = [0.0] * nFeats

    nToScore = len(featsToScore)
    if self.scoreMode == FeatMapScoreMode.Closest:
      defScore = 1000.0
    else:
      defScore = 0.0

    if featsScoreVect is not None:
      if len(featsScoreVect) != nToScore:
        raise ValueError('if provided, len(featsScoreVect) should equal len(featsToScore)')
      for i in range(nToScore):
        featsScoreVect[i] = defScore
    else:
      featsScoreVect = [defScore] * nToScore

    if featsToFeatMapIdx is not None:  # Initialize a 2D-empty array
      if len(featsToFeatMapIdx) != nToScore:
        raise ValueError('if provided, len(featsToFeatMapIdx) should equal len(featsToScore)')
    else:
      featsToFeatMapIdx = [None] * nToScore

    for i in range(nToScore):
      if self.scoreMode != FeatMapScoreMode.All:
        featsToFeatMapIdx[i] = [-1]
      else:
        featsToFeatMapIdx[i] = []

    for oIdx, oFeat in enumerate(featsToScore):
      for sIdx, sFeat in self._loopOverMatchingFeats(oFeat):
        if self.scoreMode == FeatMapScoreMode.Closest:
          d = sFeat.GetDist2(oFeat)
          if d < featsScoreVect[oIdx]:
            featsScoreVect[oIdx] = d
            featsToFeatMapIdx[oIdx][0] = sIdx
        else:
          lScore = self.GetFeatFeatScore(sFeat, oFeat, typeMatch=False)
          if self.scoreMode == FeatMapScoreMode.Best:
            if lScore > featsScoreVect[oIdx]:
              featsScoreVect[oIdx] = lScore
              featsToFeatMapIdx[oIdx][0] = sIdx
          elif self.scoreMode == FeatMapScoreMode.All:
            featsScoreVect[oIdx] += lScore
            mapScoreVect[sIdx] += lScore
            featsToFeatMapIdx[oIdx].append(sIdx)
          else:
            raise ValueError('bad score mode')

    totScore = 0.0
    if self.scoreMode == FeatMapScoreMode.Closest:
      for oIdx, oFeat in enumerate(featsToScore):
        sIdx = featsToFeatMapIdx[oIdx][0]
        if sIdx > -1:
          lScore = self.GetFeatFeatScore(sFeat, oFeat, typeMatch=False)
          featsScoreVect[oIdx] = lScore
          mapScoreVect[sIdx] = lScore
          totScore += lScore
        else:
          featsScoreVect[oIdx] = 0
    else:
      totScore = sum(featsScoreVect)
      if self.scoreMode == FeatMapScoreMode.Best:
        for oIdx, lScore in enumerate(featsScoreVect):
          sIdx = featsToFeatMapIdx[oIdx][0]
          if sIdx > -1:
            mapScoreVect[sIdx] = lScore

    # replace placeholders:
    if self.scoreMode != FeatMapScoreMode.All:
      for elem in featsToFeatMapIdx:
        if elem == [-1]:
          elem.pop()
    return totScore

  def __str__(self):
    res = ''
    for i, feat in enumerate(self._feats):
      weight = feat.weight
      pos = feat.GetPos()
      res += '% 3d % 12s % 6.4f % 6.4f % 6.4f % 6.4f\n' % (i + 1, feat.GetFamily(), pos.x, pos.y,
                                                           pos.z, weight)
    return res
