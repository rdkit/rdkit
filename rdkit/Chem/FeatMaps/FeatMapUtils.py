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
import copy

from rdkit.Chem.FeatMaps import FeatMaps


class MergeMethod(object):
  # Put the new point at the weighted average position of the two fused points
  WeightedAverage = 0
  # Put the new point at the un-weighted average position of the two fused points
  Average = 1
  # Put the new point at the position of the larger (by weight) of the two points
  UseLarger = 2

  @classmethod
  def valid(cls, mergeMethod):
    """ Check that mergeMethod is valid """
    if mergeMethod not in (cls.WeightedAverage, cls.Average, cls.UseLarger):
      raise ValueError('unrecognized mergeMethod')


class MergeMetric(object):
  # Do not merge points
  NoMerge = 0
  # merge two points if they come within a threshold distance
  Distance = 1
  # merge two points if their percent overlap exceeds a threshold
  Overlap = 2

  @classmethod
  def valid(cls, mergeMetric):
    """ Check that mergeMetric is valid """
    if mergeMetric not in (cls.NoMerge, cls.Distance, cls.Overlap):
      raise ValueError('unrecognized mergeMetric')


class DirMergeMode(object):
  # Do not merge directions (i.e. keep all direction vectors)
  NoMerge = 0
  # Sum direction vectors
  Sum = 1

  @classmethod
  def valid(cls, dirMergeMode):
    """ Check that dirMergeMode is valid """
    if dirMergeMode not in (cls.NoMerge, cls.Sum):
      raise ValueError('unrecognized dirMergeMode')


def __copyAll(res, fm1, fm2):
  """ no user-serviceable parts inside """
  for feat in fm1.GetFeatures():
    res.AddFeatPoint(copy.deepcopy(feat))
  for feat in fm2.GetFeatures():
    res.AddFeatPoint(copy.deepcopy(feat))


def GetFeatFeatDistMatrix(fm, mergeMetric, mergeTol, dirMergeMode, compatFunc):
  """

    NOTE that mergeTol is a max value for merging when using distance-based
    merging and a min value when using score-based merging.

  """
  MergeMetric.valid(mergeMetric)
  numFeatures = fm.GetNumFeatures()
  dists = [[1e8] * numFeatures for _ in range(numFeatures)]
  if mergeMetric == MergeMetric.NoMerge:
    return dists

  # Setup distance matrix, depending on mergeMetric.
  benchmarkDict = {MergeMetric.Distance: mergeTol * mergeTol, MergeMetric.Overlap: mergeTol}
  benchmark = benchmarkDict[mergeMetric]

  def assignMatrix(matrix, i, j, value, constraint):
    if value < constraint:
      matrix[i][j] = value
      matrix[j][i] = value

  getFeature = fm.GetFeature
  for i in range(numFeatures):
    ptI = getFeature(i)
    for j in range(i + 1, numFeatures):
      ptJ = getFeature(j)
      if compatFunc(ptI, ptJ):
        if mergeMetric == MergeMetric.Distance:
          dist2 = ptI.GetDist2(ptJ)
          assignMatrix(matrix=dists, i=i, j=j, value=dist2, constraint=benchmark)
        elif mergeMetric == MergeMetric.Overlap:
          score = fm.GetFeatFeatScore(ptI, ptJ, typeMatch=False) * (-1 * ptJ.weight)
          assignMatrix(matrix=dists, i=i, j=j, value=score, constraint=benchmark)

  return dists


def familiesMatch(f1, f2):
  return f1.GetFamily() == f2.GetFamily()


def feq(v1, v2, tol=1e-4):
  return abs(v1 - v2) < tol


def MergeFeatPoints(fm, mergeMetric=MergeMetric.NoMerge, mergeTol=1.5,
                    dirMergeMode=DirMergeMode.NoMerge, mergeMethod=MergeMethod.WeightedAverage,
                    compatFunc=familiesMatch):
  """

    NOTE that mergeTol is a max value for merging when using distance-based
    merging and a min value when using score-based merging.

    returns whether or not any points were actually merged

  """
  MergeMetric.valid(mergeMetric)
  MergeMethod.valid(mergeMethod)
  DirMergeMode.valid(dirMergeMode)

  res = False
  if mergeMetric == MergeMetric.NoMerge:
    return res
  dists = GetFeatFeatDistMatrix(fm, mergeMetric, mergeTol, dirMergeMode, compatFunc)
  distOrders = [None] * len(dists)
  for i, distV in enumerate(dists):
    distOrders[i] = []
    for j, dist in enumerate(distV):
      if dist < mergeTol:
        distOrders[i].append((dist, j))
    distOrders[i].sort()

  # print('distOrders:')
  # print(distOrders)

  # we now know the "distances" and have rank-ordered list of
  # each point's neighbors. Work with that.

  # progressively merge nearest neighbors until there
  # are no more points left to merge
  featsInPlay = list(range(fm.GetNumFeatures()))
  featsToRemove = []
  # print '--------------------------------'
  while featsInPlay:
    # find two features who are mutual nearest neighbors:
    fipCopy = featsInPlay[:]
    for fi in fipCopy:
      # print('>>>',fi,fipCopy,featsInPlay)
      # print('\t',distOrders[fi])
      mergeThem = False
      if not distOrders[fi]:
        featsInPlay.remove(fi)
        continue
      dist, nbr = distOrders[fi][0]
      if nbr not in featsInPlay:
        continue
      if distOrders[nbr][0][1] == fi:
        # print 'direct:',fi,nbr
        mergeThem = True
      else:
        # it may be that there are several points at about the same distance,
        # check for that now
        if feq(distOrders[nbr][0][0], dist):
          for distJ, nbrJ in distOrders[nbr][1:]:
            if feq(dist, distJ):
              if nbrJ == fi:
                # print 'indirect: ',fi,nbr
                mergeThem = True
                break
            else:
              break
      # print '    bottom:',mergeThem
      if mergeThem:
        break

    if mergeThem:
      res = True
      featI = fm.GetFeature(fi)
      nbrFeat = fm.GetFeature(nbr)

      if mergeMethod == MergeMethod.WeightedAverage:
        newPos = featI.GetPos() * featI.weight + nbrFeat.GetPos() * nbrFeat.weight
        newPos /= (featI.weight + nbrFeat.weight)
        newWeight = (featI.weight + nbrFeat.weight) / 2
      elif mergeMethod == MergeMethod.Average:
        newPos = featI.GetPos() + nbrFeat.GetPos()
        newPos /= 2
        newWeight = (featI.weight + nbrFeat.weight) / 2
      elif mergeMethod == MergeMethod.UseLarger:
        if featI.weight > nbrFeat.weight:
          newPos = featI.GetPos()
          newWeight = featI.weight
        else:
          newPos = nbrFeat.GetPos()
          newWeight = nbrFeat.weight

      featI.SetPos(newPos)
      featI.weight = newWeight

      # nbr and fi are no longer valid targets:
      # print 'nbr done:',nbr,featsToRemove,featsInPlay
      featsToRemove.append(nbr)
      featsInPlay.remove(fi)
      featsInPlay.remove(nbr)
      for nbrList in distOrders:
        try:
          nbrList.remove(fi)
        except ValueError:
          pass
        try:
          nbrList.remove(nbr)
        except ValueError:
          pass
    else:
      # print ">>>> Nothing found, abort"
      break
  featsToRemove.sort()
  for i, fIdx in enumerate(featsToRemove):
    fm.DropFeature(fIdx - i)
  return res


def CombineFeatMaps(fm1, fm2, mergeMetric=MergeMetric.NoMerge, mergeTol=1.5,
                    dirMergeMode=DirMergeMode.NoMerge):
  """
     the parameters will be taken from fm1
  """
  res = FeatMaps.FeatMap(params=fm1.params)

  __copyAll(res, fm1, fm2)
  if mergeMetric != MergeMetric.NoMerge:
    MergeFeatPoints(res, mergeMetric=mergeMetric, mergeTol=mergeTol)
  return res
