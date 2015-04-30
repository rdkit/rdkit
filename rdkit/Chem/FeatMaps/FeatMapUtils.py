## Automatically adapted for numpy.oldnumeric Jun 27, 2008 by -c

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
from rdkit.six.moves import range
from rdkit.Chem.FeatMaps import FeatMaps

class MergeMethod(object):
  WeightedAverage=0
  """ Put the new point at the weighted average position of the two
   fused points
  """

  Average=1
  """ Put the new point at the un-weighted average position of the two
   fused points
  """

  UseLarger=2
  """ Put the new point at the position of the larger (by weight)
   of the two points
  """

class MergeMetric(object):
  NoMerge=0
  """ Do not merge points """

  Distance=1
  """ merge two points if they come within a threshold distance """

  Overlap=2
  """ merge two points if their percent overlap exceeds a threshold """

class DirMergeMode(object):
  NoMerge=0
  """ Do not merge directions (i.e. keep all direction vectors) """

  Sum=1
  """ Sum direction vectors """

def __copyAll(res,fm1,fm2):
  """ no user-serviceable parts inside """
  for feat in fm1.GetFeatures():
    res.AddFeatPoint(copy.deepcopy(feat))
  for feat in fm2.GetFeatures():
    res.AddFeatPoint(copy.deepcopy(feat))


def GetFeatFeatDistMatrix(fm,mergeMetric,mergeTol,dirMergeMode,compatFunc):
  """

    NOTE that mergeTol is a max value for merging when using distance-based
    merging and a min value when using score-based merging.
    
  """
  dists = [[1e8]*fm.GetNumFeatures() for x in range(fm.GetNumFeatures())]
  if mergeMetric==MergeMetric.NoMerge:
    return dists
  elif mergeMetric==MergeMetric.Distance:
    mergeTol2 = mergeTol*mergeTol
    for i in range(fm.GetNumFeatures()):
      ptI = fm.GetFeature(i)
      for j in range(i+1,fm.GetNumFeatures()):
        ptJ = fm.GetFeature(j)
        if compatFunc(ptI,ptJ):
          dist2 = ptI.GetDist2(ptJ)
          if dist2<mergeTol2:
            dists[i][j]=dist2
            dists[j][i]=dist2
  elif mergeMetric==MergeMetric.Overlap:
    for i in range(fm.GetNumFeatures()):
      ptI = fm.GetFeature(i)
      for j in range(i+1,fm.GetNumFeatures()):
        ptJ = fm.GetFeature(j)
        if compatFunc(ptI,ptJ):
          score = fm.GetFeatFeatScore(ptI,ptJ,typeMatch=False)
          score *= -1*ptJ.weight
          if score<mergeTol:
            dists[i][j]=score
            dists[j][i]=score
  else:
    raise ValueError('unrecognized mergeMetric')

  return dists

def familiesMatch(f1,f2):
  return f1.GetFamily()==f2.GetFamily()

def feq(v1,v2,tol=1e-4):
  return abs(v1-v2)<tol

def MergeFeatPoints(fm,mergeMetric=MergeMetric.NoMerge,mergeTol=1.5,
                    dirMergeMode=DirMergeMode.NoMerge,
                    mergeMethod=MergeMethod.WeightedAverage,
                    compatFunc=familiesMatch):
  """

    NOTE that mergeTol is a max value for merging when using distance-based
    merging and a min value when using score-based merging.

    returns whether or not any points were actually merged
    
  """
  res=False
  if mergeMetric==MergeMetric.NoMerge:
    return res
  dists = GetFeatFeatDistMatrix(fm,mergeMetric,mergeTol,dirMergeMode,compatFunc)
  distOrders = [None]*len(dists)
  for i in range(len(dists)):
    distV = dists[i]
    distOrders[i] = []
    for j,dist in enumerate(distV):
      if dist<mergeTol:
        distOrders[i].append((dist,j))
    distOrders[i].sort()

  #print 'distOrders:'
  #print distOrders

  # we now know the "distances" and have rank-ordered list of
  # each point's neighbors. Work with that.

  # progressively merge nearest neighbors until there
  # are no more points left to merge
  featsInPlay=list(range(fm.GetNumFeatures()))
  featsToRemove = []
  #print '--------------------------------'
  while featsInPlay:
    # find two features who are mutual nearest neighbors:
    fipCopy=featsInPlay[:]
    for fi in fipCopy:
      #print '>>>',fi,fipCopy,featsInPlay
      #print '\t',distOrders[fi]
      mergeThem=False
      if not distOrders[fi]:
        featsInPlay.remove(fi)
        continue
      dist,nbr = distOrders[fi][0]
      if nbr not in featsInPlay:
        continue
      if distOrders[nbr][0][1]==fi:
        #print 'direct:',fi,nbr
        mergeThem=True
      else:
        # it may be that there are several points at about the same distance,
        # check for that now
        if(feq(distOrders[nbr][0][0],dist)):
          for distJ,nbrJ in distOrders[nbr][1:]:
            if feq(dist,distJ):
              if nbrJ==fi:
                #print 'indirect: ',fi,nbr
                mergeThem=True
                break
            else:
              break
      #print '    bottom:',mergeThem
      if mergeThem: break
    if mergeThem:
      res=True
      featI = fm.GetFeature(fi)
      nbrFeat = fm.GetFeature(nbr)
      
      if mergeMethod==MergeMethod.WeightedAverage:
        newPos = featI.GetPos()*featI.weight+nbrFeat.GetPos()*nbrFeat.weight
        newPos /= (featI.weight+nbrFeat.weight)
        newWeight = (featI.weight+nbrFeat.weight)/2
      elif mergeMethod==MergeMethod.Average:
        newPos = featI.GetPos()+nbrFeat.GetPos()
        newPos /= 2
        newWeight = (featI.weight+nbrFeat.weight)/2
      elif mergeMethod==MergeMethod.UseLarger:
        if featI.weight>nbrFeat.weight:
          newPos=featI.GetPos()
          newWeight = featI.weight
        else:
          newPos=nbrFeat.GetPos()
          newWeight = nbrFeat.weight
      else:
        raise ValueError("bad mergeMethod")

      featI.SetPos(newPos)
      featI.weight = newWeight
      
      # nbr and fi are no longer valid targets:
      #print 'nbr done:',nbr,featsToRemove,featsInPlay
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
      #print ">>>> Nothing found, abort"
      break
  featsToRemove.sort()
  for i,fIdx in enumerate(featsToRemove):
    fm.DropFeature(fIdx-i)
  return res

def CombineFeatMaps(fm1,fm2,mergeMetric=MergeMetric.NoMerge,mergeTol=1.5,
                    dirMergeMode=DirMergeMode.NoMerge):
  """
     the parameters will be taken from fm1
  """
  res = FeatMaps.FeatMap(params=fm1.params)

  __copyAll(res,fm1,fm2)
  if mergeMetric!=MergeMetric.NoMerge:
    MergeFeatPoints(res,mergeMetric=mergeMetric,mergeTol=mergeTol)
  return res

#------------------------------------
#
#  doctest boilerplate
#
def _test():
  import doctest,sys
  return doctest.testmod(sys.modules["__main__"])

if __name__ == '__main__':
  import sys
  failed,tried = _test()
  sys.exit(failed)
