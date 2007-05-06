# $Id$
#
# Copyright (C) 2007 Greg Landrum
#   All Rights Reserved
#
""" Implementation of the clustering algorithm published in:
  Butina JCICS 39 747-750 (1999)

"""
from Numeric import *
import RDLogger
logger=RDLogger.logger()

def ClusterData(data,nPts,distThresh,isDistData=False):
  """  clusters the data points passed in and returns the list of clusters

    **Arguments**

      - data: a list of lists (or array, or whatever) with the input
        data (see discussion of _isDistData_ argument for the exception)

      - nPts: the number of points to be used

      - distThresh: elements within this range of each other are considered
        to be neighbors            

      - isDistData: set this toggle when the data passed in is a
          distance matrix.  The distance matrix should be stored
          symmetrically so that _LookupDist (above) can retrieve
          the results:
            for i<j: d_ij = dists[j*(j-1)/2 + i]

    **Returns**

      - a tuple of tuples containing information about the clusters:
         ( (cluster1_elem1, cluster1_elem2, ...),
           (cluster2_elem1, cluster2_elem2, ...),
           ...
         )  
         The first element for each cluster is its centroid.

  """
  data = array(data)
  if isDistData and len(data)>(nPts*(nPts-1)/2):
    logger.warning("Distance matrix is too long")
  nbrLists = [None]*nPts
  for i in range(nPts): nbrLists[i] = []

  for i in range(nPts):
    for j in range(i+1,nPts):
      if not isDistData:
        dv = data[i]-data[j]
        dij = sqrt(dv*dv)
      else:
        dij = data[(j*(j-1))/2+i]
        #print i,j,dij
      if dij<=distThresh:
        nbrLists[i].append(j)
        nbrLists[j].append(i)
  #print nbrLists
  # sort by the number of neighbors:
  tLists = [(len(y),x) for x,y in enumerate(nbrLists)]
  tLists.sort()
  tLists.reverse()

  res = []
  seen = [0]*nPts
  while tLists:
    nNbrs,idx = tLists.pop(0)
    if seen[idx]:
      continue
    tRes = [idx]
    for nbr in nbrLists[idx]:
      if not seen[nbr]:
        tRes.append(nbr)
        seen[nbr]=1
    res.append(tuple(tRes))
  return tuple(res)
    
