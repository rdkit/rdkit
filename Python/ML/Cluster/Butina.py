# $Id$
#
# Copyright (C) 2007-2008 Greg Landrum
#   All Rights Reserved
#
""" Implementation of the clustering algorithm published in:
  Butina JCICS 39 747-750 (1999)

"""
import numpy
from rdkit import RDLogger
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
          symmetrically. An example of how to do this:

            dists = []
            for i in range(nPts):
              for j in range(i):
                dists.append( distfunc(i,j) )
          
    **Returns**

      - a tuple of tuples containing information about the clusters:
         ( (cluster1_elem1, cluster1_elem2, ...),
           (cluster2_elem1, cluster2_elem2, ...),
           ...
         )  
         The first element for each cluster is its centroid.

  """
  data = numpy.array(data)
  if isDistData and len(data)>(nPts*(nPts-1)/2):
    logger.warning("Distance matrix is too long")
  nbrLists = [None]*nPts
  for i in range(nPts): nbrLists[i] = []

  dmIdx=0
  for i in range(nPts):
    for j in range(i):
      if not isDistData:
        dv = data[i]-data[j]
        dij = numpy.sqrt(dv*dv)
      else:
        dij = data[dmIdx]
        dmIdx+=1
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
    
