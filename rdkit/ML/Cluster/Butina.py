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
logger = RDLogger.logger()


def EuclideanDist(pi, pj):
  dv = numpy.array(pi) - numpy.array(pj)
  return numpy.sqrt(dv * dv)


def ClusterData(data, nPts, distThresh, isDistData=False, distFunc=EuclideanDist, reordering=False):
  """  clusters the data points passed in and returns the list of clusters

    **Arguments**

      - data: a list of items with the input data
        (see discussion of _isDistData_ argument for the exception)

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

      - distFunc: a function to calculate distances between points.
           Receives 2 points as arguments, should return a float

      - reordering: if this toggle is set, the number of neighbors is updated
           for the unassigned molecules after a new cluster is created such
           that always the molecule with the largest number of unassigned
           neighbors is selected as the next cluster center.

    **Returns**

      - a tuple of tuples containing information about the clusters:
         ( (cluster1_elem1, cluster1_elem2, ...),
           (cluster2_elem1, cluster2_elem2, ...),
           ...
         )
         The first element for each cluster is its centroid.

  """
  if isDistData and len(data) > (nPts * (nPts - 1) / 2):
    logger.warning("Distance matrix is too long")
  nbrLists = [None] * nPts
  for i in range(nPts):
    nbrLists[i] = []

  dmIdx = 0
  for i in range(nPts):
    for j in range(i):
      if not isDistData:
        dij = distFunc(data[i], data[j])
      else:
        dij = data[dmIdx]
        dmIdx += 1
      if dij <= distThresh:
        nbrLists[i].append(j)
        nbrLists[j].append(i)
  # sort by the number of neighbors:
  tLists = [(len(y), x) for x, y in enumerate(nbrLists)]
  tLists.sort(reverse=True)

  res = []
  seen = [0] * nPts
  while tLists:
    _, idx = tLists.pop(0)
    if seen[idx]:
      continue
    tRes = [idx]
    for nbr in nbrLists[idx]:
      if not seen[nbr]:
        tRes.append(nbr)
        seen[nbr] = 1
    # update the number of neighbors:
    # remove all members of the new cluster from the list of
    # neighbors and reorder the tLists
    if reordering:
      # get the list of affected molecules, i.e. all molecules
      # which have at least one of the members of the new cluster
      # as a neighbor
      nbrNbr = [nbrLists[t] for t in tRes]
      nbrNbr = frozenset().union(*nbrNbr)
      # loop over all remaining molecules in tLists but only
      # consider unassigned and affected compounds
      for x, y in enumerate(tLists):
        y1 = y[1]
        if seen[y1] or (y1 not in nbrNbr):
          continue
        # update the number of neighbors
        nbrLists[y1] = set(nbrLists[y1]).difference(tRes)
        tLists[x] = (len(nbrLists[y1]), y1)
      # now reorder the list
      tLists.sort(reverse=True)
    res.append(tuple(tRes))
  return tuple(res)
