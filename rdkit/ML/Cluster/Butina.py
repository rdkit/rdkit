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

  # Create an empty matrix of booleans with a point for each i,j pair of molecule numbers.
  # Also tally the count of all pairs, skipping duplicates, and all self-comparison.
  # For efficiency, we never change the size of this bit bit-table.  We just flip bits.
  matrix_passing_threshold: numpy.array = numpy.zeros([nPts, nPts], dtype=bool)
  counts: numpy.array = numpy.zeros([nPts], dtype=numpy.int32)
  dmIdx = 0
  for i in range(nPts):
    for j in range(i):
      if not isDistData:
        dij = distFunc(data[i], data[j])
      else:
        dij = data[dmIdx]
        dmIdx += 1
      if dij <= distThresh:
        matrix_passing_threshold[i, j] = True
        matrix_passing_threshold[j, i] = True
        counts[i] += 1
        counts[j] += 1
      else:
        pass

  # sort by the number of neighbors:
  tLists = [(counts[i], i) for i in range(nPts)]
  tLists.sort(reverse=True)

  res = []
  seen = numpy.zeros([nPts], dtype=bool)
  while tLists:
    _, idx = tLists.pop(0)
    if seen[idx]:
      continue

    new_cluster = [idx]
    for jdx_other_cluster in range(nPts):
      if matrix_passing_threshold[idx, jdx_other_cluster]:
        if not seen[jdx_other_cluster]:
          new_cluster.append(jdx_other_cluster)
          seen[jdx_other_cluster] = True

    # update the number of neighbors:
    # remove all members of the new cluster from the list of
    # neighbors and reorder the tLists
    if reordering:
      # get the list of affected molecules, i.e. all molecules
      # which have at least one of the members of the new cluster
      # as a neighbor
      neighbors = set()
      for idx_new_cluster in new_cluster:
          for jdx_other_cluster in range(nPts):
              if matrix_passing_threshold[idx_new_cluster, jdx_other_cluster]:
                  neighbors.add(jdx_other_cluster)
      neighbors = frozenset(neighbors)

      # loop over all remaining molecules in tLists but only
      # consider unassigned and affected compounds
      for tlist_n, tlist_count_and_idx in enumerate(tLists):
        idx_remaining_cluster = tlist_count_and_idx[1]
        if seen[idx_remaining_cluster] or (idx_remaining_cluster not in neighbors):
          continue
        # update the number of neighbors
        for jdx_other_cluster in new_cluster:
            matrix_passing_threshold[idx_remaining_cluster, jdx_other_cluster] = False
        tLists[tlist_n] = (numpy.count_nonzero(matrix_passing_threshold[idx_remaining_cluster, :]), idx_remaining_cluster)

      # now reorder the list
      tLists.sort(reverse=True)
    res.append(tuple(new_cluster))
  return tuple(res)
