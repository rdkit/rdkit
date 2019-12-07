# $Id$
#
# Copyright (C) 2001-2008  greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""utility functions for clustering

"""


def GetNodeList(cluster):
  """returns an ordered list of all nodes below cluster

  the ordering is done using the lengths of the child nodes

   **Arguments**

     - cluster: the cluster in question

   **Returns**

     - a list of the leaves below this cluster

  """
  if len(cluster) == 1:
    return [cluster]
  else:
    children = cluster.GetChildren()
    children.sort(key=lambda x: len(x), reverse=True)
    res = []
    for child in children:
      res += GetNodeList(child)
    res += [cluster]
    return res


def GetNodesDownToCentroids(cluster, above=1):
  """returns an ordered list of all nodes below cluster


  """
  if hasattr(cluster, '_isCentroid'):
    cluster._aboveCentroid = 0
    above = -1
  else:
    cluster._aboveCentroid = above
  if len(cluster) == 1:
    return [cluster]
  else:
    res = []
    children = cluster.GetChildren()
    children.sort(key=lambda x: len(x), reverse=True)
    #     children.sort(lambda x, y: cmp(len(y), len(x)))
    for child in children:
      res = res + GetNodesDownToCentroids(child, above)
    res = res + [cluster]
    return res


def FindClusterCentroidFromDists(cluster, dists):
  """ find the point in a cluster which has the smallest summed
     Euclidean distance to all others

   **Arguments**

     - cluster: the cluster to work with

     - dists: the distance matrix to use for the points

   **Returns**

     - the index of the centroid point

  """
  children = cluster.GetPoints()
  pts = [x.GetData() for x in children]

  best = 1e24
  bestIdx = -1
  for pt in pts:
    dAccum = 0.0
    # loop over others and add'em up
    for other in pts:
      if other != pt:
        if other > pt:
          row, col = pt, other
        else:
          row, col = other, pt
        dAccum += dists[col * (col - 1) / 2 + row]
        if dAccum >= best:
          # minor efficiency hack
          break
    if dAccum < best:
      best = dAccum
      bestIdx = pt
  for i in range(len(pts)):
    pt = pts[i]
    if pt != bestIdx:
      if pt > bestIdx:
        row, col = bestIdx, pt
      else:
        row, col = pt, bestIdx
      children[i]._distToCenter = dists[col * (col - 1) / 2 + row]
    else:
      children[i]._distToCenter = 0.0
    children[i]._clustCenter = bestIdx
  cluster._clustCenter = bestIdx
  cluster._distToCenter = 0.0

  return bestIdx


def _BreadthFirstSplit(cluster, n):
  """  *Internal Use Only*

  """
  if len(cluster) < n:
    raise ValueError('Cannot split cluster of length %d into %d pieces' % (len(cluster), n))
  if len(cluster) == n:
    return cluster.GetPoints()
  clusters = [cluster]
  nxtIdx = 0
  for _ in range(n - 1):
    while nxtIdx < len(clusters) and len(clusters[nxtIdx]) == 1:
      nxtIdx += 1
    assert nxtIdx < len(clusters)

    children = clusters[nxtIdx].GetChildren()
    children.sort(key=lambda x: x.GetMetric(), reverse=True)
    for child in children:
      clusters.append(child)
    del clusters[nxtIdx]
  return clusters


def _HeightFirstSplit(cluster, n):
  """  *Internal Use Only*

  """
  if len(cluster) < n:
    raise ValueError('Cannot split cluster of length %d into %d pieces' % (len(cluster), n))
  if len(cluster) == n:
    return cluster.GetPoints()
  clusters = [cluster]
  for _ in range(n - 1):
    nxtIdx = 0
    while nxtIdx < len(clusters) and len(clusters[nxtIdx]) == 1:
      nxtIdx += 1
    assert nxtIdx < len(clusters)

    children = clusters[nxtIdx].GetChildren()
    for child in children:
      clusters.append(child)
    del clusters[nxtIdx]
    clusters.sort(key=lambda x: x.GetMetric(), reverse=True)
  return clusters


def SplitIntoNClusters(cluster, n, breadthFirst=True):
  """  splits a cluster tree into a set of branches

    **Arguments**

      - cluster: the root of the cluster tree

      - n: the number of clusters to include in the split

      - breadthFirst: toggles breadth first (vs depth first) cleavage
        of the cluster tree.

    **Returns**

      - a list of sub clusters

  """
  if breadthFirst:
    return _BreadthFirstSplit(cluster, n)
  else:
    return _HeightFirstSplit(cluster, n)
