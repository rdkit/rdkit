#
# Copyright (C) 2001-2019  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Interface to the C++ Murtagh hierarchic clustering code

"""
import numpy

from rdkit.ML.Cluster import Clusters

try:
  from rdkit.ML.Cluster.Clustering import MurtaghCluster, MurtaghDistCluster
except ImportError:
  MurtaghCluster = None
  MurtaghDistCluster = None

# constants to select the clustering algorithm
WARDS = 1
SLINK = 2
CLINK = 3
UPGMA = 4
MCQUITTY = 5
GOWER = 6
CENTROID = 7

# descriptions of the methods:
methods = [
  ("Ward's Minimum Variance", WARDS, "Ward's Minimum Variance"),
  ('Average Linkage', UPGMA, 'Group Average Linkage (UPGMA)'),
  ('Single Linkage', SLINK, 'Single Linkage (SLINK)'),
  ('Complete Linkage', CLINK, 'Complete Linkage (CLINK)'),
  #  ("McQuitty",MCQUITTY,"McQuitty's method"),
  #  ("Gower",GOWER,"Gower's median method"),
  ("Centroid", CENTROID, "Centroid method"),
]


def _LookupDist(dists, i, j, n):
  """ *Internal Use Only*

     returns the distance between points i and j in the symmetric
     distance matrix _dists_

    """
  if i == j:
    return 0.0
  if i > j:
    i, j = j, i
  return dists[j * (j - 1) // 2 + i]


def _ToClusters(data, nPts, ia, ib, crit, isDistData=0):
  """ *Internal Use Only*

      Converts the results of the Murtagh clustering code into
      a cluster tree, which is returned in a single-entry list

    """
  cs = [None] * nPts
  for i in range(nPts):
    cs[i] = Clusters.Cluster(metric=0.0, data=i, index=(i + 1))

  nClus = len(ia) - 1
  for i in range(nClus):
    idx1 = ia[i] - 1
    idx2 = ib[i] - 1
    c1 = cs[idx1]
    c2 = cs[idx2]
    newClust = Clusters.Cluster(metric=crit[i], children=[c1, c2], index=nPts + i + 1)
    cs[idx1] = newClust

  return [newClust]


def ClusterData(data, nPts, method, isDistData=0):
  """  clusters the data points passed in and returns the cluster tree

      **Arguments**

        - data: a list of lists (or array, or whatever) with the input
          data (see discussion of _isDistData_ argument for the exception)

        - nPts: the number of points to be used

        - method: determines which clustering algorithm should be used.
            The defined constants for these are:
            'WARDS, SLINK, CLINK, UPGMA'

        - isDistData: set this toggle when the data passed in is a
            distance matrix.  The distance matrix should be stored
            symmetrically so that _LookupDist (above) can retrieve
            the results:
              for i<j: d_ij = dists[j*(j-1)//2 + i]


      **Returns**

        - a single entry list with the cluster tree
    """
  data = numpy.array(data)
  if not isDistData:
    sz = data.shape[1]
    ia, ib, crit = MurtaghCluster(data, nPts, sz, method)
  else:
    ia, ib, crit = MurtaghDistCluster(data, nPts, method)
  c = _ToClusters(data, nPts, ia, ib, crit, isDistData=isDistData)

  return c
