# $Id$
#
# Copyright (C) 2001-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" defines the ClusterModel

"""
from Numeric import *
from ML.Cluster import Classify

class ClusterModel(object):
  """ this is a model suitable for inclusion in a _ML.Composite.Composite_

  """
  def ClassifyExample(self,example):
    """ returns a classification for the example passed in

    """
    d = [array(example[1:-1])]
    res = Classify.FindClosestClusters(self.clusters,d,clusterPos=self.clusterPos)
    return self.clusters[res[0]]._classification
  
  def __init__(self,clusters=[]):
    """ Constructor

      **Arguments**

        - clusters: the list of clusters to be used here

      **Notes**

        - _clusters_ is not copied, so beware.

    """
    self.clusters = clusters
    self.clusterPos = array(map(lambda x:x.GetPosition(),self.clusters))
