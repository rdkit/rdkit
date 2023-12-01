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
""" contains the Cluster class for representing hierarchical cluster trees

"""


def cmp(t1, t2):
  return (t1 < t2) * -1 or (t1 > t2) * 1


CMPTOL = 1e-6


class Cluster(object):
  """a class for storing clusters/data

      **General Remarks**

       - It is assumed that the bottom of any cluster hierarchy tree is composed of
         the individual data points which were clustered.

       - Clusters objects store the following pieces of data, most are
          accessible via standard Setters/Getters:

         - Children: *Not Settable*, the list of children.  You can add children
           with the _AddChild()_ and _AddChildren()_ methods.

           **Note** this can be of arbitrary length,
           but the current algorithms I have only produce trees with two children
           per cluster

         - Metric: the metric for this cluster (i.e. how far apart its children are)

         - Index: the order in which this cluster was generated

         - Points: *Not Settable*, the list of original points in this cluster
              (calculated recursively from the children)

         - PointsPositions: *Not Settable*, the list of positions of the original
            points in this cluster (calculated recursively from the children)

         - Position: the location of the cluster **Note** for a cluster this
           probably means the location of the average of all the Points which are
           its children.

         - Data: a data field.  This is used with the original points to store their
           data value (i.e. the value we're using to classify)

         - Name: the name of this cluster

    """

  def __init__(self, metric=0.0, children=None, position=None, index=-1, name=None, data=None):
    """Constructor

          **Arguments**

            see the class documentation for the meanings of these arguments

            *my wrists are tired*

        """
    if children is None:
      children = []
    if position is None:
      position = []
    self.metric = metric
    self.children = children
    self._UpdateLength()
    self.pos = position
    self.index = index
    self.name = name
    self._points = None
    self._pointsPositions = None
    self.data = data

  def SetMetric(self, metric):
    self.metric = metric

  def GetMetric(self):
    return self.metric

  def SetIndex(self, index):
    self.index = index

  def GetIndex(self):
    return self.index

  def SetPosition(self, pos):
    self.pos = pos

  def GetPosition(self):
    return self.pos

  def GetPointsPositions(self):
    if self._pointsPositions is not None:
      return self._pointsPositions
    else:
      self._GenPoints()
      return self._pointsPositions

  def GetPoints(self):
    if self._points is not None:
      return self._points
    else:
      self._GenPoints()
      return self._points

  def FindSubtree(self, index):
    """ finds and returns the subtree with a particular index
        """
    res = None
    if index == self.index:
      res = self
    else:
      for child in self.children:
        res = child.FindSubtree(index)
        if res:
          break
    return res

  def _GenPoints(self):
    """ Generates the _Points_ and _PointsPositions_ lists

         *intended for internal use*

        """
    if len(self) == 1:
      self._points = [self]
      self._pointsPositions = [self.GetPosition()]
      return self._points
    else:
      res = []
      children = self.GetChildren()
      children.sort(key=lambda x: len(x), reverse=True)
      for child in children:
        res += child.GetPoints()
      self._points = res
      self._pointsPositions = [x.GetPosition() for x in res]

  def AddChild(self, child):
    """Adds a child to our list

          **Arguments**

            - child: a Cluster

        """
    self.children.append(child)
    self._GenPoints()
    self._UpdateLength()

  def AddChildren(self, children):
    """Adds a bunch of children to our list

          **Arguments**

            - children: a list of Clusters

        """
    self.children += children
    self._GenPoints()
    self._UpdateLength()

  def RemoveChild(self, child):
    """Removes a child from our list

          **Arguments**

            - child: a Cluster

        """
    self.children.remove(child)
    self._UpdateLength()

  def GetChildren(self):
    self.children.sort(key=lambda x: x.GetMetric())
    return self.children

  def SetData(self, data):
    self.data = data

  def GetData(self):
    return self.data

  def SetName(self, name):
    self.name = name

  def GetName(self):
    if self.name is None:
      return 'Cluster(%d)' % (self.GetIndex())
    else:
      return self.name

  def Print(self, level=0, showData=0, offset='\t'):
    if not showData or self.GetData() is None:
      print('%s%s%s Metric: %f' % ('  ' * level, self.GetName(), offset, self.GetMetric()))
    else:
      print('%s%s%s Data: %f\t Metric: %f' %
            ('  ' * level, self.GetName(), offset, self.GetData(), self.GetMetric()))

    for child in self.GetChildren():
      child.Print(level=level + 1, showData=showData, offset=offset)

  def Compare(self, other, ignoreExtras=1):
    """ not as choosy as self==other

        """
    tv1, tv2 = str(type(self)), str(type(other))
    tv = cmp(tv1, tv2)
    if tv:
      return tv
    tv1, tv2 = len(self), len(other)
    tv = cmp(tv1, tv2)
    if tv:
      return tv

    if not ignoreExtras:
      m1, m2 = self.GetMetric(), other.GetMetric()
      if abs(m1 - m2) > CMPTOL:
        return cmp(m1, m2)

      if cmp(self.GetName(), other.GetName()):
        return cmp(self.GetName(), other.GetName())

      sP = self.GetPosition()
      oP = other.GetPosition()
      try:
        r = cmp(len(sP), len(oP))
      except Exception:
        pass
      else:
        if r:
          return r

      try:
        r = cmp(sP, oP)
      except Exception:
        r = sum(sP - oP)
      if r:
        return r

    c1, c2 = self.GetChildren(), other.GetChildren()
    if cmp(len(c1), len(c2)):
      return cmp(len(c1), len(c2))
    for i in range(len(c1)):
      t = c1[i].Compare(c2[i], ignoreExtras=ignoreExtras)
      if t:
        return t

    return 0

  def _UpdateLength(self):
    """ updates our length

         *intended for internal use*

        """
    self._len = sum(len(c) for c in self.children) + 1

  def IsTerminal(self):
    return self._len <= 1

  def __len__(self):
    """ allows _len(cluster)_ to work

        """
    return self._len

  def __cmp__(self, other):
    """ allows _cluster1 == cluster2_ to work

        """
    if cmp(type(self), type(other)):
      return cmp(type(self), type(other))

    m1, m2 = self.GetMetric(), other.GetMetric()
    if abs(m1 - m2) > CMPTOL:
      return cmp(m1, m2)

    if cmp(self.GetName(), other.GetName()):
      return cmp(self.GetName(), other.GetName())

    c1, c2 = self.GetChildren(), other.GetChildren()
    return cmp(c1, c2)


if __name__ == '__main__':  # pragma: nocover
  from rdkit.ML.Cluster import ClusterUtils
  root = Cluster(index=1, metric=1000)
  c1 = Cluster(index=10, metric=100)
  c1.AddChild(Cluster(index=30, metric=10))
  c1.AddChild(Cluster(index=31, metric=10))
  c1.AddChild(Cluster(index=32, metric=10))

  c2 = Cluster(index=11, metric=100)
  c2.AddChild(Cluster(index=40, metric=10))
  c2.AddChild(Cluster(index=41, metric=10))

  root.AddChild(c1)
  root.AddChild(c2)

  nodes = ClusterUtils.GetNodeList(root)

  indices = [x.GetIndex() for x in nodes]
  print('XXX:', indices)
