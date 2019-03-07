# $Id$
#
#  Copyright (C) 2003-2013 Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

import bisect


class TopNContainer(object):
  """ maintains a sorted list of a particular number of data elements.

  """

  def __init__(self, size, mostNeg=-1e99):
    """
    if size is negative, all entries will be kept in sorted order
    """
    self._size = size
    if (size >= 0):
      self.best = [mostNeg] * self._size
      self.extras = [None] * self._size
    else:
      self.best = []
      self.extras = []

  def Insert(self, val, extra=None):
    """ only does the insertion if val fits """
    if self._size >= 0:
      if val > self.best[0]:
        idx = bisect.bisect(self.best, val)
        # insert the new element
        if idx == self._size:
          self.best.append(val)
          self.extras.append(extra)
        else:
          self.best.insert(idx, val)
          self.extras.insert(idx, extra)
        # and pop off the head
        self.best.pop(0)
        self.extras.pop(0)
    else:
      idx = bisect.bisect(self.best, val)
      self.best.insert(idx, val)
      self.extras.insert(idx, extra)

  def GetPts(self):
    """ returns our set of points """
    return self.best

  def GetExtras(self):
    """ returns our set of extras """
    return self.extras

  def __len__(self):
    if self._size >= 0:
      return self._size
    else:
      return len(self.best)

  def __getitem__(self, which):
    return self.best[which], self.extras[which]

  def reverse(self):
    self.best.reverse()
    self.extras.reverse()


def _exampleCode():
  import random
  random.seed(0)
  pts = [int(100 * random.random()) for _ in range(10)]

  c = TopNContainer(4)
  for pt in pts:
    c.Insert(pt, extra=str(pt))
  print(c.GetPts())
  print(c.GetExtras())


if __name__ == '__main__':  # pragma: nocover
  _exampleCode()
