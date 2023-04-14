# $Id$
#
# Copyright (C) 2003-2006 Greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" class definitions for similarity screening

See _SimilarityScreener_ for overview of required API

"""
from rdkit import DataStructs
from rdkit.DataStructs import TopNContainer


class SimilarityScreener(object):
  """  base class

     important attributes:
        probe: the probe fingerprint against which we screen.

        metric: a function that takes two arguments and returns a similarity
                measure between them

        dataSource: the source pool from which to draw, needs to support
                a next() method

        fingerprinter: a function that takes a molecule and returns a
               fingerprint of the appropriate format


      **Notes**
         subclasses must support either an iterator interface
         or __len__ and __getitem__
    """

  def __init__(self, probe=None, metric=None, dataSource=None, fingerprinter=None):
    self.metric = metric
    self.dataSource = dataSource
    self.fingerprinter = fingerprinter
    self.probe = probe

  def Reset(self):
    """ used to reset screeners that behave as iterators """
    pass

  # FIX: add setters/getters for attributes
  def SetProbe(self, probeFingerprint):
    """ sets our probe fingerprint """
    self.probe = probeFingerprint

  def GetSingleFingerprint(self, probe):
    """ returns a fingerprint for a single probe object

         This is potentially useful in initializing our internal
         probe object.

        """
    return self.fingerprinter(probe)


class ThresholdScreener(SimilarityScreener):
  """ Used to return all compounds that have a similarity
      to the probe beyond a threshold value

     **Notes**:

       - This is as lazy as possible, so the data source isn't
         queried until the client asks for a hit.

       - In addition to being lazy, this class is as thin as possible.
         (Who'd have thought it was possible!)
         Hits are *not* stored locally, so if a client resets
         the iteration and starts over, the same amount of work must
         be done to retrieve the hits.

       - The thinness and laziness forces us to support only forward
         iteration (not random access)

    """

  def __init__(self, threshold, **kwargs):
    SimilarityScreener.__init__(self, **kwargs)
    self.threshold = threshold
    self.dataIter = iter(self.dataSource)

  # FIX: add setters/getters for attributes

  def _nextMatch(self):
    """ *Internal use only* """
    done = 0
    res = None
    sim = 0
    while not done:
      # this is going to crap out when the data source iterator finishes,
      #  that's how we stop when no match is found
      obj = next(self.dataIter)
      fp = self.fingerprinter(obj)
      sim = DataStructs.FingerprintSimilarity(fp, self.probe, self.metric)
      if sim >= self.threshold:
        res = obj
        done = 1
    return sim, res

  def Reset(self):
    """ used to reset our internal state so that iteration
          starts again from the beginning
        """
    self.dataSource.reset()
    self.dataIter = iter(self.dataSource)

  def __iter__(self):
    """ returns an iterator for this screener
        """
    self.Reset()
    return self

  def next(self):
    """ required part of iterator interface """
    return self._nextMatch()

  __next__ = next


class TopNScreener(SimilarityScreener):
  """ A screener that only returns the top N hits found

      **Notes**

        - supports forward iteration and getitem

    """

  def __init__(self, num, **kwargs):
    SimilarityScreener.__init__(self, **kwargs)
    self.numToGet = num
    self.topN = None
    self._pos = 0

  def Reset(self):
    self._pos = 0

  def __iter__(self):
    if self.topN is None:
      self._initTopN()
    self.Reset()
    return self

  def next(self):
    if self._pos >= self.numToGet:
      raise StopIteration
    else:
      res = self.topN[self._pos]
      self._pos += 1
      return res

  __next__ = next

  def _initTopN(self):
    self.topN = TopNContainer.TopNContainer(self.numToGet)
    for obj in self.dataSource:
      fp = self.fingerprinter(obj)
      sim = DataStructs.FingerprintSimilarity(fp, self.probe, self.metric)
      self.topN.Insert(sim, obj)

  def __len__(self):
    if self.topN is None:
      self._initTopN()
    return self.numToGet

  def __getitem__(self, idx):
    if self.topN is None:
      self._initTopN()
    return self.topN[idx]
