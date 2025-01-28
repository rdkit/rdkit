#
# Copyright (C) 2004-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#


class ExcludedVolume(object):

  def __init__(self, featInfo, index=-1, exclusionDist=3.0):
    """
    featInfo should be a sequence of ([indices],min,max) tuples

    """
    self.index = index
    try:
      _ = len(featInfo)
    except TypeError:
      raise ValueError('featInfo argument must be a sequence of sequences')

    if not len(featInfo):
      raise ValueError('featInfo argument must non-empty')

    try:
      _, _, _ = featInfo[0]
    except (TypeError, ValueError):
      raise ValueError('featInfo elements must be 3-sequences')

    self.featInfo = featInfo[:]
    self.exclusionDist = exclusionDist
    self.pos = None
