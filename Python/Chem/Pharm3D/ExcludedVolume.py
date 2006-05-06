# $Id: ExcludedVolume.py 5022 2006-03-02 01:34:20Z glandrum $
#
# Copyright (C) 2004-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
class ExcludedVolume(object):
  def __init__(self, featInfo,index=-1,exclusionDist=3.0):
    """
    featInfo should be a sequence of ([indices],min,max) tuples

    """
    self.index = index
    try:
      l = len(featInfo)
    except AttributeError:
      raise ValueError,'featInfo argument must be a sequence of sequences'

    if not len(featInfo):
      raise ValueError,'featInfo argument must non-empty'
      
    try:
      a,b,c = featInfo[0]
    except Type:
      raise ValueError,'featInfo elements must be 3-sequences'
    except ValueError:
      raise ValueError,'featInfo elements must be 3-sequences'

    self.featInfo = featInfo[:]
    self.exclusionDist = exclusionDist
    self.pos = None
    
