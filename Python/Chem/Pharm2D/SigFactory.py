# $Id$
#
# Copyright (C) 2003-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" contains factory class for producing signatures


"""
from Chem.Pharm2D.Signature import Pharm2DSig as Signature

class SigFactory(object):
  """

    SigFactory's are used by creating one, setting the relevant
    parameters, then calling the GetSignature() method each time a
    signature is required.

  """
  def __init__(self):
    self._patts = None
    self._lables = None
    self._bins = None
    self._minCnt = 2
    self._maxCnt = 4
    self._shortestPathsOnly = 1
    self._includeBondOrder = 0
    self._labels = None
    
  def SetPatterns(self,patts):
    self._patts = patts[:]
  def SetPatternsFromSmarts(self,smarts):
    import Chem
    self._patts = [None]*len(smarts)
    for i in range(len(smarts)):
      p = Chem.MolFromSmarts(smarts[i])
      self._patts[i] = p
  def GetPatterns(self):
    return self._patts
  def GetNumPatterns(self):
    return len(self._patts)
  
  def SetLabels(self,labels):
    self._labels = labels[:]
  def GetLabel(self,which):
    return self._labels[which]
  def GetLabels(self):
    return self._labels

  def SetBins(self,bins):
    """ bins should be a list of 2-tuples """
    self._bins = bins[:]
  def GetBins(self):
    return self._bins
  def GetNumBins(self):
    return len(self._bins)

  def SetMinCount(self,min):
    self._minCnt = min
  def GetMinCount(self):
    return self._minCnt

  def SetMaxCount(self,max):
    self._maxCnt = max
  def GetMaxCount(self):
    return self._maxCnt

  def SetShortestPathsOnly(self,val):
    if not val:
      raise ValueError,'only shortest paths signatures are currently supported'
    self._shortestPathsOnly = val
  def GetShortestPathsOnly(self):
    return self._shortestPathsOnly
    
  def SetIncludeBondOrder(self,val):
    self._includeBondOrder = val
  def GetIncludeBondOrder(self):
    return self._includeBondOrder
    
  
  def GetSignature(self,initialize=1):
    sig = Signature(patts=self.GetPatterns(),
                    bins=self.GetBins(),
                    labels = self.GetLabels(),
                    minCnt=self.GetMinCount(),
                    maxCnt=self.GetMaxCount(),
                    shortestPathsOnly=self.GetShortestPathsOnly(),
                    includeBondOrder=self.GetIncludeBondOrder())
    if initialize:
      sig.Init()
    return sig
    
