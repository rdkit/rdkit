#
#  Copyright (C) 2004 Rational Discovery LLC
#   All Rights Reserved
#
from Numeric import *
from DataStructs.TopNContainer import TopNContainer
import SVDPack
import PySVD

def showMat(mat):
  for row in mat:
    for col in row:
      print '% 6.3f'%col,
    print  


"""
throughout we use the notation of:
 Deerwester et al. J. Am. Soc. Inf. Sci. 41 391-407 (1990)
 
"""
class SimilarityCalculator(object):
  def __init__(self):
    self._vects = None
    self._idMap = None
    self._tForm = None
    self._DS = None
    self._T = None

  def SetVects(self,vects):
    """

    vects is a sequence of *sorted* sequences of bit IDs


    >>> calc = SimilarityCalculator()
    >>> calc.SetVects( ((1,2),(3,100),(1,2,2,4)) )
    >>> calc._vects
    ((0, 1), (2, 4), (0, 1, 3))
    >>> calc._vals
    (1, 1, 1, 1, 1, 2, 1)
    >>> calc._idMap[100]
    4
    >>> calc._idMap[4]
    3

    """
    self._idMap = {}
    self._vects = []
    self._vals = []
    self._tForm = None
    self._DS = None
    self._T = None

    tmpD = {}
    for vect in vects:
      for bit in vect:
        if not tmpD.has_key(bit):
          tmpD[bit] = 1
    ks = tmpD.keys()
    ks.sort()
    for i in range(len(ks)):
      self._idMap[ks[i]] = i
    ks = None
    tmpD = None
    
    for vect in vects:
      tmp = []
      i = 0
      nBits = len(vect)
      while i<nBits:
        bit = vect[i]
        # we need to do a few things with this bit:
        #  1) find its location in the global idMap
        #     (mapping bitID -> reduced space ID),
        #     adding a new location if necessary
        #  2) add an entry to the reduced-space vector
        #  3) add an entry to the reduced-space value array
        idx = self._idMap.get(bit)
        if idx < 0:
          idx = len(self._idMap)
          self._idMap[bit] = idx
        # update the reduced-space vector:
        tmp.append(idx)  
        # add an entry to the value array:
        self._vals.append(1)
        # now grab duplicates:
        j=i+1
        while(j<nBits and vect[j]==vect[i]):
          self._vals[-1] += 1
          j += 1
        i=j
      #if len(self._vects)==0:  
      #  print vect
      #  print tuple(tmp)
      self._vects.append(tuple(tmp))
    self._vects = tuple(self._vects)
    self._vals = tuple(self._vals)
          
  def UpdateSingularValues1(self,k=-1,cleanup=1):
    """
    >>> calc = SimilarityCalculator()
    >>> try:
    ...   calc.UpdateSingularValues()
    ... except ValueError:
    ...   ok=1
    ... else:
    ...   ok=0
    >>> ok
    1
    >>> calc.SetVects( ((0,2),(1,3),(0,1,2)) )
    >>> calc.UpdateSingularValues()
    >>> calc._S.shape[0]
    3

    Unless the optional cleanup argument is unset,the local vects
      (untransformed data points) will be destroyed after we're done
      with them here.  This can save significant memory:
    >>> try:
    ...   calc.UpdateSingularValues(2)
    ... except ValueError:
    ...   ok=1
    ... else:
    ...   ok=0
    >>> ok
    1
      
    Have to call SetVects again:
    >>> calc.SetVects( ((0,2),(1,3),(0,1,2)) )
    >>> calc.UpdateSingularValues(2)
    >>> calc._S.shape[0]
    2
    >>> print '%.4f'%calc._S[0]
    2.1889
    >>> print '%.4f'%calc._S[1]
    1.4142

    """
    if not self._vects:
      raise ValueError,"SetVects() not called"
    
    nRows = len(self._vects)
    nCols = len(self._idMap)
    if k==-1:
      k = min(nRows,nCols)
    #print self._vects
    #print self._vals
    #print 'shape:',nRows,nCols,len(self._vals)
    D,s,T =  PySVD.SparseSVD(self._vects,self._vals,nRows,nCols,k,1)
    T = transpose(T)
    D = transpose(D)
    # sometimes the rank of the matrix is smaller than we think
    # it should be:
    k = min(k,len(s))
    self._S = s[:k]
    #print 'k:',k,'T:',T.shape,'S:',s.shape,'D:',D.shape
    self._T = T
    self._DS = D*s

    if 0:
      print 'T:'
      showMat(self._T)

      print 'S:'
      print s

      print 'D:'
      showMat(D)


      print 'DS:'
      showMat(self._DS)

      print '------------------------'
      showMat(matrixmultiply(transpose(T),T))

      print '------------------------'
      showMat(matrixmultiply(transpose(D),D))

      print '------------------------'
      print self._T.shape,self._DS.shape
      showMat(matrixmultiply(self._T,transpose(self._DS)))
    

    # save some memory:
    if cleanup:
      self._vects = None
      self._vals = None
    
  def UpdateSingularValues(self,k=-1,cleanup=1):
    """


    >>> calc = SimilarityCalculator()
    >>> try:
    ...   calc.UpdateSingularValues()
    ... except ValueError:
    ...   ok=1
    ... else:
    ...   ok=0
    >>> ok
    1
    >>> calc.SetVects( ((0,2),(1,3),(0,1,2)) )
    >>> calc.UpdateSingularValues()
    >>> calc._S.shape[0]
    3

    Unless the optional cleanup argument is unset,the local vects
      (untransformed data points) will be destroyed after we're done
      with them here.  This can save significant memory:
    >>> try:
    ...   calc.UpdateSingularValues(2)
    ... except ValueError:
    ...   ok=1
    ... else:
    ...   ok=0
    >>> ok
    1
      
    Have to call SetVects again:
    >>> calc.SetVects( ((0,2),(1,3),(0,1,2)) )
    >>> calc.UpdateSingularValues(2)
    >>> calc._S.shape[0]
    2
    >>> print '%.4f'%calc._S[0]
    2.1889
    >>> print '%.4f'%calc._S[1]
    1.4142

    """
    if not self._vects:
      raise ValueError,"SetVects() not called"
    
    nRows = len(self._vects)
    nCols = len(self._idMap)
    if k==-1:
      k = min(nRows,nCols)

    SVDPack.DoSVD(self,k)
    
    # save some memory:
    if cleanup:
      self._vects = None
      self._vals = None
  def ForceSingularValues(self,k,T,D,s,cleanup=1):
    k = min(k,len(s))
    self._S = s[:k]
    self._T = T
    self._DS = D*s
    if cleanup:
      self._vects = None
      self._vals = None
    
  


  def PackPoint(self,pt):
    """

    converts a point from the normal space to the reduced (packed) space
    
    >>> calc = SimilarityCalculator()
    >>> calc.SetVects( ((1,2),(3,100),(1,2,2,4)) )
    >>> calc.PackPoint( (1,2) )
    array([ 1.,  1.,  0.,  0.,  0.])
    >>> calc.PackPoint( (1,2,2) )
    array([ 1.,  2.,  0.,  0.,  0.])
    >>> calc.PackPoint( (1,2,5) )
    array([ 1.,  1.,  0.,  0.,  0.])
    
    """
    if not self._idMap:
      raise ValueError,"SetVects() not called"
    res = zeros((len(self._idMap),),Float)
    for val in pt:
      idx = self._idMap.get(val,-1)
      if idx>=0:
        res[idx] += 1
    return res    

  def TransformPoint(self,pt):
    """ Transforms a point into the SVD space

    >>> calc = SimilarityCalculator()
    >>> calc.SetVects( ((0,2),(1,3),(0,1,2)) )
    >>> calc.UpdateSingularValues(k=3)

    if we pass in a point used for the SVD, we should just get the
    transformed version of that point back:
    >>> v2 = calc.TransformPoint( (0,2) )

    #>>> v1 = transpose(calc._singularVects[0])
    #>>> abs(max(v1-v2))<1e-6
    #1
    #>>> v1 = transpose(calc._singularVects[1])
    #>>> abs(max(v1-v2))>1e-6
    #1
    
    
    """
    if not self._T:
      raise ValueError,"UpdateSingularValues() not called"
    packed = self.PackPoint(pt)
    #print packed.shape,self._T.shape,self._DS.shape
    v = matrixmultiply(packed,self._T)
    return v
  
  def ScorePoint(self,pt,against=None,isTransformed=0,topN=0,
                 threshold=-1.0,
                 excludeThese=[]):
    """

    return value is a sequence of 2-tuples: (score, index)

    >>> calc = SimilarityCalculator()
    >>> calc.SetVects( ((0,2),(1,3),(0,1,2)) )
    >>> calc.UpdateSingularValues(k=3)
    >>> r = calc.ScorePoint((0,2),against=[0])[0]
    >>> print '%.2f'%r[0], r[1]
    1.00 0

    can transform the point in advance:
    >>> pt = calc.TransformPoint( (0,2) )
    >>> r = calc.ScorePoint(pt,against=[0],isTransformed=1)[0]
    >>> print '%.2f'%r[0], r[1]
    1.00 0

    default is to score against a variety of vectors at once:
    >>> [abs(x[0])>1e-4 for x in calc.ScorePoint(pt,isTransformed=1)]
    [True, False, True]

    >>> [abs(x[0])>1e-4 for x in calc.ScorePoint((0,3,6))]
    [True, True, True]

    >>> [abs(x[0])>1e-4 for x in calc.ScorePoint((0,3,6),topN=2)]
    [True, True]

    you can also put a threshold on the similarity metric:
    >>> len(calc.ScorePoint((0,3,6)))
    3
    >>> len(calc.ScorePoint((0,3,6),threshold=0.50))
    2



    "extra" bits (those that weren't in the training vectors) are
    ignored:
    >>> [abs(x[0])>1e-4 for x in calc.ScorePoint((0,2,12))]
    [True, False, True]
    
    # look at the indices:
    >>> v = [x[1] for x in calc.ScorePoint((0,3,6),topN=2)]
    >>> v.sort()
    >>> v
    [0, 1]
    >>> [x[1] for x in calc.ScorePoint((0,3,6),topN=2,excludeThese=[1])]
    [2, 0]

    
    """
    if not self._T or not self._DS:
      raise ValueError,"UpdateSingularValues() not called"
    if not isTransformed:
      pt = self.TransformPoint(pt)
    if against is None:
      against = range(self._DS.shape[0])
    try:
      nPts = len(against)
    except AttributeError:
      against = [against]
      nPts =1
    if topN:
      res = TopNContainer(topN)
    else:
      res = []
    ptSize = sqrt(dot(pt,pt))
    #indices = range(nPts)
    for idx in excludeThese:
      try:
        against.remove(idx)
      except ValueError:
        pass

    #print 'PT:',pt
    for i in against:
      v = self._DS[i]
      vSize = sqrt(dot(v,v))
      numer = dot(v,pt)
      denom = vSize*ptSize

      if denom != 0.0:
        simVal = numer/denom
      else:
        simVal = 0.0
      if simVal>threshold:
        if not topN:
          res.append((simVal,i))
        else:
          res.Insert(simVal,i)
    return res
  
#------------------------------------
#
#  doctest boilerplate
#
molTest="""

This is a nice test because not only is it molecular with vaguely
sensible results, but the matrix is only of rank 2 (despite being
3x3), so it can catch boundary conditions in the solver.

>>> import Chem
>>> from Chem.AtomPairs import Torsions,Utils
>>> m1 = Chem.MolFromSmiles('CCN(CCO)CC')
>>> m2 = Chem.MolFromSmiles('CCN(CCO)CCO')
>>> m3 = Chem.MolFromSmiles('OCCN(CCO)CCO')
>>> fp1 = Torsions.GetTopologicalTorsionFingerprintAsIds(m1)
>>> print fp1
>>> fp2 = Torsions.GetTopologicalTorsionFingerprintAsIds(m2)
>>> fp3 = Torsions.GetTopologicalTorsionFingerprintAsIds(m3)
>>> calc = SimilarityCalculator()
>>> calc.SetVects((fp1,fp2,fp3))
>>> calc.UpdateSingularValues()
>>> scores = calc.ScorePoint(fp1)
>>> print '%.3f'%scores[0][0],scores[0][1]
1.000 0
>>> print '%.3f'%scores[1][0],scores[1][1]
0.802 1
>>> print '%.3f'%scores[2][0],scores[2][1]
0.488 2
>>> scores = calc.ScorePoint(fp2)
>>> print '%.3f'%scores[0][0],scores[0][1]
0.802 0
>>> print '%.3f'%scores[1][0],scores[1][1]
1.000 1
>>> print '%.3f'%scores[2][0],scores[2][1]
0.913 2
>>> scores = calc.ScorePoint(fp3)
>>> print '%.3f'%scores[0][0],scores[0][1]
0.488 0
>>> print '%.3f'%scores[1][0],scores[1][1]
0.913 1
>>> print '%.3f'%scores[2][0],scores[2][1]
1.000 2
>>> scores = calc.ScorePoint(fp1,topN=2)
>>> scores.reverse()
>>> print '%.3f'%scores[0][0],scores[0][1]
1.000 0
>>> print '%.3f'%scores[1][0],scores[1][1]
0.802 1


"""

__test__={'molTest':molTest}

def _test():
  import doctest,sys
  return doctest.testmod(sys.modules["__main__"])


if __name__ == '__main__':
  import sys
  failed,tried = _test()
  sys.exit(failed)
  
