# $Id$
#
# Copyright (C) 2002-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" data structures for holding 2D pharmacophore signatures


  See Docs/Chem/Pharm2D.triangles.jpg for an illustration of the way
  pharmacophores are broken into triangles and labelled.

  See Docs/Chem/Pharm2D.signatures.jpg for an illustration of bit
  numbering

"""
import Chem
from Chem.Pharm2D import Utils
from DataStructs import SparseBitVect as BitVect

_verbose = 0

_sigVersion=100

class Pharm2DSig(object):
  """ 

  **Notes on Use**

    - After any changes, the Init() method should be called

  **Important Attributes**

   - patterns: a list of SmartsPatterns used to determine the features
      patterns can be initialized directly from the list or from a list of SMARTS
      strings

   - bins: the list of distance bins to be used.

   - minCount/maxCount: the minimum/maximum number of points to be included
      in a pharmacophore

   - shortestPathsOnly: toggles limiting the path-discovery algorithms to
      only find the shortest paths between 2 features.

   - includeBondOrder: toggles inclusion of bond order in consideration of
      shortest paths.

  **Limitations of Current Implementation**

   - All distances have the same numbers of bins

   - Only shortest-path matches are implemented

   
  """
  def __init__(self,patts=None,bins=None,labels=None,minCnt=2,maxCnt=4,
               shortestPathsOnly=1,includeBondOrder=0):
    self._sigVersion = _sigVersion
    self._bv = None
    self._patts = None
    if patts is not None:
      self._patts = patts[:]
    self._labels = None
    if labels is not None:
      self._labels = labels[:]
      
    self.bins = None
    if bins is not None:
      self.bins = bins[:]
    self._minCnt = minCnt
    self._maxCnt = maxCnt
    self._shortestPathsOnly = shortestPathsOnly
    self._includeBondOrder = includeBondOrder
    self._initLocals()

  def _initLocals(self):
    """ Internal use only

    """
    self._bv = None
    self._size = -1
    self._starts = {}
    self._scaffolds = []
    
    
  def __getstate__(self):
    """ used by the pickling machinery

    """
    res = {'_minCnt':self._minCnt,
           '_maxCnt':self._maxCnt,
           '_shortestPathsOnly':self._shortestPathsOnly,
           '_includeBondOrder':self._includeBondOrder,
           'bins': self.bins,
           '_bv':self._bv,
           '_labels':self._labels,
           '_sigVersion':self._sigVersion,
           }
    res['_patts'] = [Chem.MolToSmarts(x) for x in self._patts]

    return res
  def __setstate__(self,state):
    """ used by the pickling machinery

    """
    self.__dict__ = state
    patts = state['_patts']
    self.SetPatternsFromSmarts(patts)
    bv = self._bv
    self._initLocals()
    self._bv = bv
    try:
      self._sigVersion
    except AttributeError:
      self._sigVersion = _sigVersion
    self.Init(createBitVect=0)
  def __len__(self):
    return self.GetSize()
  def __getitem__(self,idx):
    if idx < 0 or idx >= self.GetSize():
      raise IndexError,'Index %d invalid'%(idx)
    return self._bv[idx]


  def SetPatterns(self,patts):
    self._patts = patts[:]
  def SetPatternsFromSmarts(self,smarts):
    import Chem
    self._patts = [None]*len(smarts)
    for i in range(len(smarts)):
      p = Chem.MolFromSmarts(smarts[i])
      self._patts[i] = p
  def GetPattern(self,which):
    return self._patts[which]
  def GetNumPatterns(self):
    return len(self._patts)

  def SetLabels(self,labels):
    self._labels = labels[:]
  def GetLabel(self,which):
    return self._labels[which]
    
  
  def SetBins(self,bins):
    """ bins should be a list of 2-tuples """
    self.bins = bins[:]
  def GetBin(self,which):
    return self.bins[which]
  def GetNumBins(self):
    return len(self.bins)
  def GetMinDist(self):
    return self.bins[0][0]
  def GetMaxDist(self):
    return self.bins[-1][1]
  

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
    
  def GetSize(self):
    return self._size
  
  def Init(self,createBitVect=1):
    """ Initializes internal parameters.  This **must** be called after
      making any changes to the signature

    **Arguments**

      - createBitVect: (optional) if this is nonzero, the bit vector
        used to store the on bits will be allocated.  Otherwise the
        existing bit vect will be reused (it better be big enough)

    """
    accum = 0
    self._scaffolds = [0]*(len(Utils.nPointDistDict[self.GetMaxCount()+1]))
    for i in range(self.GetMinCount(),self.GetMaxCount()+1):
      self._starts[i] = accum
      nDistsHere = len(Utils.nPointDistDict[i])
      scaffoldsHere = Utils.GetPossibleScaffolds(i,self.bins)
      nBitsHere = len(scaffoldsHere)
      self._scaffolds[nDistsHere] = scaffoldsHere
      pointsHere = Utils.NumCombinations(self.GetNumPatterns(),i) * nBitsHere
      
      accum += pointsHere
    self._size = accum
    if createBitVect:
      self._bv = BitVect(self._size)

    
  def _findBinIdx(self,dists,bins,scaffolds):
    """ OBSOLETE: this has been rewritten in C++
    Internal use only
     Returns the index of a bin defined by a set of distances.

    **Arguments**

      - dists: a sequence of distances (not binned)

      - bins: a sorted sequence of distance bins (2-tuples)

      - scaffolds: a list of possible scaffolds (bin combinations)

    **Returns**

      an integer bin index
      
    **Note**

      the value returned here is not an index in the overall
      signature.  It is, rather, an offset of a scaffold in the
      possible combinations of distance bins for a given
      proto-pharmacophore.
    
    """
    nBins = len(bins)
    nDists = len(dists)
    whichBins = [0]*nDists

    # This would be a ton easier if we had contiguous bins
    # i.e. if we could maintain the bins as a list of bounds)
    # because then we could use Python's bisect module.
    # Since we can't do that, we've got to do our own binary
    # search here.
    for i in range(nDists):
      dist = dists[i]
      where = -1

      # do a simple binary search:
      startP,endP = 0,len(bins)
      while startP<endP:
        midP = (startP+endP) // 2
        begBin,endBin = bins[midP]
        if dist < begBin:
          endP = midP
        elif dist >= endBin:
          startP = midP+1
        else:
          where = midP
          break
      if where < 0:
        return None
      whichBins[i] = where
    res = scaffolds.index(tuple(whichBins))
    return res

  def GetBitIdx(self,patts,dists,checkPatts=1):
    """ returns the index for a pharmacophore described using a set of
      patterns and distances

    **Arguments***

      - patts: a sequence of pattern indices

      - dists: a sequence of distance between the patterns, only the
        unique distances should be included, and they should be in the
        order defined in Utils.

      - checkPatts: (optional) if nonzero, the pattern vector is
        checked to ensure it is sorted


    **Returns**

      the integer bit index
      
    """
    nPoints = len(patts)
    if nPoints < self._minCnt: raise IndexError,'bad number of patterns'
    if nPoints > self._maxCnt: raise IndexError,'bad number of patterns'

    # this is the start of the nPoint-point pharmacophores
    startIdx = self._starts[nPoints]

    #
    # now we need to map the pattern indices to an offset from startIdx
    # 
    nPatts = len(self._patts)
    if checkPatts:
      tmp = list(patts)
      tmp.sort()
      if tmp!=list(patts):
        raise ValueError,'pattern vector not sorted'
    #patts = list(patts)
    #patts.sort()
    if patts[0]<0: raise IndexError,'bad pattern index'
    if max(patts)>=nPatts: raise IndexError,'bad pattern index'
    offset = Utils.CountUpTo(nPatts,nPoints,patts)
    if _verbose: print 'offset for patts %s: %d'%(str(patts),offset)
    offset *= len(self._scaffolds[len(dists)])
      
    try:
      if _verbose:
        print '>>>>>>>>>>>>>>>>>>>>>>>'
        print '\tScaffolds:',repr(self._scaffolds[len(dists)]),type(self._scaffolds[len(dists)])
        print '\tDists:',repr(dists),type(dists)
        print '\tbins:',repr(self.bins),type(self.bins)
      bin = self._findBinIdx(dists,self.bins,self._scaffolds[len(dists)])
    except ValueError:
      raise IndexError,'distance bin not found'

    return startIdx + offset + bin
    
  def SetBit(self,patts,dists,checkPatts=1):
    """ sets the bit defined by a collection of patterns and distances

    **Arguments***

      - patts: a sequence of pattern indices

      - dists: a sequence of distance between the patterns, only the
        unique distances should be included, and they should be in the
        order defined in Utils.

      - checkPatts: (optional) if nonzero, the pattern vector is
        checked to ensure it is sorted

    **Returns**

      the original status of the bit (whether or not it was set)
      
    """
    idx = self.GetBitIdx(patts,dists,checkPatts=checkPatts)
    if _verbose:
      print '*--> setting bit: %d'%(idx)
      print '\tfrom patts: %s and dists: %s\n'%(repr(patts),repr(dists))
    if idx >= self.GetSize():
      raise IndexError,'bad index (%d) calculated. %d is the max'%(idx,self.GetSize())
    return self._bv.SetBit(idx)

  def GetBit(self,patts,dists,checkPatts=1):
    """ returns the value of a 

    **Arguments***

      - patts: a sequence of pattern indices

      - dists: a sequence of distance between the patterns, only the
        unique distances should be included, and they should be in the
        order defined in Utils.

      - checkPatts: (optional) if nonzero, the pattern vector is
        checked to ensure it is sorted

    **Returns**

       whether or not the bit is set
       
    """
    idx = self.GetBitIdx(patts,dists,checkPatts=checkPatts)
    if idx >= self.GetSize():
      raise IndexError,'bad index (%d) calculated. %d is the max'%(idx,self.GetSize())
    return self._bv.GetBit(idx)

  def GetOnBits(self):
    """ returns our on bits 

    """
    return self._bv.GetOnBits()

  def GetBitInfo(self,idx):
    """ returns information about the given bit

     **Arguments**

       - idx: the bit index to be considered

     **Returns**

       a 3-tuple:

         1) the number of points in the pharmacophore

         2) the proto-pharmacophore (tuple of pattern indices)

         3) the scaffold (tuple of distance indices)
     
    """
    if idx >= self.GetSize():
      raise IndexError,'bad index (%d) queried. %d is the max'%(idx,self.GetSize())
    # first figure out how many points are in the p'cophore
    nPts = self.GetMinCount()
    while nPts < self.GetMaxCount() and self._starts[nPts+1]<=idx:
      nPts+=1

    # how far are we in from the start point?
    offsetFromStart = idx - self._starts[nPts]
    if _verbose:
      print '\t %d Points, %d offset'%(nPts,offsetFromStart)

    # lookup the number of scaffolds
    nDists = len(Utils.nPointDistDict[nPts])
    scaffolds = self._scaffolds[nDists]
      
    nScaffolds = len(scaffolds)

    # figure out to which proto-pharmacophore we belong:
    protoIdx = offsetFromStart / nScaffolds
    indexCombos = Utils.GetIndexCombinations(self.GetNumPatterns(),nPts)
    combo = indexCombos[protoIdx]
    if _verbose:
      print '\t combo: %s'%(str(combo))

    # and which scaffold:
    scaffoldIdx = offsetFromStart % nScaffolds
    scaffold = scaffolds[scaffoldIdx]
    if _verbose:
      print '\t scaffold: %s'%(str(scaffold))
    
    return nPts,combo,scaffold

  def GetBitDescription(self,bitIdx,includeBins=0,fullPage=1):
    """  returns HTML with a description of the bit

    **Arguments**

      - bitIdx: an integer bit index

      - includeBins: (optional) if nonzero, information about the bins will be
        included as well

      - fullPage: (optional) if nonzero, html headers and footers will
        be included (so as to make the output a complete page)

    **Returns**

      a string with the HTML

    """
    nPts,combo,scaffold = self.GetBitInfo(bitIdx)
    labels = [self._labels[x] for x in combo]
    dMat = zeros((nPts,nPts),Int)
    dVect = Utils.nPointDistDict[nPts]
    for idx in range(len(dVect)):
      i,j = dVect[idx]
      dMat[i,j] = scaffold[idx]
      dMat[j,i] = scaffold[idx]
    if fullPage:
      lines = ['<html><body>']
    else:
      lines = []
    lines.append("""<h2>Bit %d</h2>
    <p><b>Num Points:</b> %d
    """%(bitIdx,nPts))
    lines.append('<p><b>Distances</b><table border=1>')
    hdr = ' '.join(['<th>%s</th>'%x for x in labels])
    lines.append('<tr><td></td>%s</tr>'%(hdr))
    for i in range(nPts):
      row = ' '.join(['<td>%s</td>'%(str(dMat[i,x])) for x in range(nPts)])
      lines.append('<tr><th>%s</th>%s</tr>'%(labels[i],row))
    lines.append('</table>')

    if includeBins:
      lines.append('<p> <b>Distance Bin Information</b>')
      lines.append('<table border=1>')
      lines.append('<tr><td>bin</td><td>begin</td><td>end</td></tr>')
      for idx in range(self.GetNumBins()):
        beg,end = self.GetBin(idx)
        lines.append('<tr><td>%d</td><td>%d</td><td>%d</td></tr>'%(idx,beg,end))
      lines.append('</table>')
    if fullPage:
      lines.append("</body></html>")
    return '\n'.join(lines)



try:
  from Chem.Pharmacophores import cUtils
except ImportError:
  pass
else:
  Pharm2DSig._findBinIdx = cUtils.FindBinIdx


if __name__=='__main__':
  def test1():
    sig = Pharm2DSig()
    sig.SetPatternsFromSmarts(['O','N'])
    sig.SetBins([(1,2),(2,4),(4,8)])
    sig.SetMinCount(2)
    sig.SetMaxCount(3)
    sig.Init()
    print sig.GetSize()

  def test2():
    sig = Pharm2DSig()
    sig.SetPatternsFromSmarts(['O','N'])
    sig.SetBins([(0,2),(2,4),(4,8)])
    sig.SetMinCount(2)
    sig.SetMaxCount(3)
    sig.Init()
    vs = [((0,0),[1]),((1,1),[1]),((0,0),[2]),((0,0),[6]),((0,1),[1])]
    for patts,dist in vs:
      idx = sig.GetBitIdx(patts,dist)
      print patts,dist,idx
      
  def test4():
    import Chem
    import Generate
    sig = Pharm2DSig()
    sig.SetPatternsFromSmarts(['O'])
    sig.SetBins([(1,3),(3,4),(4,8)])
    sig.SetMinCount(2)
    sig.SetMaxCount(3)
    sig.Init()
    #print '---------c'
    #patts,dist = (0,0),[4]
    #idx = sig.GetBitIdx(patts,dist)
    #print patts,dist,idx
    mol = Chem.MolFromSmiles('OCCC1COCCO1')
    try:
      Generate.Gen2DFingerprint(mol,sig)
    except TypeError:
      import traceback
      traceback.print_exc()
    print '---------c'
    patts,dist = [0,0],[4]
    #idx = sig.GetBitIdx(patts,dist)
    #print patts,dist,idx
    sig.SetBit(patts,dist)




  def test3():
    sig = Pharm2DSig()
    sig.SetPatternsFromSmarts(['[OD1]','[OD2]','[ND2]','[N]'])
    sig.SetBins([(0,2),(2,4),(4,6),(6,8),(8,100)])
    sig.SetMinCount(2)
    sig.SetMaxCount(4)
    sig.Init()
    vs = [((0,0),[1]),((1,1),[1]),((0,0),[2]),((0,0),[6]),((0,1),[1]),((0,0,0),[1,1,1]),((0,0,0),[1,1,3]),
          ((0,0,0),[3,1,2]),((0,0,1),[1,1,1]),]
    for patts,dist in vs:
      print patts,dist,sig.GetBitIdx(patts,dist)

  test2()


    
