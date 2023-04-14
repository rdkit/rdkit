#
# Copyright (C) 2003-2008 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" contains factory class for producing signatures


"""
import copy

import numpy

from rdkit.Chem.Pharm2D import Utils
from rdkit.DataStructs import (IntSparseIntVect, LongSparseIntVect,
                               SparseBitVect)

_verbose = False


class SigFactory(object):
  """

      SigFactory's are used by creating one, setting the relevant
      parameters, then calling the GetSignature() method each time a
      signature is required.

    """

  def __init__(self, featFactory, useCounts=False, minPointCount=2, maxPointCount=3,
               shortestPathsOnly=True, includeBondOrder=False, skipFeats=None,
               trianglePruneBins=True):
    self.featFactory = featFactory
    self.useCounts = useCounts
    self.minPointCount = minPointCount
    self.maxPointCount = maxPointCount
    self.shortestPathsOnly = shortestPathsOnly
    self.includeBondOrder = includeBondOrder
    self.trianglePruneBins = trianglePruneBins
    if skipFeats is None:
      self.skipFeats = []
    else:
      self.skipFeats = skipFeats
    self._bins = None
    self.sigKlass = None

  def SetBins(self, bins):
    """ bins should be a list of 2-tuples """
    self._bins = copy.copy(bins)
    self.Init()

  def GetBins(self):
    return self._bins

  def GetNumBins(self):
    return len(self._bins)

  def GetSignature(self):
    return self.sigKlass(self._sigSize)

  def _GetBitSummaryData(self, bitIdx):
    nPts, combo, scaffold = self.GetBitInfo(bitIdx)
    fams = self.GetFeatFamilies()
    labels = [fams[x] for x in combo]
    dMat = numpy.zeros((nPts, nPts), dtype=numpy.int64)
    dVect = Utils.nPointDistDict[nPts]
    for idx in range(len(dVect)):
      i, j = dVect[idx]
      dMat[i, j] = scaffold[idx]
      dMat[j, i] = scaffold[idx]

    return nPts, combo, scaffold, labels, dMat

  def GetBitDescriptionAsText(self, bitIdx, includeBins=0, fullPage=1):
    """  returns text with a description of the bit

        **Arguments**

          - bitIdx: an integer bit index

          - includeBins: (optional) if nonzero, information about the bins will be
            included as well

          - fullPage: (optional) if nonzero, html headers and footers will
            be included (so as to make the output a complete page)

        **Returns**

          a string with the HTML

        """
    raise NotImplementedError('Missing implementation')

  def GetBitDescription(self, bitIdx):
    """  returns a text description of the bit

        **Arguments**

          - bitIdx: an integer bit index

        **Returns**

          a string

        """
    _, _, _, labels, dMat = self._GetBitSummaryData(bitIdx)
    res = " ".join(labels) + " "
    for row in dMat:
      res += "|" + " ".join([str(x) for x in row])
    res += "|"
    return res

  def _findBinIdx(self, dists, bins, scaffolds):
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

    whichBins = [0] * len(dists)
    # This would be a ton easier if we had contiguous bins
    # i.e. if we could maintain the bins as a list of bounds)
    # because then we could use Python's bisect module.
    # Since we can't do that, we've got to do our own binary
    # search here.
    for i, dist in enumerate(dists):
      where = -1
      # do a simple binary search:
      startP, endP = 0, len(bins)
      while startP < endP:
        midP = (startP + endP) // 2
        begBin, endBin = bins[midP]
        if dist < begBin:
          endP = midP
        elif dist >= endBin:
          startP = midP + 1
        else:
          where = midP
          break
      if where < 0:
        return None
      whichBins[i] = where
    res = scaffolds.index(tuple(whichBins))
    if _verbose:
      print('----- _fBI  -----------')
      print(' scaffolds:', scaffolds)
      print(' bins:', whichBins)
      print(' res:', res)
    return res

  def GetFeatFamilies(self):
    fams = [fam for fam in self.featFactory.GetFeatureFamilies() if fam not in self.skipFeats]
    fams.sort()
    return fams

  def GetMolFeats(self, mol):
    featFamilies = self.GetFeatFamilies()
    featMatches = {}
    for fam in featFamilies:
      feats = self.featFactory.GetFeaturesForMol(mol, includeOnly=fam)
      featMatches[fam] = [feat.GetAtomIds() for feat in feats]

    return [featMatches[x] for x in featFamilies]

  def GetBitIdx(self, featIndices, dists, sortIndices=True):
    """ returns the index for a pharmacophore described using a set of
          feature indices and distances

        **Arguments***

          - featIndices: a sequence of feature indices

          - dists: a sequence of distance between the features, only the
            unique distances should be included, and they should be in the
            order defined in Utils.

          - sortIndices : sort the indices

        **Returns**

          the integer bit index

        """
    nPoints = len(featIndices)
    if nPoints > 3:
      raise NotImplementedError('>3 points not supported')
    if nPoints < self.minPointCount:
      raise IndexError('bad number of points')
    if nPoints > self.maxPointCount:
      raise IndexError('bad number of points')

    # this is the start of the nPoint-point pharmacophores
    startIdx = self._starts[nPoints]

    #
    # now we need to map the pattern indices to an offset from startIdx
    #
    if sortIndices:
      tmp = list(featIndices)
      tmp.sort()
      featIndices = tmp

    if featIndices[0] < 0:
      raise IndexError('bad feature index')
    if max(featIndices) >= self._nFeats:
      raise IndexError('bad feature index')

    if nPoints == 3:
      featIndices, dists = Utils.OrderTriangle(featIndices, dists)

    offset = Utils.CountUpTo(self._nFeats, nPoints, featIndices)
    if _verbose:
      print(f'offset for feature {str(featIndices)}: {offset}')
    offset *= len(self._scaffolds[len(dists)])

    try:
      if _verbose:
        print('>>>>>>>>>>>>>>>>>>>>>>>')
        print('\tScaffolds:', repr(self._scaffolds[len(dists)]), type(self._scaffolds[len(dists)]))
        print('\tDists:', repr(dists), type(dists))
        print('\tbins:', repr(self._bins), type(self._bins))
      bin_ = self._findBinIdx(dists, self._bins, self._scaffolds[len(dists)])
    except ValueError:
      fams = self.GetFeatFamilies()
      fams = [fams[x] for x in featIndices]
      raise IndexError('distance bin not found: feats: %s; dists=%s; bins=%s; scaffolds: %s' %
                       (fams, dists, self._bins, self._scaffolds))

    return startIdx + offset + bin_

  def GetBitInfo(self, idx):
    """ returns information about the given bit

         **Arguments**

           - idx: the bit index to be considered

         **Returns**

           a 3-tuple:

             1) the number of points in the pharmacophore

             2) the proto-pharmacophore (tuple of pattern indices)

             3) the scaffold (tuple of distance indices)

        """
    if idx >= self._sigSize:
      raise IndexError(f'bad index ({idx}) queried. {self._sigSize} is the max')
    # first figure out how many points are in the p'cophore
    nPts = self.minPointCount
    while nPts < self.maxPointCount and self._starts[nPts + 1] <= idx:
      nPts += 1

    # how far are we in from the start point?
    offsetFromStart = idx - self._starts[nPts]
    if _verbose:
      print(f'\t {nPts} Points, {offsetFromStart} offset')

    # lookup the number of scaffolds
    nDists = len(Utils.nPointDistDict[nPts])
    scaffolds = self._scaffolds[nDists]

    nScaffolds = len(scaffolds)

    # figure out to which proto-pharmacophore we belong:
    protoIdx = offsetFromStart // nScaffolds
    indexCombos = Utils.GetIndexCombinations(self._nFeats, nPts)
    combo = tuple(indexCombos[protoIdx])
    if _verbose:
      print(f'\t combo: {str(combo)}')

    # and which scaffold:
    scaffoldIdx = offsetFromStart % nScaffolds
    scaffold = scaffolds[scaffoldIdx]
    if _verbose:
      print(f'\t scaffold: {str(scaffold)}')
    return nPts, combo, scaffold

  def Init(self):
    """ Initializes internal parameters.  This **must** be called after
          making any changes to the signature parameters

        """
    accum = 0
    self._scaffolds = [0] * (len(Utils.nPointDistDict[self.maxPointCount + 1]))
    self._starts = {}
    if not self.skipFeats:
      self._nFeats = len(self.featFactory.GetFeatureFamilies())
    else:
      self._nFeats = 0
      for fam in self.featFactory.GetFeatureFamilies():
        if fam not in self.skipFeats:
          self._nFeats += 1
    for i in range(self.minPointCount, self.maxPointCount + 1):
      self._starts[i] = accum
      nDistsHere = len(Utils.nPointDistDict[i])
      scaffoldsHere = Utils.GetPossibleScaffolds(i, self._bins,
                                                 useTriangleInequality=self.trianglePruneBins)
      self._scaffolds[nDistsHere] = scaffoldsHere
      accum += (Utils.NumCombinations(self._nFeats, i) * len(scaffoldsHere))
    self._sigSize = accum
    if not self.useCounts:
      self.sigKlass = SparseBitVect
    elif self._sigSize < 2**31:
      self.sigKlass = IntSparseIntVect
    else:
      self.sigKlass = LongSparseIntVect

  def GetSigSize(self):
    return self._sigSize


try:
  from rdkit.Chem.Pharmacophores import cUtils
except ImportError:
  pass
else:
  SigFactory._findBinIdx = cUtils.FindBinIdx
