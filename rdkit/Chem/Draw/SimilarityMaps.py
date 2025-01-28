#
#  Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
#  Copyright (c) 2021, Greg Landrum
#  All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Created by Sereina Riniker, Aug 2013

import copy
import math

try:
  import matplotlib
  from matplotlib import cm
  from matplotlib.colors import LinearSegmentedColormap
except ImportError:
  cm = None
except RuntimeError:
  cm = None

import numpy

from rdkit import Chem, DataStructs, Geometry
from rdkit.Chem import Draw, rdDepictor
from rdkit.Chem import rdMolDescriptors as rdMD
from rdkit.Chem.Draw import rdMolDraw2D


def _DeleteFpInfoAttr(mol):
  if hasattr(mol, '_fpInfo'):
    delattr(mol, '_fpInfo')


def GetAtomicWeightsForFingerprint(refMol, probeMol, fpFunction, metric=DataStructs.DiceSimilarity):
  """
    Calculates the atomic weights for the probe molecule
    based on a fingerprint function and a metric.

    Parameters:
      refMol -- the reference molecule
      probeMol -- the probe molecule
      fpFunction -- the fingerprint function
      metric -- the similarity metric

    Note:
      If fpFunction needs additional parameters, use a lambda construct
    """
  _DeleteFpInfoAttr(probeMol)
  _DeleteFpInfoAttr(refMol)

  refFP = fpFunction(refMol, -1)
  baseSimilarity = metric(refFP, fpFunction(probeMol, -1))

  # loop over atoms
  weights = [
    baseSimilarity - metric(refFP, fpFunction(probeMol, atomId))
    for atomId in range(probeMol.GetNumAtoms())
  ]

  _DeleteFpInfoAttr(probeMol)
  _DeleteFpInfoAttr(refMol)
  return weights


def GetAtomicWeightsForModel(probeMol, fpFunction, predictionFunction):
  """
    Calculates the atomic weights for the probe molecule based on
    a fingerprint function and the prediction function of a ML model.

    Parameters:
      probeMol -- the probe molecule
      fpFunction -- the fingerprint function
      predictionFunction -- the prediction function of the ML model
    """

  _DeleteFpInfoAttr(probeMol)
  baseProba = predictionFunction(fpFunction(probeMol, -1))

  # loop over atoms
  weights = [
    baseProba - predictionFunction(fpFunction(probeMol, atomId))
    for atomId in range(probeMol.GetNumAtoms())
  ]

  _DeleteFpInfoAttr(probeMol)
  return weights


def GetStandardizedWeights(weights):
  """
    Normalizes the weights,
    such that the absolute maximum weight equals 1.0.

    Parameters:
      weights -- the list with the atomic weights
    """

  maxAbsWeight = max(math.fabs(w) for w in weights)

  if maxAbsWeight > 0:
    return [w / maxAbsWeight for w in weights], maxAbsWeight
  return weights, maxAbsWeight


def GetSimilarityMapFromWeights(mol, weights, draw2d, colorMap=None, scale=-1, size=(250, 250),
                                sigma=None, coordScale=1.5, step=0.01, colors='k', contourLines=10,
                                alpha=0.5, **kwargs):
  """
    Generates the similarity map for a molecule given the atomic weights.

    Parameters:
      mol -- the molecule of interest
      weights -- the weights to use
      draw2d -- the RDKit drawing object
      colorMap -- the matplotlib color map scheme, default is custom PiWG color map
      scale -- the scaling: scale < 0 -> the absolute maximum weight is used as maximum scale
                            scale = double -> this is the maximum scale
      size -- the size of the figure
      sigma -- the sigma for the Gaussians
      coordScale -- scaling factor for the coordinates
      step -- the step for calcAtomGaussian
      colors -- color of the contour lines
      contourLines -- if integer number N: N contour lines are drawn
                      if list(numbers): contour lines at these numbers are drawn
      alpha -- the alpha blending value for the contour lines
      kwargs -- additional arguments for drawing
    """
  if mol.GetNumAtoms() < 2:
    raise ValueError("too few atoms")

  if draw2d is None:
    raise ValueError("the draw2d argument must be provided")
  mol = rdMolDraw2D.PrepareMolForDrawing(mol, addChiralHs=False)
  if not mol.GetNumConformers():
    rdDepictor.Compute2DCoords(mol)
  if sigma is None:
    if mol.GetNumBonds() > 0:
      bond = mol.GetBondWithIdx(0)
      idx1 = bond.GetBeginAtomIdx()
      idx2 = bond.GetEndAtomIdx()
      sigma = 0.3 * (mol.GetConformer().GetAtomPosition(idx1) -
                     mol.GetConformer().GetAtomPosition(idx2)).Length()
    else:
      sigma = 0.3 * (mol.GetConformer().GetAtomPosition(0) -
                     mol.GetConformer().GetAtomPosition(1)).Length()
    sigma = round(sigma, 2)

  sigmas = [sigma] * mol.GetNumAtoms()
  locs = []
  for i in range(mol.GetNumAtoms()):
    p = mol.GetConformer().GetAtomPosition(i)
    locs.append(Geometry.Point2D(p.x, p.y))
  draw2d.ClearDrawing()
  ps = Draw.ContourParams()
  ps.fillGrid = True
  ps.gridResolution = 0.1
  ps.extraGridPadding = 0.5

  if colorMap is not None:
    if cm is not None and isinstance(colorMap, type(cm.Blues)):
      # it's a matplotlib colormap:
      clrs = [tuple(x) for x in colorMap([0, 0.5, 1])]
    elif type(colorMap) == str:
      if cm is None:
        raise ValueError("cannot provide named colormaps unless matplotlib is installed")
      clrs = [tuple(x) for x in matplotlib.colormaps[colorMap]([0, 0.5, 1])]
    else:
      clrs = [colorMap[0], colorMap[1], colorMap[2]]
    ps.setColourMap(clrs)

  Draw.ContourAndDrawGaussians(draw2d, locs, weights, sigmas, nContours=contourLines, params=ps)
  draw2d.drawOptions().clearBackground = False
  draw2d.DrawMolecule(mol)
  return draw2d


def GetSimilarityMapForFingerprint(refMol, probeMol, fpFunction, draw2d,
                                   metric=DataStructs.DiceSimilarity, **kwargs):
  """
    Generates the similarity map for a given reference and probe molecule,
    fingerprint function and similarity metric.

    Parameters:
      refMol -- the reference molecule
      probeMol -- the probe molecule
      fpFunction -- the fingerprint function
      metric -- the similarity metric.
      kwargs -- additional arguments for drawing
    """

  weights = GetAtomicWeightsForFingerprint(refMol, probeMol, fpFunction, metric)
  weights, maxWeight = GetStandardizedWeights(weights)
  draw2d = GetSimilarityMapFromWeights(probeMol, weights, draw2d, **kwargs)
  return draw2d, maxWeight


def GetSimilarityMapForModel(probeMol, fpFunction, predictionFunction, draw2d, **kwargs):
  """
    Generates the similarity map for a given ML model and probe molecule,
    and fingerprint function.

    Parameters:
      probeMol -- the probe molecule
      fpFunction -- the fingerprint function
      predictionFunction -- the prediction function of the ML model
      kwargs -- additional arguments for drawing
    """

  weights = GetAtomicWeightsForModel(probeMol, fpFunction, predictionFunction)
  weights, maxWeight = GetStandardizedWeights(weights)
  fig = GetSimilarityMapFromWeights(probeMol, weights, draw2d, **kwargs)
  return fig, maxWeight


apDict = {}
apDict['normal'] = lambda m, bits, minl, maxl, bpe, ia, **kwargs: rdMD.GetAtomPairFingerprint(
  m, minLength=minl, maxLength=maxl, ignoreAtoms=ia, **kwargs)
apDict['hashed'] = lambda m, bits, minl, maxl, bpe, ia, **kwargs: rdMD.GetHashedAtomPairFingerprint(
  m, nBits=bits, minLength=minl, maxLength=maxl, ignoreAtoms=ia, **kwargs)
apDict[
  'bv'] = lambda m, bits, minl, maxl, bpe, ia, **kwargs: rdMD.GetHashedAtomPairFingerprintAsBitVect(
    m, nBits=bits, minLength=minl, maxLength=maxl, nBitsPerEntry=bpe, ignoreAtoms=ia, **kwargs)


# usage:   lambda m,i: GetAPFingerprint(m, i, fpType, nBits, minLength, maxLength, nBitsPerEntry)
def GetAPFingerprint(mol, atomId=-1, fpType='normal', nBits=2048, minLength=1, maxLength=30,
                     nBitsPerEntry=4, **kwargs):
  """
    Calculates the atom pairs fingerprint with the torsions of atomId removed.

    Parameters:
      mol -- the molecule of interest
      atomId -- the atom to remove the pairs for (if -1, no pair is removed)
      fpType -- the type of AP fingerprint ('normal', 'hashed', 'bv')
      nBits -- the size of the bit vector (only for fpType='bv')
      minLength -- the minimum path length for an atom pair
      maxLength -- the maxmimum path length for an atom pair
      nBitsPerEntry -- the number of bits available for each pair
    """
  if fpType not in ['normal', 'hashed', 'bv']:
    raise ValueError("Unknown Atom pairs fingerprint type")
  if atomId < 0:
    return apDict[fpType](mol, nBits, minLength, maxLength, nBitsPerEntry, 0, **kwargs)
  if atomId >= mol.GetNumAtoms():
    raise ValueError("atom index greater than number of atoms")
  return apDict[fpType](mol, nBits, minLength, maxLength, nBitsPerEntry, [atomId], **kwargs)


ttDict = {}
ttDict['normal'] = lambda m, bits, ts, bpe, ia, **kwargs: rdMD.GetTopologicalTorsionFingerprint(
  m, targetSize=ts, ignoreAtoms=ia, **kwargs)
ttDict[
  'hashed'] = lambda m, bits, ts, bpe, ia, **kwargs: rdMD.GetHashedTopologicalTorsionFingerprint(
    m, nBits=bits, targetSize=ts, ignoreAtoms=ia, **kwargs)
ttDict[
  'bv'] = lambda m, bits, ts, bpe, ia, **kwargs: rdMD.GetHashedTopologicalTorsionFingerprintAsBitVect(
    m, nBits=bits, targetSize=ts, nBitsPerEntry=bpe, ignoreAtoms=ia, **kwargs)


# usage:   lambda m,i: GetTTFingerprint(m, i, fpType, nBits, targetSize)
def GetTTFingerprint(mol, atomId=-1, fpType='normal', nBits=2048, targetSize=4, nBitsPerEntry=4,
                     **kwargs):
  """
    Calculates the topological torsion fingerprint with the pairs of atomId removed.

    Parameters:
      mol -- the molecule of interest
      atomId -- the atom to remove the torsions for (if -1, no torsion is removed)
      fpType -- the type of TT fingerprint ('normal', 'hashed', 'bv')
      nBits -- the size of the bit vector (only for fpType='bv')
      minLength -- the minimum path length for an atom pair
      maxLength -- the maxmimum path length for an atom pair
      nBitsPerEntry -- the number of bits available for each torsion

    any additional keyword arguments will be passed to the fingerprinting function.

    """
  if fpType not in ['normal', 'hashed', 'bv']:
    raise ValueError("Unknown Topological torsion fingerprint type")
  if atomId < 0:
    return ttDict[fpType](mol, nBits, targetSize, nBitsPerEntry, 0, **kwargs)
  if atomId >= mol.GetNumAtoms():
    raise ValueError("atom index greater than number of atoms")
  return ttDict[fpType](mol, nBits, targetSize, nBitsPerEntry, [atomId], **kwargs)


# usage:   lambda m,i: GetMorganFingerprint(m, i, radius, fpType, nBits, useFeatures)
def GetMorganFingerprint(mol, atomId=-1, radius=2, fpType='bv', nBits=2048, useFeatures=False,
                         **kwargs):
  """
    Calculates the Morgan fingerprint with the environments of atomId removed.

    Parameters:
      mol -- the molecule of interest
      radius -- the maximum radius
      fpType -- the type of Morgan fingerprint: 'count' or 'bv'
      atomId -- the atom to remove the environments for (if -1, no environments is removed)
      nBits -- the size of the bit vector (only for fpType = 'bv')
      useFeatures -- if false: ConnectivityMorgan, if true: FeatureMorgan

    any additional keyword arguments will be passed to the fingerprinting function.
    """
  if fpType not in ['bv', 'count']:
    raise ValueError("Unknown Morgan fingerprint type")

  isBitVect = fpType == 'bv'
  if not hasattr(mol, '_fpInfo'):
    info = {}
    # get the fingerprint & construct the bitmap template
    if isBitVect:
      molFp = rdMD.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits, useFeatures=useFeatures,
                                                 bitInfo=info, **kwargs)
      bitmap = [DataStructs.ExplicitBitVect(nBits) for _ in range(mol.GetNumAtoms())]

    else:
      molFp = rdMD.GetMorganFingerprint(mol, radius, useFeatures=useFeatures, bitInfo=info,
                                        **kwargs)
      bitmap = [[] for _ in range(mol.GetNumAtoms())]

    # construct the bit map
    for bit, es in info.items():
      for at1, rad in es:
        if rad == 0:  # for radius 0
          if isBitVect:
            bitmap[at1][bit] = 1
          else:
            bitmap[at1].append(bit)
        else:  # for radii > 0
          env = Chem.FindAtomEnvironmentOfRadiusN(mol, rad, at1)
          amap = {}
          Chem.PathToSubmol(mol, env, atomMap=amap)
          for at2 in amap.keys():
            if isBitVect:
              bitmap[at2][bit] = 1
            else:
              bitmap[at2].append(bit)
    mol._fpInfo = (molFp, bitmap)

  if atomId < 0:
    return mol._fpInfo[0]
  # remove the bits of atomId
  if atomId >= mol.GetNumAtoms():
    raise ValueError("atom index greater than number of atoms")
  if len(mol._fpInfo) != 2:
    raise ValueError("_fpInfo not set")
  if isBitVect:
    molFp = mol._fpInfo[0] ^ mol._fpInfo[1][atomId]  # xor
  else:  # count
    molFp = copy.deepcopy(mol._fpInfo[0])
    # delete the bits with atomId
    for bit in mol._fpInfo[1][atomId]:
      molFp[bit] -= 1
  return molFp


# usage:   lambda m,i: GetRDKFingerprint(m, i, fpType, nBits, minPath, maxPath, nBitsPerHash)
def GetRDKFingerprint(mol, atomId=-1, fpType='bv', nBits=2048, minPath=1, maxPath=5, nBitsPerHash=2,
                      **kwargs):
  """
    Calculates the RDKit fingerprint with the paths of atomId removed.

    Parameters:
      mol -- the molecule of interest
      atomId -- the atom to remove the paths for (if -1, no path is removed)
      fpType -- the type of RDKit fingerprint: 'bv'
      nBits -- the size of the bit vector
      minPath -- minimum path length
      maxPath -- maximum path length
      nBitsPerHash -- number of to set per path
    """
  if fpType not in ['bv', '']:
    raise ValueError("Unknown RDKit fingerprint type")

  fpType = 'bv'
  if not hasattr(mol, '_fpInfo'):
    info = []  # list with bits for each atom
    # get the fingerprint
    molFp = Chem.RDKFingerprint(mol, fpSize=nBits, minPath=minPath, maxPath=maxPath,
                                nBitsPerHash=nBitsPerHash, atomBits=info, **kwargs)
    mol._fpInfo = (molFp, info)

  if atomId < 0:
    return mol._fpInfo[0]

  # remove the bits of atomId
  if atomId >= mol.GetNumAtoms():
    raise ValueError("atom index greater than number of atoms")
  if len(mol._fpInfo) != 2:
    raise ValueError("_fpInfo not set")

  molFp = copy.deepcopy(mol._fpInfo[0])
  molFp.UnSetBitsFromList(mol._fpInfo[1][atomId])
  return molFp
