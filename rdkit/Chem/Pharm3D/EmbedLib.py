#
# Copyright (C) 2004-2008 Greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

import math
import sys
import time

import numpy

import rdkit.DistanceGeometry as DG
from rdkit import Chem
from rdkit import RDLogger as logging
from rdkit.Chem import ChemicalFeatures, ChemicalForceFields
from rdkit.Chem import rdDistGeom as MolDG
from rdkit.Chem.Pharm3D import ExcludedVolume
from rdkit.ML.Data import Stats

_times = {}

logger = logging.logger()
defaultFeatLength = 2.0


def GetAtomHeavyNeighbors(atom):
  """ returns a list of the heavy-atom neighbors of the
  atom passed in:

  >>> m = Chem.MolFromSmiles('CCO')
  >>> l = GetAtomHeavyNeighbors(m.GetAtomWithIdx(0))
  >>> len(l)
  1
  >>> isinstance(l[0],Chem.Atom)
  True
  >>> l[0].GetIdx()
  1

  >>> l = GetAtomHeavyNeighbors(m.GetAtomWithIdx(1))
  >>> len(l)
  2
  >>> l[0].GetIdx()
  0
  >>> l[1].GetIdx()
  2

  """
  return [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]


def ReplaceGroup(match, bounds, slop=0.01, useDirs=False, dirLength=defaultFeatLength):
  """ Adds an entry at the end of the bounds matrix for a point at
   the center of a multi-point feature

   returns a 2-tuple:
     new bounds mat
     index of point added

   >>> boundsMat = numpy.array([[0.0, 2.0, 2.0],[1.0, 0.0, 2.0],[1.0, 1.0, 0.0]])
   >>> match = [0, 1, 2]
   >>> bm,idx = ReplaceGroup(match, boundsMat, slop=0.0)

   the index is at the end:

   >>> idx == 3
   True

   and the matrix is one bigger:

   >>> bm.shape == (4, 4)
   True

   but the original bounds mat is not altered:

   >>> boundsMat.shape == (3, 3)
   True


   We make the assumption that the points of the
   feature form a regular polygon, are listed in order
   (i.e. pt 0 is a neighbor to pt 1 and pt N-1)
   and that the replacement point goes at the center:

   >>> print(', '.join([f'{x:.3f}' for x in bm[-1]]))
   0.577, 0.577, 0.577, 0.000
   >>> print(', '.join([f'{x:.3f}' for x in bm[:,-1]]))
   1.155, 1.155, 1.155, 0.000

   The slop argument (default = 0.01) is fractional:

   >>> bm, idx = ReplaceGroup(match, boundsMat)
   >>> print(', '.join([f'{x:.3f}' for x in bm[-1]]))
   0.572, 0.572, 0.572, 0.000
   >>> print(', '.join([f'{x:.3f}' for x in bm[:,-1]]))
   1.166, 1.166, 1.166, 0.000

  """
  maxVal = -1000.0
  minVal = 1e8
  nPts = len(match)
  for i in range(nPts):
    idx0 = match[i]
    if i < nPts - 1:
      idx1 = match[i + 1]
    else:
      idx1 = match[0]
    if idx1 < idx0:
      idx0, idx1 = idx1, idx0
    minVal = min(minVal, bounds[idx1, idx0])
    maxVal = max(maxVal, bounds[idx0, idx1])
  maxVal *= (1 + slop)
  minVal *= (1 - slop)

  scaleFact = 1.0 / (2.0 * math.sin(math.pi / nPts))
  minVal *= scaleFact
  maxVal *= scaleFact

  replaceIdx = bounds.shape[0]
  enhanceSize: int = int(bool(useDirs))  # Whether to increase the size of the bounds matrix by one
  bm = numpy.zeros((bounds.shape[0] + 1 + enhanceSize, bounds.shape[1] + 1 + enhanceSize),
                   dtype=numpy.float64)
  bm[:bounds.shape[0], :bounds.shape[1]] = bounds
  bm[:replaceIdx, replaceIdx] = 1000.

  if useDirs:
    bm[:replaceIdx + 1, replaceIdx + 1] = 1000.
    # set the feature - direction point bounds:
    bm[replaceIdx, replaceIdx + 1] = dirLength + slop
    bm[replaceIdx + 1, replaceIdx] = dirLength - slop

  for idx1 in match:
    bm[idx1, replaceIdx] = maxVal
    bm[replaceIdx, idx1] = minVal
    if useDirs:
      # set the point - direction point bounds:
      bm[idx1, replaceIdx + 1] = numpy.sqrt(bm[replaceIdx, replaceIdx + 1]**2 + maxVal**2)
      bm[replaceIdx + 1, idx1] = numpy.sqrt(bm[replaceIdx + 1, replaceIdx]**2 + minVal**2)
  return bm, replaceIdx


def EmbedMol(mol, bm, atomMatch=None, weight=2.0, randomSeed=-1, excludedVolumes=None):
  """  Generates an embedding for a molecule based on a bounds matrix and adds
  a conformer (id 0) to the molecule

  if the optional argument atomMatch is provided, it will be used to provide
  supplemental weights for the embedding routine (used in the optimization
  phase to ensure that the resulting geometry really does satisfy the
  pharmacophore).

  if the excludedVolumes is provided, it should be a sequence of
  ExcludedVolume objects

  >>> m = Chem.MolFromSmiles('c1ccccc1C')
  >>> bounds = MolDG.GetMoleculeBoundsMatrix(m)
  >>> bounds.shape == (7, 7)
  True
  >>> m.GetNumConformers()
  0
  >>> EmbedMol(m,bounds,randomSeed=23)
  >>> m.GetNumConformers()
  1


  """
  nAts = mol.GetNumAtoms()
  weights = []
  if atomMatch:
    atomMatchSize = len(atomMatch)
    weights = [(i, j, weight) for i in range(atomMatchSize) for j in range(i + 1, atomMatchSize)]

  if excludedVolumes:
    for vol in excludedVolumes:
      idx = vol.index
      # excluded volumes affect every other atom:
      for i in range(nAts):
        weights.append((i, idx, weight))
  coords = DG.EmbedBoundsMatrix(bm, weights=weights, numZeroFail=1, randomSeed=randomSeed)
  # for row in coords:
  #  print(', '.join(['%.2f'%x for x in row]))

  conf = Chem.Conformer(nAts)
  conf.SetId(0)
  for i in range(nAts):
    conf.SetAtomPosition(i, list(coords[i]))
  if excludedVolumes:
    for vol in excludedVolumes:
      vol.pos = numpy.array(coords[vol.index])

    print('   % 7.4f   % 7.4f   % 7.4f Ar  0  0  0  0  0  0  0  0  0  0  0  0' % tuple(coords[-1]),
          file=sys.stderr)
  mol.AddConformer(conf)


def AddExcludedVolumes(bm, excludedVolumes, smoothIt=True):
  """ Adds a set of excluded volumes to the bounds matrix
  and returns the new matrix

  excludedVolumes is a list of ExcludedVolume objects


   >>> boundsMat = numpy.array([[0.0, 2.0, 2.0],[1.0, 0.0, 2.0],[1.0, 1.0, 0.0]])
   >>> ev1 = ExcludedVolume.ExcludedVolume(([(0, ), 0.5, 1.0], ), exclusionDist=1.5)
   >>> bm = AddExcludedVolumes(boundsMat, (ev1, ))

   the results matrix is one bigger:

   >>> bm.shape == (4, 4)
   True

   and the original bounds mat is not altered:

   >>> boundsMat.shape == (3, 3)
   True

   >>> print(', '.join([f'{x:.3f}' for x in bm[-1]]))
   0.500, 1.500, 1.500, 0.000
   >>> print(', '.join([f'{x:.3f}' for x in bm[:,-1]]))
   1.000, 3.000, 3.000, 0.000

  """
  oDim = bm.shape[0]
  dim = oDim + len(excludedVolumes)
  res = numpy.zeros((dim, dim), dtype=numpy.float64)
  res[:oDim, :oDim] = bm
  for i, vol in enumerate(excludedVolumes):
    bmIdx = oDim + i
    vol.index = bmIdx

    # set values to all the atoms:
    res[bmIdx, :bmIdx] = vol.exclusionDist
    res[:bmIdx, bmIdx] = 1000.0

    # set values to our defining features:
    for indices, minV, maxV in vol.featInfo:
      for index in indices:
        try:
          res[bmIdx, index] = minV
          res[index, bmIdx] = maxV
        except IndexError:
          logger.error(f'BAD INDEX: res[{bmIdx},{index}], shape is {str(res.shape)}')
          raise IndexError

    # set values to other excluded volumes:
    for j in range(bmIdx + 1, dim):
      res[bmIdx, j:dim] = 0.0
      res[j:dim, bmIdx] = 1000.0

  if smoothIt:
    DG.DoTriangleSmoothing(res)
  return res


def UpdatePharmacophoreBounds(bm, atomMatch, pcophore, useDirs=False, dirLength=defaultFeatLength,
                              mol=None):
  """ loops over a distance bounds matrix and replaces the elements
  that are altered by a pharmacophore

  **NOTE** this returns the resulting bounds matrix, but it may also
  alter the input matrix

  atomMatch is a sequence of sequences containing atom indices
  for each of the pharmacophore's features.

    >>> from rdkit import Geometry
    >>> from rdkit.Chem.Pharm3D import Pharmacophore
    >>> feats = [
    ...   ChemicalFeatures.FreeChemicalFeature('HBondAcceptor', 'HAcceptor1',
    ...                                        Geometry.Point3D(0.0, 0.0, 0.0)),
    ...   ChemicalFeatures.FreeChemicalFeature('HBondDonor', 'HDonor1',
    ...                                        Geometry.Point3D(2.65, 0.0, 0.0)),
    ...   ]
    >>> pcophore = Pharmacophore.Pharmacophore(feats)
    >>> pcophore.setLowerBound(0,1, 1.0)
    >>> pcophore.setUpperBound(0,1, 2.0)

    >>> boundsMat = numpy.array([[0.0, 3.0, 3.0],[2.0, 0.0, 3.0],[2.0, 2.0, 0.0]])
    >>> atomMatch = ((0, ), (1, ))
    >>> bm = UpdatePharmacophoreBounds(boundsMat, atomMatch, pcophore)


     In this case, there are no multi-atom features, so the result matrix
     is the same as the input:

     >>> bm is boundsMat
     True

     this means, of course, that the input boundsMat is altered:

     >>> print(', '.join([f'{x:.3f}' for x in boundsMat[0]]))
     0.000, 2.000, 3.000
     >>> print(', '.join([f'{x:.3f}' for x in boundsMat[1]]))
     1.000, 0.000, 3.000
     >>> print(', '.join([f'{x:.3f}' for x in boundsMat[2]]))
     2.000, 2.000, 0.000

  """
  replaceMap = {}
  for i, matchI in enumerate(atomMatch):
    if len(matchI) > 1:
      bm, replaceMap[i] = ReplaceGroup(matchI, bm, useDirs=useDirs)

  for i, matchI in enumerate(atomMatch):
    mi = replaceMap.get(i, matchI[0])
    for j in range(i + 1, len(atomMatch)):
      mj = replaceMap.get(j, atomMatch[j][0])
      if mi < mj:
        idx0, idx1 = mi, mj
      else:
        idx0, idx1 = mj, mi
      bm[idx0, idx1] = pcophore.getUpperBound(i, j)
      bm[idx1, idx0] = pcophore.getLowerBound(i, j)

  return bm


def EmbedPharmacophore(mol, atomMatch, pcophore, randomSeed=-1, count=10, smoothFirst=True,
                       silent=False, bounds=None, excludedVolumes=None, targetNumber=-1,
                       useDirs=False):
  """ Generates one or more embeddings for a molecule that satisfy a pharmacophore

  atomMatch is a sequence of sequences containing atom indices
  for each of the pharmacophore's features.

    - count: is the maximum number of attempts to make a generating an embedding
    - smoothFirst: toggles triangle smoothing of the molecular bounds matix
    - bounds: if provided, should be the molecular bounds matrix. If this isn't
       provided, the matrix will be generated.
    - targetNumber: if this number is positive, it provides a maximum number
       of embeddings to generate (i.e. we'll have count attempts to generate
       targetNumber embeddings).

  returns: a 3 tuple:
    1) the molecular bounds matrix adjusted for the pharmacophore
    2) a list of embeddings (molecules with a single conformer)
    3) the number of failed attempts at embedding

    >>> from rdkit import Geometry
    >>> from rdkit.Chem.Pharm3D import Pharmacophore
    >>> m = Chem.MolFromSmiles('OCCN')
    >>> feats = [
    ...   ChemicalFeatures.FreeChemicalFeature('HBondAcceptor', 'HAcceptor1',
    ...                                        Geometry.Point3D(0.0, 0.0, 0.0)),
    ...   ChemicalFeatures.FreeChemicalFeature('HBondDonor', 'HDonor1',
    ...                                        Geometry.Point3D(2.65, 0.0, 0.0)),
    ...   ]
    >>> pcophore=Pharmacophore.Pharmacophore(feats)
    >>> pcophore.setLowerBound(0,1, 2.5)
    >>> pcophore.setUpperBound(0,1, 3.5)
    >>> atomMatch = ((0, ), (3, ))

    >>> bm,embeds,nFail = EmbedPharmacophore(m, atomMatch, pcophore, randomSeed=23, silent=1)
    >>> len(embeds)
    10
    >>> nFail
    0

    Set up a case that can't succeed:

    >>> pcophore = Pharmacophore.Pharmacophore(feats)
    >>> pcophore.setLowerBound(0,1, 2.0)
    >>> pcophore.setUpperBound(0,1, 2.1)
    >>> atomMatch = ((0, ), (3, ))

    >>> bm, embeds, nFail = EmbedPharmacophore(m, atomMatch, pcophore, randomSeed=23, silent=1)
    >>> len(embeds)
    0
    >>> nFail
    10

  """
  global _times
  if not hasattr(mol, '_chiralCenters'):
    mol._chiralCenters = Chem.FindMolChiralCenters(mol)

  if bounds is None:
    bounds = MolDG.GetMoleculeBoundsMatrix(mol)
  if smoothFirst:
    DG.DoTriangleSmoothing(bounds)

  # print '------------'
  # print 'initial'
  # for row in bm:
  #  print ' ',' '.join(['% 4.2f'%x for x in row])
  # print '------------'

  bm = UpdatePharmacophoreBounds(bounds.copy(), atomMatch, pcophore, useDirs=useDirs, mol=mol)
  if excludedVolumes:
    bm = AddExcludedVolumes(bm, excludedVolumes, smoothIt=False)

  if not DG.DoTriangleSmoothing(bm):
    raise ValueError("could not smooth bounds matrix")

  # print '------------'
  # print 'post replace and smooth'
  # for row in bm:
  #  print ' ',' '.join(['% 4.2f'%x for x in row])
  # print '------------'

  if targetNumber <= 0:
    targetNumber = count
  nFailed = 0
  res = []
  for i in range(count):
    m2 = Chem.Mol(mol)
    t1 = time.perf_counter()
    try:
      if randomSeed <= 0:
        seed = i * 10 + 1
      else:
        seed = i * 10 + randomSeed
      EmbedMol(m2, bm.copy(), atomMatch, randomSeed=seed, excludedVolumes=excludedVolumes)
    except ValueError:
      if not silent:
        logger.info('Embed failed')
      nFailed += 1
    else:
      t2 = time.perf_counter()
      _times['embed'] = _times.get('embed', 0) + t2 - t1
      keepIt = True
      for idx, stereo in mol._chiralCenters:
        if stereo in ('R', 'S'):
          vol = ComputeChiralVolume(m2, idx)
          if (stereo == 'R' and vol >= 0) or (stereo == 'S' and vol <= 0):
            keepIt = False
            break
      if keepIt:
        res.append(m2)
      else:
        logger.debug('Removed embedding due to chiral constraints.')
      if len(res) == targetNumber:
        break
  return bm, res, nFailed


def isNaN(v):
  """ provides an OS independent way of detecting NaNs
  This is intended to be used with values returned from the C++
  side of things.

  We can't actually test this from Python (which traps
  zero division errors), but it would work something like
  this if we could:

  >>> isNaN(0)
  False

  #>>> isNan(1/0)
  #True

  """
  if v != v and sys.platform == 'win32':
    return True
  elif v == 0 and v == 1 and sys.platform != 'win32':
    return True
  return False


def OptimizeMol(mol, bm, atomMatches=None, excludedVolumes=None, forceConstant=1200.0, maxPasses=5,
                verbose=False):
  """  carries out a UFF optimization for a molecule optionally subject
  to the constraints in a bounds matrix

    - atomMatches, if provided, is a sequence of sequences
    - forceConstant is the force constant of the spring used to enforce
      the constraints

   returns a 2-tuple:
     1) the energy of the initial conformation
     2) the energy post-embedding
   NOTE that these energies include the energies of the constraints

    >>> from rdkit import Geometry
    >>> from rdkit.Chem.Pharm3D import Pharmacophore
    >>> m = Chem.MolFromSmiles('OCCN')
    >>> feats = [
    ...  ChemicalFeatures.FreeChemicalFeature('HBondAcceptor', 'HAcceptor1',
    ...                                       Geometry.Point3D(0.0, 0.0, 0.0)),
    ...  ChemicalFeatures.FreeChemicalFeature('HBondDonor', 'HDonor1',
    ...                                       Geometry.Point3D(2.65, 0.0, 0.0)),
    ...  ]
    >>> pcophore=Pharmacophore.Pharmacophore(feats)
    >>> pcophore.setLowerBound(0,1, 2.5)
    >>> pcophore.setUpperBound(0,1, 2.8)
    >>> atomMatch = ((0, ), (3, ))
    >>> bm, embeds, nFail = EmbedPharmacophore(m, atomMatch, pcophore, randomSeed=23, silent=1)
    >>> len(embeds)
    10
    >>> testM = embeds[0]

    Do the optimization:

    >>> e1, e2 = OptimizeMol(testM,bm,atomMatches=atomMatch)

    Optimizing should have lowered the energy:

    >>> e2 < e1
    True

    Check the constrained distance:

    >>> conf = testM.GetConformer(0)
    >>> p0 = conf.GetAtomPosition(0)
    >>> p3 = conf.GetAtomPosition(3)
    >>> d03 = p0.Distance(p3)
    >>> bool(d03 >= pcophore.getLowerBound(0,1) - 0.01)
    True
    >>> bool(d03 <= pcophore.getUpperBound(0,1) + 0.01)
    True

    If we optimize without the distance constraints (provided via the atomMatches
    argument) we're not guaranteed to get the same results, particularly in a case
    like the current one where the pharmacophore brings the atoms uncomfortably
    close together:

    >>> testM = embeds[1]
    >>> e1, e2 = OptimizeMol(testM,bm)
    >>> e2 < e1
    True
    >>> conf = testM.GetConformer(0)
    >>> p0 = conf.GetAtomPosition(0)
    >>> p3 = conf.GetAtomPosition(3)
    >>> d03 = p0.Distance(p3)
    >>> bool(d03 >= pcophore.getLowerBound(0, 1) - 0.01)
    True
    >>> bool(d03 <= pcophore.getUpperBound(0, 1) + 0.01)
    False

  """
  try:
    ff = ChemicalForceFields.UFFGetMoleculeForceField(mol)
  except Exception:
    logger.info('Problems building molecular forcefield', exc_info=True)
    return -1.0, -1.0

  weights = []
  if atomMatches:
    weights = [(i, j) for k in range(len(atomMatches)) for i in atomMatches[k]
               for l in range(k + 1, len(atomMatches)) for j in atomMatches[l]]

  for i, j in weights:
    if j < i:
      i, j = j, i
    ff.AddDistanceConstraint(i, j, bm[j, i], bm[i, j], forceConstant)

  if excludedVolumes:
    nAts = mol.GetNumAtoms()
    conf = mol.GetConformer()
    idx = nAts
    for exVol in excludedVolumes:
      assert exVol.pos is not None
      logger.debug(f'ff.AddExtraPoint({exVol.pos[0]:.4f},{exVol.pos[1]:.4f},{exVol.pos[2]:.4f}')
      ff.AddExtraPoint(exVol.pos[0], exVol.pos[1], exVol.pos[2], True)

      indices = []
      for localIndices, _, _ in exVol.featInfo:
        indices.extend(list(localIndices))
      indicesSet = set(indices)
      del indices

      for i in range(nAts):
        v = numpy.array(conf.GetAtomPosition(i)) - numpy.array(exVol.pos)
        d = numpy.sqrt(numpy.dot(v, v))
        if i not in indicesSet:
          if d < 5.0:
            logger.debug(
              f'ff.AddDistanceConstraint({i},{idx},{exVol.exclusionDist:.3f},1000,{forceConstant:.0f})'
            )
            ff.AddDistanceConstraint(i, idx, exVol.exclusionDist, 1000, forceConstant)

        else:
          logger.debug(f'ff.AddDistanceConstraint({i},{idx},{bm[exVol.index, i]:.3f},'
                       f'{bm[i, exVol.index]:.3f},{forceConstant:.0f})')
          ff.AddDistanceConstraint(i, idx, bm[exVol.index, i], bm[i, exVol.index], forceConstant)
      idx += 1

  ff.Initialize()
  e1 = ff.CalcEnergy()
  if isNaN(e1):
    raise ValueError('bogus energy')
  if verbose:
    print(Chem.MolToMolBlock(mol))
    for i, _ in enumerate(excludedVolumes):
      pos = ff.GetExtraPointPos(i)
      print('   % 7.4f   % 7.4f   % 7.4f As  0  0  0  0  0  0  0  0  0  0  0  0' % tuple(pos),
            file=sys.stderr)

  needsMore = ff.Minimize()
  nPasses = 0
  while needsMore and nPasses < maxPasses:
    needsMore = ff.Minimize()
    nPasses += 1

  e2 = ff.CalcEnergy()
  if isNaN(e2):
    raise ValueError('bogus energy')
  if verbose:
    print('--------')
    print(Chem.MolToMolBlock(mol))
    for i, _ in enumerate(excludedVolumes):
      pos = ff.GetExtraPointPos(i)
      print('   % 7.4f   % 7.4f   % 7.4f Sb  0  0  0  0  0  0  0  0  0  0  0  0' % tuple(pos),
            file=sys.stderr)
  ff = None
  return e1, e2


def EmbedOne(mol, name, match, pcophore, count=1, silent=0, **kwargs):
  """ generates statistics for a molecule's embeddings

  Four energies are computed for each embedding:
      1) E1: the energy (with constraints) of the initial embedding
      2) E2: the energy (with constraints) of the optimized embedding
      3) E3: the energy (no constraints) the geometry for E2
      4) E4: the energy (no constraints) of the optimized free-molecule
         (starting from the E3 geometry)

  Returns a 9-tuple:
      1) the mean value of E1
      2) the sample standard deviation of E1
      3) the mean value of E2
      4) the sample standard deviation of E2
      5) the mean value of E3
      6) the sample standard deviation of E3
      7) the mean value of E4
      8) the sample standard deviation of E4
      9) The number of embeddings that failed

  """
  global _times
  atomMatch = [list(x.GetAtomIds()) for x in match]
  bm, ms, nFailed = EmbedPharmacophore(mol, atomMatch, pcophore, count=count, silent=silent,
                                       **kwargs)
  e1s = []
  e2s = []
  e3s = []
  e4s = []
  d12s = []
  d23s = []
  d34s = []
  for m in ms:
    t1 = time.perf_counter()
    try:
      e1, e2 = OptimizeMol(m, bm, atomMatch)
    except ValueError:
      pass
    else:
      t2 = time.perf_counter()
      _times['opt1'] = _times.get('opt1', 0) + t2 - t1

      e1s.append(e1)
      e2s.append(e2)

      d12s.append(e1 - e2)

      t1 = time.perf_counter()
      try:
        e3, e4 = OptimizeMol(m, bm)
      except ValueError:
        pass
      else:
        t2 = time.perf_counter()
        _times['opt2'] = _times.get('opt2', 0) + t2 - t1
        e3s.append(e3)
        e4s.append(e4)
        d23s.append(e2 - e3)
        d34s.append(e3 - e4)
    count += 1
  try:
    e1, e1d = Stats.MeanAndDev(e1s)
  except Exception:
    e1 = -1.0
    e1d = -1.0
  try:
    e2, e2d = Stats.MeanAndDev(e2s)
  except Exception:
    e2 = -1.0
    e2d = -1.0
  try:
    e3, e3d = Stats.MeanAndDev(e3s)
  except Exception:
    e3 = -1.0
    e3d = -1.0

  try:
    e4, e4d = Stats.MeanAndDev(e4s)
  except Exception:
    e4 = -1.0
    e4d = -1.0
  if not silent:
    print(f'{name}({nFailed}): {e1:.2f}({e1d:.2f}) -> {e2:.2f}({e2d:.2f}) : '
          f'{e3:.2f}({e3d:.2f}) -> {e4:.2f}({e4d:.2f})')
  return e1, e1d, e2, e2d, e3, e3d, e4, e4d, nFailed


def MatchPharmacophoreToMol(mol, featFactory, pcophore):
  """ generates a list of all possible mappings of a pharmacophore to a molecule

  Returns a 2-tuple:
    1) a boolean indicating whether or not all features were found
    2) a list, numFeatures long, of sequences of features


    >>> import os.path
    >>> from rdkit import Geometry, RDConfig
    >>> from rdkit.Chem.Pharm3D import Pharmacophore
    >>> fdefFile = os.path.join(RDConfig.RDCodeDir,'Chem/Pharm3D/test_data/BaseFeatures.fdef')
    >>> featFactory = ChemicalFeatures.BuildFeatureFactory(fdefFile)
    >>> activeFeats = [
    ...  ChemicalFeatures.FreeChemicalFeature('Acceptor', Geometry.Point3D(0.0, 0.0, 0.0)),
    ...  ChemicalFeatures.FreeChemicalFeature('Donor',Geometry.Point3D(0.0, 0.0, 0.0))]
    >>> pcophore= Pharmacophore.Pharmacophore(activeFeats)
    >>> m = Chem.MolFromSmiles('FCCN')
    >>> match, mList = MatchPharmacophoreToMol(m,featFactory,pcophore)
    >>> match
    True

    Two feature types:

    >>> len(mList)
    2

    The first feature type, Acceptor, has two matches:

    >>> len(mList[0])
    2
    >>> mList[0][0].GetAtomIds()
    (0,)
    >>> mList[0][1].GetAtomIds()
    (3,)

    The first feature type, Donor, has a single match:

    >>> len(mList[1])
    1
    >>> mList[1][0].GetAtomIds()
    (3,)

  """
  return MatchFeatsToMol(mol, featFactory, pcophore.getFeatures())


def _getFeatDict(mol, featFactory, features):
  """ **INTERNAL USE ONLY**

    >>> import os.path
    >>> from rdkit import Geometry, RDConfig, Chem
    >>> fdefFile = os.path.join(RDConfig.RDCodeDir, 'Chem/Pharm3D/test_data/BaseFeatures.fdef')
    >>> featFactory = ChemicalFeatures.BuildFeatureFactory(fdefFile)
    >>> activeFeats = [
    ...  ChemicalFeatures.FreeChemicalFeature('Acceptor', Geometry.Point3D(0.0, 0.0, 0.0)),
    ...  ChemicalFeatures.FreeChemicalFeature('Donor', Geometry.Point3D(0.0, 0.0, 0.0))]
    >>> m = Chem.MolFromSmiles('FCCN')
    >>> d = _getFeatDict(m, featFactory, activeFeats)
    >>> sorted(list(d.keys()))
    ['Acceptor', 'Donor']
    >>> donors = d['Donor']
    >>> len(donors)
    1
    >>> donors[0].GetAtomIds()
    (3,)
    >>> acceptors = d['Acceptor']
    >>> len(acceptors)
    2
    >>> acceptors[0].GetAtomIds()
    (0,)
    >>> acceptors[1].GetAtomIds()
    (3,)

  """
  molFeats = {}
  for feat in features:
    family = feat.GetFamily()
    if family not in molFeats:
      molFeats[family] = featFactory.GetFeaturesForMol(mol, includeOnly=family)
  return molFeats


def MatchFeatsToMol(mol, featFactory, features):
  """ generates a list of all possible mappings of each feature to a molecule

  Returns a 2-tuple:
    1) a boolean indicating whether or not all features were found
    2) a list, numFeatures long, of sequences of features


    >>> import os.path
    >>> from rdkit import RDConfig, Geometry
    >>> fdefFile = os.path.join(RDConfig.RDCodeDir, 'Chem/Pharm3D/test_data/BaseFeatures.fdef')
    >>> featFactory = ChemicalFeatures.BuildFeatureFactory(fdefFile)
    >>> activeFeats = [
    ...  ChemicalFeatures.FreeChemicalFeature('Acceptor', Geometry.Point3D(0.0, 0.0, 0.0)),
    ...  ChemicalFeatures.FreeChemicalFeature('Donor', Geometry.Point3D(0.0, 0.0, 0.0))]
    >>> m = Chem.MolFromSmiles('FCCN')
    >>> match, mList = MatchFeatsToMol(m, featFactory, activeFeats)
    >>> match
    True

    Two feature types:

    >>> len(mList)
    2

    The first feature type, Acceptor, has two matches:

    >>> len(mList[0])
    2
    >>> mList[0][0].GetAtomIds()
    (0,)
    >>> mList[0][1].GetAtomIds()
    (3,)

    The first feature type, Donor, has a single match:

    >>> len(mList[1])
    1
    >>> mList[1][0].GetAtomIds()
    (3,)

  """
  molFeats = _getFeatDict(mol, featFactory, features)
  res = []
  for feat in features:
    matches = molFeats.get(feat.GetFamily(), [])
    if len(matches) == 0:
      return False, None
    res.append(matches)
  return True, res


def CombiEnum(sequence):
  """ This generator takes a sequence of sequences as an argument and
  provides all combinations of the elements of the subsequences:

  >>> gen = CombiEnum(((1, 2), (10, 20)))
  >>> next(gen)
  [1, 10]
  >>> next(gen)
  [1, 20]

  >>> [x for x in CombiEnum(((1, 2), (10,20)))]
  [[1, 10], [1, 20], [2, 10], [2, 20]]

  >>> [x for x in CombiEnum(((1, 2),(10, 20), (100, 200)))]
  [[1, 10, 100], [1, 10, 200], [1, 20, 100], [1, 20, 200], [2, 10, 100],
   [2, 10, 200], [2, 20, 100], [2, 20, 200]]

  """
  if not len(sequence):
    yield []
  elif len(sequence) == 1:
    for entry in sequence[0]:
      yield [entry]
  else:
    for entry in sequence[0]:
      for subVal in CombiEnum(sequence[1:]):
        yield [entry] + subVal


def DownsampleBoundsMatrix(bm, indices, maxThresh=4.0):
  """ Removes rows from a bounds matrix that are that are greater
  than a threshold value away from a set of other points

  Returns the modfied bounds matrix

  The goal of this function is to remove rows from the bounds matrix
  that correspond to atoms (atomic index) that are likely to be quite far from
  the pharmacophore we're interested in. Because the bounds smoothing
  we eventually have to do is N^3, this can be a big win

   >>> boundsMat = numpy.array([[0.0, 3.0, 4.0],[2.0, 0.0, 3.0],[2.0, 2.0, 0.0]])
   >>> bm = DownsampleBoundsMatrix(boundsMat,(0, ), 3.5)
   >>> bm.shape == (2, 2)
   True

   we don't touch the input matrix:

   >>> boundsMat.shape == (3, 3)
   True

   >>> print(', '.join([f'{x:.3f}' for x in bm[0]]))
   0.000, 3.000
   >>> print(', '.join([f'{x:.3f}' for x in bm[1]]))
   2.000, 0.000

   if the threshold is high enough, we don't do anything:

   >>> boundsMat = numpy.array([[0.0, 4.0, 3.0],[2.0, 0.0, 3.0],[2.0, 2.0, 0.0]])
   >>> bm = DownsampleBoundsMatrix(boundsMat, (0, ), 5.0)
   >>> bm.shape == (3, 3)
   True

   If there's a max value that's close enough to *any* of the indices
   we pass in, we'll keep it:

   >>> boundsMat = numpy.array([[0.0, 4.0, 3.0],[2.0, 0.0, 3.0],[2.0, 2.0, 0.0]])
   >>> bm = DownsampleBoundsMatrix(boundsMat, (0, 1), 3.5)
   >>> bm.shape == (3, 3)
   True

   However, the datatype should not be changed or uprank into np.float64 as default behaviour
   >>> boundsMat = numpy.array([[0.0, 4.0, 3.0],[2.0, 0.0, 3.0],[2.0, 2.0, 0.0]], dtype=numpy.float32)
   >>> bm = DownsampleBoundsMatrix(boundsMat,(0, 1), 3.5)
   >>> bm.dtype == numpy.float64
   False
   >>> bm.dtype == numpy.float32 or numpy.issubdtype(bm.dtype, numpy.float32)
   True
   >>> bm.dtype == boundsMat.dtype or numpy.issubdtype(bm.dtype, boundsMat.dtype)
   True

  """
  nPts = bm.shape[0]
  if len(indices) == 0:
    return numpy.zeros(shape=tuple([0] * len(bm.shape)), dtype=bm.dtype)
  indicesSet = list(set(indices))
  maskMatrix = numpy.zeros(nPts, dtype=numpy.uint8)
  maskMatrix[indicesSet] = 1
  for idx in indicesSet:
    maskMatrix[numpy.nonzero(bm[idx, idx + 1:] < maxThresh)[0] + (idx + 1)] = 1

  keep = numpy.nonzero(maskMatrix)[0]
  if keep.shape[0] == nPts:
    return bm.copy()
  return bm[numpy.ix_(keep, keep)]


def CoarseScreenPharmacophore(atomMatch, bounds, pcophore, verbose=False):
  """
  >>> from rdkit import Geometry
  >>> from rdkit.Chem.Pharm3D import Pharmacophore
  >>> feats = [
  ...   ChemicalFeatures.FreeChemicalFeature('HBondAcceptor', 'HAcceptor1',
  ...                                        Geometry.Point3D(0.0, 0.0, 0.0)),
  ...   ChemicalFeatures.FreeChemicalFeature('HBondDonor', 'HDonor1',
  ...                                        Geometry.Point3D(2.65, 0.0, 0.0)),
  ...   ChemicalFeatures.FreeChemicalFeature('Aromatic', 'Aromatic1',
  ...                                        Geometry.Point3D(5.12, 0.908, 0.0)),
  ...   ]
  >>> pcophore = Pharmacophore.Pharmacophore(feats)
  >>> pcophore.setLowerBound(0, 1, 1.1)
  >>> pcophore.setUpperBound(0, 1, 1.9)
  >>> pcophore.setLowerBound(0, 2, 2.1)
  >>> pcophore.setUpperBound(0, 2, 2.9)
  >>> pcophore.setLowerBound(1, 2, 2.1)
  >>> pcophore.setUpperBound(1, 2, 3.9)

  >>> bounds = numpy.array([[0, 2, 3],[1, 0, 4],[2, 3, 0]], dtype=numpy.float64)
  >>> CoarseScreenPharmacophore(((0, ),(1, )),bounds, pcophore)
  True

  >>> CoarseScreenPharmacophore(((0, ),(2, )),bounds, pcophore)
  False

  >>> CoarseScreenPharmacophore(((1, ),(2, )),bounds, pcophore)
  False

  >>> CoarseScreenPharmacophore(((0, ),(1, ),(2, )),bounds, pcophore)
  True

  >>> CoarseScreenPharmacophore(((1, ),(0, ),(2, )),bounds, pcophore)
  False

  >>> CoarseScreenPharmacophore(((2, ),(1, ),(0, )),bounds, pcophore)
  False

  # we ignore the point locations here and just use their definitions:

  >>> feats = [
  ...   ChemicalFeatures.FreeChemicalFeature('HBondAcceptor', 'HAcceptor1',
  ...                                        Geometry.Point3D(0.0, 0.0, 0.0)),
  ...   ChemicalFeatures.FreeChemicalFeature('HBondDonor', 'HDonor1',
  ...                                        Geometry.Point3D(2.65, 0.0, 0.0)),
  ...   ChemicalFeatures.FreeChemicalFeature('Aromatic', 'Aromatic1',
  ...                                        Geometry.Point3D(5.12, 0.908, 0.0)),
  ...   ChemicalFeatures.FreeChemicalFeature('HBondDonor', 'HDonor1',
  ...                                        Geometry.Point3D(2.65, 0.0, 0.0)),
  ...                ]
  >>> pcophore=Pharmacophore.Pharmacophore(feats)
  >>> pcophore.setLowerBound(0,1, 2.1)
  >>> pcophore.setUpperBound(0,1, 2.9)
  >>> pcophore.setLowerBound(0,2, 2.1)
  >>> pcophore.setUpperBound(0,2, 2.9)
  >>> pcophore.setLowerBound(0,3, 2.1)
  >>> pcophore.setUpperBound(0,3, 2.9)
  >>> pcophore.setLowerBound(1,2, 1.1)
  >>> pcophore.setUpperBound(1,2, 1.9)
  >>> pcophore.setLowerBound(1,3, 1.1)
  >>> pcophore.setUpperBound(1,3, 1.9)
  >>> pcophore.setLowerBound(2,3, 1.1)
  >>> pcophore.setUpperBound(2,3, 1.9)
  >>> bounds = numpy.array([[0, 3, 3, 3],
  ...                       [2, 0, 2, 2],
  ...                       [2, 1, 0, 2],
  ...                       [2, 1, 1, 0]],
  ...                      dtype=numpy.float64)

  >>> CoarseScreenPharmacophore(((0, ), (1, ), (2, ), (3, )), bounds, pcophore)
  True

  >>> CoarseScreenPharmacophore(((0, ), (1, ), (3, ), (2, )), bounds, pcophore)
  True

  >>> CoarseScreenPharmacophore(((1, ), (0, ), (3, ), (2, )), bounds, pcophore)
  False

  """
  atomMatchSize = len(atomMatch)
  for k in range(atomMatchSize):
    if len(atomMatch[k]) == 1:
      for l in range(k + 1, atomMatchSize):
        if len(atomMatch[l]) == 1:
          if atomMatch[l][0] < atomMatch[k][0]:
            idx0, idx1 = atomMatch[l][0], atomMatch[k][0]
          else:
            idx0, idx1 = atomMatch[k][0], atomMatch[l][0]

          if bounds[idx1, idx0] >= pcophore.getUpperBound(k, l) or \
            bounds[idx0, idx1] <= pcophore.getLowerBound(k, l):
            if verbose:
              print(f'\t  ({idx1},{idx0}) [{k},{l}] fail')
              print(f'\t    {bounds[idx1, idx0]},{pcophore.getUpperBound(k, l)} - '
                    f'{bounds[idx0, idx1]},{pcophore.getLowerBound(k, l)}')
            # logger.debug('\t >%s'%str(atomMatch))
            # logger.debug()
            # logger.debug('\t    %f,%f - %f,%f'%(bounds[idx1,idx0],pcophore.getUpperBound(k,l),
            #                                    bounds[idx0,idx1],pcophore.getLowerBound(k,l)))
            return False
  return True


def Check2DBounds(atomMatch, mol, pcophore):
  """ checks to see if a particular mapping of features onto
  a molecule satisfies a pharmacophore's 2D restrictions

    >>> from rdkit import Geometry
    >>> from rdkit.Chem.Pharm3D import Pharmacophore
    >>> activeFeats = [
    ...  ChemicalFeatures.FreeChemicalFeature('Acceptor', Geometry.Point3D(0.0, 0.0, 0.0)),
    ...  ChemicalFeatures.FreeChemicalFeature('Donor', Geometry.Point3D(0.0, 0.0, 0.0))]
    >>> pcophore= Pharmacophore.Pharmacophore(activeFeats)
    >>> pcophore.setUpperBound2D(0, 1, 3)
    >>> m = Chem.MolFromSmiles('FCC(N)CN')
    >>> Check2DBounds(((0, ), (3, )), m, pcophore)
    True
    >>> Check2DBounds(((0, ), (5, )), m, pcophore)
    False

  """
  dm = Chem.GetDistanceMatrix(mol, False, False, False)
  nFeats = len(atomMatch)
  for i in range(nFeats):
    for j in range(i + 1, nFeats):
      lowerB = pcophore._boundsMat2D[j, i]  # lowerB = pcophore.getLowerBound2D(i,j)
      upperB = pcophore._boundsMat2D[i, j]  # upperB = pcophore.getUpperBound2D(i,j)
      dij = 10000
      for atomI in atomMatch[i]:
        for atomJ in atomMatch[j]:
          try:
            dij = min(dij, dm[atomI, atomJ])
          except IndexError:
            print('bad indices:', atomI, atomJ)
            print('  shape:', dm.shape)
            print('  match:', atomMatch)
            print('    mol:')
            print(Chem.MolToMolBlock(mol))
            raise IndexError
      if dij < lowerB or dij > upperB:
        return False
  return True


def _checkMatch(match, mol, bounds, pcophore, use2DLimits):
  """ **INTERNAL USE ONLY**

  checks whether a particular atom match can be satisfied by
  a molecule

  """
  atomMatch = ChemicalFeatures.GetAtomMatch(match)
  if not atomMatch:
    return None
  elif use2DLimits:
    if not Check2DBounds(atomMatch, mol, pcophore):
      return None
  if not CoarseScreenPharmacophore(atomMatch, bounds, pcophore):
    return None
  return atomMatch


def ConstrainedEnum(matches, mol, pcophore, bounds, use2DLimits=False, index=0, soFar=[]):
  """ Enumerates the list of atom mappings a molecule
  has to a particular pharmacophore.
  We do check distance bounds here.


  """
  nMatches = len(matches)
  if index >= nMatches:
    yield soFar, []
  elif index == nMatches - 1:
    for entry in matches[index]:
      nextStep = soFar + [entry]
      if index != 0:
        atomMatch = _checkMatch(nextStep, mol, bounds, pcophore, use2DLimits)
      else:
        atomMatch = ChemicalFeatures.GetAtomMatch(nextStep)
      if atomMatch:
        yield soFar + [entry], atomMatch
  else:
    for entry in matches[index]:
      nextStep = soFar + [entry]
      if index != 0:
        atomMatch = _checkMatch(nextStep, mol, bounds, pcophore, use2DLimits)
        if not atomMatch:
          continue
      for val in ConstrainedEnum(matches, mol, pcophore, bounds, use2DLimits=use2DLimits,
                                 index=index + 1, soFar=nextStep):
        if val:
          yield val


def MatchPharmacophore(matches, bounds, pcophore, useDownsampling=False, use2DLimits=False,
                       mol=None, excludedVolumes=None, useDirs=False):
  """

  if use2DLimits is set, the molecule must also be provided and topological
  distances will also be used to filter out matches

  """
  for match, atomMatch in ConstrainedEnum(matches, mol, pcophore, bounds, use2DLimits=use2DLimits):
    bm = UpdatePharmacophoreBounds(bounds.copy(), atomMatch, pcophore, useDirs=useDirs, mol=mol)

    if excludedVolumes:
      localEvs = []
      for eV in excludedVolumes:
        featInfo = []
        for i, entry in enumerate(atomMatch):
          info = list(eV.featInfo[i])
          info[0] = entry
          featInfo.append(info)
        localEvs.append(ExcludedVolume.ExcludedVolume(featInfo, eV.index, eV.exclusionDist))
      bm = AddExcludedVolumes(bm, localEvs, smoothIt=False)

    sz = bm.shape[0]
    if useDownsampling:
      indices = []
      for entry in atomMatch:
        indices.extend(entry)
      if excludedVolumes:
        for vol in localEvs:
          indices.append(vol.index)
      bm = DownsampleBoundsMatrix(bm, indices)
    if DG.DoTriangleSmoothing(bm):
      return 0, bm, match, (sz, bm.shape[0])

  return 1, None, None, None


def GetAllPharmacophoreMatches(matches, bounds, pcophore, useDownsampling=0, progressCallback=None,
                               use2DLimits=False, mol=None, verbose=False):
  res = []
  nDone = 0
  for match in CombiEnum(matches):
    atomMatch = ChemicalFeatures.GetAtomMatch(match)
    if atomMatch and use2DLimits and mol:
      pass2D = Check2DBounds(atomMatch, mol, pcophore)
      if verbose:
        print('..', atomMatch)
        print('  ..Pass2d:', pass2D)
    else:
      pass2D = True
    if atomMatch and pass2D and \
       CoarseScreenPharmacophore(atomMatch, bounds, pcophore, verbose=verbose):
      if verbose:
        print('  ..CoarseScreen: Pass')

      if verbose:
        print('pre update:')
        for row in bm:
          print(' ', ' '.join(['% 4.2f' % x for x in row]))

      bm = UpdatePharmacophoreBounds(bounds.copy(), atomMatch, pcophore)
      if verbose:
        print('pre downsample:')
        for row in bm:
          print(' ', ' '.join(['% 4.2f' % x for x in row]))

      if useDownsampling:
        indices = []
        for entry in atomMatch:
          indices.extend(list(entry))
        bm = DownsampleBoundsMatrix(bm, indices)
      if verbose:
        print('post downsample:')
        for row in bm:
          print(' ', ' '.join(['% 4.2f' % x for x in row]))

      if DG.DoTriangleSmoothing(bm):
        res.append(match)
      elif verbose:
        print('cannot smooth')
      nDone += 1
      if progressCallback:
        progressCallback(nDone)
  return res


def ComputeChiralVolume(mol, centerIdx, confId=-1):
  """ Computes the chiral volume of an atom

  We're using the chiral volume formula from Figure 7 of
  Blaney and Dixon, Rev. Comp. Chem. V, 299-335 (1994)

    >>> import os.path
    >>> from rdkit import RDConfig
    >>> dataDir = os.path.join(RDConfig.RDCodeDir,'Chem/Pharm3D/test_data')

    R configuration atoms give negative volumes:

    >>> mol = Chem.MolFromMolFile(os.path.join(dataDir, 'mol-r.mol'))
    >>> Chem.AssignStereochemistry(mol)
    >>> mol.GetAtomWithIdx(1).GetProp('_CIPCode')
    'R'
    >>> ComputeChiralVolume(mol, 1) < 0
    True

    S configuration atoms give positive volumes:

    >>> mol = Chem.MolFromMolFile(os.path.join(dataDir, 'mol-s.mol'))
    >>> Chem.AssignStereochemistry(mol)
    >>> mol.GetAtomWithIdx(1).GetProp('_CIPCode')
    'S'
    >>> ComputeChiralVolume(mol, 1) > 0
    True

    Non-chiral (or non-specified) atoms give zero volume:

    >>> ComputeChiralVolume(mol, 0) == 0.0
    True

    We also work on 3-coordinate atoms (with implicit Hs):

    >>> mol = Chem.MolFromMolFile(os.path.join(dataDir, 'mol-r-3.mol'))
    >>> Chem.AssignStereochemistry(mol)
    >>> mol.GetAtomWithIdx(1).GetProp('_CIPCode')
    'R'
    >>> ComputeChiralVolume(mol, 1) < 0
    True

    >>> mol = Chem.MolFromMolFile(os.path.join(dataDir, 'mol-s-3.mol'))
    >>> Chem.AssignStereochemistry(mol)
    >>> mol.GetAtomWithIdx(1).GetProp('_CIPCode')
    'S'
    >>> ComputeChiralVolume(mol, 1) > 0
    True


  """
  conf = mol.GetConformer(confId)
  Chem.AssignStereochemistry(mol)
  center = mol.GetAtomWithIdx(centerIdx)
  if not center.HasProp('_CIPCode'):
    return 0.0

  nbrs = center.GetNeighbors()
  nbrRanks = [(int(nbr.GetProp('_CIPRank')), conf.GetAtomPosition(nbr.GetIdx())) for nbr in nbrs]

  # if we only have three neighbors (i.e. the determining H isn't present)
  # then use the central atom as the fourth point:
  if len(nbrRanks) == 3:
    nbrRanks.append((-1, conf.GetAtomPosition(centerIdx)))
  nbrRanks.sort()

  ps = [x[1] for x in nbrRanks]
  v1 = ps[0] - ps[3]
  v2 = ps[1] - ps[3]
  v3 = ps[2] - ps[3]
  return v1.DotProduct(v2.CrossProduct(v3))


# ------------------------------------
#
#  doctest boilerplate
#
def _runDoctests(verbose=None):  # pragma: nocover
  import doctest
  failed, _ = doctest.testmod(optionflags=doctest.ELLIPSIS + doctest.NORMALIZE_WHITESPACE,
                              verbose=verbose)
  sys.exit(failed)


if __name__ == '__main__':  # pragma: nocover
  _runDoctests()
