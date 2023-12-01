# $Id$
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

import numpy

from rdkit import Chem, Geometry

# BIG NOTE: we are going assume atom IDs starting from 0 instead of 1
# for all the functions in this file. This is so that they
# are reasonably independent of the combicode. However when using
# with combicode the caller needs to make sure the atom IDs from combicode
# are corrected before feeding them in here.


def cross(v1, v2):
  return numpy.array(
    [v1[1] * v2[2] - v1[2] * v2[1], -v1[0] * v2[2] + v1[2] * v2[0], v1[0] * v2[1] - v1[1] * v2[0]],
    dtype=numpy.float64)


def findNeighbors(atomId, adjMat):
  """
  Find the IDs of the neighboring atom IDs
  
  ARGUMENTS:
  atomId - atom of interest
  adjMat - adjacency matrix for the compound
  """
  return [i for i, nid in enumerate(adjMat[atomId]) if nid >= 1]


def _findAvgVec(conf, center, nbrs):
  # find the average of the normalized vectors going from the center atoms to the neighbors
  # the average vector is also normalized
  avgVec = 0
  for nbr in nbrs:
    nid = nbr.GetIdx()
    pt = conf.GetAtomPosition(nid)
    pt -= center
    pt.Normalize()
    if avgVec == 0:
      avgVec = pt
    else:
      avgVec += pt

  avgVec.Normalize()
  return avgVec


def GetAromaticFeatVects(conf, featAtoms, featLoc, scale=1.5):
  """
  Compute the direction vector for an aromatic feature
  
  ARGUMENTS:
     conf - a conformer
     featAtoms - list of atom IDs that make up the feature
     featLoc - location of the aromatic feature specified as point3d
     scale - the size of the direction vector
  """
  dirType = 'linear'
  head = featLoc
  ats = [conf.GetAtomPosition(x) for x in featAtoms]

  v1 = ats[0] - head
  v2 = ats[1] - head
  norm1 = v1.CrossProduct(v2)
  norm1.Normalize()
  norm1 *= scale
  #norm2 = norm1
  norm2 = head - norm1
  norm1 += head
  return ((head, norm1), (head, norm2)), dirType


def ArbAxisRotation(theta, ax, pt):
  theta = math.pi * theta / 180
  c = math.cos(theta)
  s = math.sin(theta)
  t = 1 - c
  X = ax.x
  Y = ax.y
  Z = ax.z
  mat = [[t * X * X + c, t * X * Y + s * Z, t * X * Z - s * Y],
         [t * X * Y - s * Z, t * Y * Y + c, t * Y * Z + s * X],
         [t * X * Z + s * Y, t * Y * Z - s * X, t * Z * Z + c]]
  mat = numpy.array(mat, dtype=numpy.float64)

  if isinstance(pt, Geometry.Point3D):
    pt = numpy.array((pt.x, pt.y, pt.z))
    tmp = numpy.dot(mat, pt)
    return Geometry.Point3D(tmp[0], tmp[1], tmp[2])

  if isinstance(pt, list) or isinstance(pt, tuple):
    res = []
    for p in pt:
      tmp = numpy.dot(mat, numpy.array((p.x, p.y, p.z)))
      res.append(Geometry.Point3D(tmp[0], tmp[1], tmp[2]))
    return res

  return None


def GetAcceptor2FeatVects(conf, featAtoms, scale=1.5):
  """
  Get the direction vectors for Acceptor of type 2
  
  This is the acceptor with two adjacent heavy atoms. We will special case a few things here.
  If the acceptor atom is an oxygen we will assume a sp3 hybridization
  the acceptor directions (two of them)
  reflect that configurations. Otherwise the direction vector in plane with the neighboring
  heavy atoms
  
  ARGUMENTS:
      featAtoms - list of atoms that are part of the feature
      scale - length of the direction vector
  """
  assert len(featAtoms) == 1
  aid = featAtoms[0]
  cpt = conf.GetAtomPosition(aid)

  mol = conf.GetOwningMol()
  nbrs = list(mol.GetAtomWithIdx(aid).GetNeighbors())
  hydrogens = []
  tmp = []
  while len(nbrs):
    nbr = nbrs.pop()
    if nbr.GetAtomicNum() == 1:
      hydrogens.append(nbr)
    else:
      tmp.append(nbr)
  nbrs = tmp
  assert len(nbrs) == 2

  bvec = _findAvgVec(conf, cpt, nbrs)
  bvec *= (-1.0 * scale)

  if mol.GetAtomWithIdx(aid).GetAtomicNum() == 8:
    # assume sp3
    # we will create two vectors by rotating bvec by half the tetrahedral angle in either directions
    v1 = conf.GetAtomPosition(nbrs[0].GetIdx())
    v1 -= cpt
    v2 = conf.GetAtomPosition(nbrs[1].GetIdx())
    v2 -= cpt
    rotAxis = v1 - v2
    rotAxis.Normalize()
    bv1 = ArbAxisRotation(54.5, rotAxis, bvec)
    bv1 += cpt
    bv2 = ArbAxisRotation(-54.5, rotAxis, bvec)
    bv2 += cpt
    return (
      (cpt, bv1),
      (cpt, bv2),
    ), 'linear'

  bvec += cpt
  return ((cpt, bvec), ), 'linear'


def _GetTetrahedralFeatVect(conf, aid, scale):
  mol = conf.GetOwningMol()
  cpt = conf.GetAtomPosition(aid)
  nbrs = mol.GetAtomWithIdx(aid).GetNeighbors()
  if not _checkPlanarity(conf, cpt, nbrs, tol=0.1):
    bvec = _findAvgVec(conf, cpt, nbrs)
    bvec *= (-1.0 * scale)
    bvec += cpt
    return ((cpt, bvec), )
  return tuple()


def GetDonor3FeatVects(conf, featAtoms, scale=1.5):
  """
  Get the direction vectors for Donor of type 3

  This is a donor with three heavy atoms as neighbors. We will assume
  a tetrahedral arrangement of these neighbors. So the direction we are seeking
  is the last fourth arm of the sp3 arrangement

  ARGUMENTS:
    featAtoms - list of atoms that are part of the feature
    scale - length of the direction vector
  """
  assert len(featAtoms) == 1
  return _GetTetrahedralFeatVect(conf, featAtoms[0], scale), 'linear'


def GetAcceptor3FeatVects(conf, featAtoms, scale=1.5):
  """
  Get the direction vectors for Donor of type 3

  This is a donor with three heavy atoms as neighbors. We will assume
  a tetrahedral arrangement of these neighbors. So the direction we are seeking
  is the last fourth arm of the sp3 arrangement

  ARGUMENTS:
    featAtoms - list of atoms that are part of the feature
    scale - length of the direction vector
  """
  assert len(featAtoms) == 1
  return _GetTetrahedralFeatVect(conf, featAtoms[0], scale), 'linear'


def _findHydAtoms(nbrs, atomNames):
  return [nid for nid in nbrs if atomNames[nid] == 'H']


def _checkPlanarity(conf, cpt, nbrs, tol=1.0e-3):
  assert len(nbrs) == 3
  v1 = conf.GetAtomPosition(nbrs[0].GetIdx())
  v1 -= cpt
  v2 = conf.GetAtomPosition(nbrs[1].GetIdx())
  v2 -= cpt
  v3 = conf.GetAtomPosition(nbrs[2].GetIdx())
  v3 -= cpt
  normal = v1.CrossProduct(v2)
  dotP = abs(v3.DotProduct(normal))
  return int(dotP <= tol)


def GetDonor2FeatVects(conf, featAtoms, scale=1.5):
  """
  Get the direction vectors for Donor of type 2

  This is a donor with two heavy atoms as neighbors. The atom may are may not have
  hydrogen on it. Here are the situations with the neighbors that will be considered here
    1. two heavy atoms and two hydrogens: we will assume a sp3 arrangement here
    2. two heavy atoms and one hydrogen: this can either be sp2 or sp3
    3. two heavy atoms and no hydrogens
    
  ARGUMENTS:
    featAtoms - list of atoms that are part of the feature
    scale - length of the direction vector
  """
  assert len(featAtoms) == 1
  aid = featAtoms[0]
  mol = conf.GetOwningMol()
  cpt = conf.GetAtomPosition(aid)

  # find the two atoms that are neighbors of this atoms
  nbrs = list(mol.GetAtomWithIdx(aid).GetNeighbors())
  assert len(nbrs) >= 2

  hydrogens = []
  heavy = []
  for nbr in nbrs:
    if nbr.GetAtomicNum() == 1:
      hydrogens.append(nbr)
    else:
      heavy.append(nbr)

  if len(nbrs) == 2:
    # there should be no hydrogens in this case
    assert len(hydrogens) == 0
    # in this case the direction is the opposite of the average vector of the two neighbors
    bvec = _findAvgVec(conf, cpt, heavy)
    bvec *= (-1.0 * scale)
    bvec += cpt
    return ((cpt, bvec), ), 'linear'

  if len(nbrs) == 3:
    assert len(hydrogens) == 1
    # this is a little more tricky we have to check if the hydrogen is in the plane of the
    # two heavy atoms (i.e. sp2 arrangement) or out of plane (sp3 arrangement)

    # one of the directions will be from hydrogen atom to the heavy atom
    hid = hydrogens[0].GetIdx()
    bvec = conf.GetAtomPosition(hid)
    bvec -= cpt
    bvec.Normalize()
    bvec *= scale
    bvec += cpt
    if _checkPlanarity(conf, cpt, nbrs, tol=1.0e-2):
      # only the hydrogen atom direction needs to be used
      return ((cpt, bvec), ), 'linear'

    # we have a non-planar configuration - we will assume sp3 and compute a second direction vector
    ovec = _findAvgVec(conf, cpt, heavy)
    ovec *= (-1.0 * scale)
    ovec += cpt
    return (
      (cpt, bvec),
      (cpt, ovec),
    ), 'linear'

  if len(nbrs) >= 4:
    # in this case we should have two or more hydrogens we will simple use there directions
    res = []
    for hid in hydrogens:
      hid = hid.GetIdx()
      bvec = conf.GetAtomPosition(hid)
      bvec -= cpt
      bvec.Normalize()
      bvec *= scale
      bvec += cpt
      res.append((cpt, bvec))
    return tuple(res), 'linear'

  return None


def GetDonor1FeatVects(conf, featAtoms, scale=1.5):
  """
  Get the direction vectors for Donor of type 1

  This is a donor with one heavy atom. It is not clear where we should we should be putting the
  direction vector for this. It should probably be a cone. In this case we will just use the
  direction vector from the donor atom to the heavy atom
    
  ARGUMENTS:
    
    featAtoms - list of atoms that are part of the feature
    scale - length of the direction vector
  """
  assert len(featAtoms) == 1
  aid = featAtoms[0]
  mol = conf.GetOwningMol()
  nbrs = mol.GetAtomWithIdx(aid).GetNeighbors()

  # find the neighboring heavy atom
  hnbr = -1
  for nbr in nbrs:
    if nbr.GetAtomicNum() != 1:
      hnbr = nbr.GetIdx()
      break

  cpt = conf.GetAtomPosition(aid)
  v1 = conf.GetAtomPosition(hnbr)
  v1 -= cpt
  v1.Normalize()
  v1 *= (-1.0 * scale)
  v1 += cpt
  return ((cpt, v1), ), 'cone'


def GetAcceptor1FeatVects(conf, featAtoms, scale=1.5):
  """
  Get the direction vectors for Acceptor of type 1

  This is a acceptor with one heavy atom neighbor. There are two possibilities we will
  consider here
  1. The bond to the heavy atom is a single bond e.g. CO
     In this case we don't know the exact direction and we just use the inversion of this bond direction
     and mark this direction as a 'cone'
  2. The bond to the heavy atom is a double bond e.g. C=O
     In this case the we have two possible direction except in some special cases e.g. SO2
     where again we will use bond direction
     
  ARGUMENTS:
    featAtoms - list of atoms that are part of the feature
    scale - length of the direction vector
  """
  assert len(featAtoms) == 1
  aid = featAtoms[0]
  mol = conf.GetOwningMol()
  nbrs = mol.GetAtomWithIdx(aid).GetNeighbors()

  cpt = conf.GetAtomPosition(aid)

  # find the adjacent heavy atom
  heavyAt = -1
  for nbr in nbrs:
    if nbr.GetAtomicNum() != 1:
      heavyAt = nbr
      break

  singleBnd = mol.GetBondBetweenAtoms(aid, heavyAt.GetIdx()).GetBondType() > Chem.BondType.SINGLE

  # special scale - if the heavy atom is a sulfur (we should proabably check phosphorous as well)
  sulfur = heavyAt.GetAtomicNum() == 16

  if singleBnd or sulfur:
    v1 = conf.GetAtomPosition(heavyAt.GetIdx())
    v1 -= cpt
    v1.Normalize()
    v1 *= (-1.0 * scale)
    v1 += cpt
    return ((cpt, v1), ), 'cone'

  # ok in this case we will assume that
  # heavy atom is sp2 hybridized and the direction vectors (two of them)
  # are in the same plane, we will find this plane by looking for one
  # of the neighbors of the heavy atom
  hvNbrs = heavyAt.GetNeighbors()
  hvNbr = -1
  for nbr in hvNbrs:
    if nbr.GetIdx() != aid:
      hvNbr = nbr
      break

  pt1 = conf.GetAtomPosition(hvNbr.GetIdx())
  v1 = conf.GetAtomPosition(heavyAt.GetIdx())
  pt1 -= v1
  v1 -= cpt
  rotAxis = v1.CrossProduct(pt1)
  rotAxis.Normalize()
  bv1 = ArbAxisRotation(120, rotAxis, v1)
  bv1.Normalize()
  bv1 *= scale
  bv1 += cpt
  bv2 = ArbAxisRotation(-120, rotAxis, v1)
  bv2.Normalize()
  bv2 *= scale
  bv2 += cpt
  return (
    (cpt, bv1),
    (cpt, bv2),
  ), 'linear'
