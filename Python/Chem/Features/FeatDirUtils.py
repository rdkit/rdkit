# $Id$
#
# Copyright (C) 2004-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
from Numeric import *
import math

# BIG NOTE: we are going assume atom IDs starting from 0 instead of 1
# for all the functions in this file. This is so that they
# are reasonably indepedent of the combicode. However when using
# with combicode the caller needs to make sure the atom IDs from combicode
# are corrected before feeding them in here.

def cross(v1,v2):
  res = array([ v1[1]*v2[2] - v1[2]*v2[1],
                -v1[0]*v2[2] + v1[2]*v2[0],
                v1[0]*v2[1] - v1[1]*v2[0]],Float)
  return res

def findNeighbors(atomId, adjMat):
  """
  Find the IDs of the neighboring atom IDs
  
  ARGUMENTS:
  atomId - atom of interest
  adjMat - adjacency matrix for the compound
  """
  nbrs = []
  for i,nid in enumerate(adjMat[atomId]):
    if nid >= 1 :
      nbrs.append(i)
  return nbrs

def _findAvgVec(coords, center, nbrs) :
  # find the average of the normalized vectors going from the center atoms to the
  # neighbors
  # the average vector is also normalized
  cpt = array(coords[3*center:3*center+3])
  avgVec = 0
  for nid in nbrs:
    pt = array(coords[3*nid:3*nid+3])
    pt -= cpt
    pt /= (sqrt(innerproduct(pt, pt)))
    if (avgVec == 0) :
      avgVec = pt
    else :
      avgVec += pt

  avgVec /= (sqrt(innerproduct(avgVec, avgVec)))
  return avgVec

def GetAromaticFeatVects(coords, featAtoms, featLoc, scale=1.5):
  """
  Compute the direction vector for an aromatic feature
  
  ARGUMENTS:
     coords - an list of atom coordinates the format is [x1, y1, z1, x2, y2, z3 ....]
     featAtoms - list of atom IDs that make up the feature
     featLoc - location of the aromatic feature specified as a numeric array
     scale - the size of the direction vector
  """
  dirType = 'linear'
  head = featLoc 
  ats = [coords[3*(x):3*(x)+3] for x in featAtoms]
  
  p0 = array(ats[0])
  p1 = array(ats[1])
  v1 = p0-head
  v2 = p1-head
  norm1 = cross(v1,v2)
  norm1 = norm1/sqrt(innerproduct(norm1,norm1))
  norm1 *= scale
  norm2 = -norm1
  norm1 += head
  norm2 += head
  return ( (head,norm1),(head,norm2) ), dirType

def ArbAxisRotation(theta,ax,pt):
  theta = math.pi*theta/180
  c = math.cos(theta)
  s = math.sin(theta)
  t = 1-c
  X = ax[0]
  Y = ax[1]
  Z = ax[2]
  mat = [ [t*X*X+c, t*X*Y+s*Z, t*X*Z-s*Y],
          [t*X*Y-s*Z,t*Y*Y+c,t*Y*Z+s*X],
          [t*X*Z+s*Y,t*Y*Z-s*X,t*Z*Z+c] ]
  mat = array(mat)
  return matrixmultiply(mat,pt)

def GetAcceptor2FeatVects(coords, adjMat, atomNames, featAtoms, scale=1.5):
  """
  Get the direction vectors for Acceptor of type 2
  
  This is the acceptor with two adjacent heavy atoms. We will special case a few things here.
  If the acceptor atom is an oxygen we will assume a sp3 hybridization
  the acceptor directions (two of them)
  reflect that configurations. Otherwise the direction vector in plane with the neighboring
  heavy atoms
  
  ARGUMENTS:
      coords - an list of atom coordinates the format is [x1, y1, z1, x2, y2, z3 ....]
      adjMat - adjacency martix of the compound
      featAtoms - list of atoms that are part of the feature
      scale - length of the direction vector
  """
  assert len(featAtoms) == 1
  aid = featAtoms[0]
  cpt = array(coords[3*aid:3*aid+3])
  nbrs = findNeighbors(aid, adjMat)
  assert len(nbrs) == 2
  
  bvec = _findAvgVec(coords, aid, nbrs)
  bvec *= (-1.0*scale)
  #bvec += cpt

  if (atomNames[aid] == 'O') :
    # assume sp3
    # we will create two vectors by rotating bvec by half the tetrahedral angle in either directions
    v1 = array(coords[3*nbrs[0]:3*nbrs[0]+3])
    v1 -= cpt
    v2 = array(coords[3*nbrs[1]:3*nbrs[1]+3])
    v2 -= cpt
    rotAxis = v1 - v2
    rotAxis /= sqrt(innerproduct(rotAxis, rotAxis))
    bv1 = ArbAxisRotation(54.5, rotAxis, bvec)
    bv1 += cpt
    bv2 = ArbAxisRotation(-54.5, rotAxis, bvec)
    bv2 += cpt
    return ((cpt, bv1), (cpt, bv2),), 'linear'
  else :
    bvec += cpt
    return ((cpt, bvec),), 'linear'

def GetDonor3FeatVects(coords, adjMat, featAtoms, scale=1.5):
  """
  Get the direction vectors for Donor of type 3

  This is a donor with three heavy atoms as neighbors. We will assume
  a tetrahedral arrangement of these neighbors. So the direction we are seeking
  is the last fourth arm of the sp3 arrangment

  ARGUMENTS:
    coords - an list of atom coordinates the format is [x1, y1, z1, x2, y2, z3 ....]
    adjMat - adjacency martix of the compound
    featAtoms - list of atoms that are part of the feature
    scale - length of the direction vector
  """
  assert len(featAtoms) == 1
  aid = featAtoms[0]
  cpt = array(coords[3*aid:3*aid+3])
  nbrs = findNeighbors(aid, adjMat)

  bvec = _findAvgVec(coords, aid, nbrs) 
  bvec *= (-1.0*scale)
  bvec += cpt
  return ((cpt, bvec),), 'linear'

def _findHydAtoms(nbrs, atomNames):
  hAtoms = []
  for nid in nbrs:
    if atomNames[nid] == 'H':
      hAtoms.append(nid)
  return hAtoms

def _checkPlanarity(coords, center, nbrs, tol=1.0e-3):
  assert len(nbrs) == 3
  cpt = array(coords[3*center:3*center+3])
  v1 = array(coords[3*nbrs[0]:3*nbrs[0]+3])
  v1 -= cpt
  v2 = array(coords[3*nbrs[1]:3*nbrs[1]+3])
  v2 -= cpt
  v3 = array(coords[3*nbrs[2]:3*nbrs[2]+3])
  v3 -= cpt
  normal = cross(v1,v2)
  dotP = abs(innerproduct(v3, normal))
  
  if (dotP <= tol) :
    return 1
  else :
    return 0
  
def GetDonor2FeatVects(coords, adjMat, atomNames, featAtoms, scale=1.5) :
  """
  Get the direction vectors for Donor of type 2

  This is a donor with two heavy atoms as neighbors. The atom may are may not have
  hydrogen on it. Here are the situations with the neighbors that will be considered here
    1. two heavy atoms and two hydrogens: we will assume a sp3 arrangement here
    2. two heavy atoms and one hydrogen: this can either be sp2 or sp3
    3. two heavy atoms and no hydrogens
    
  ARGUMENTS:
    coords - an list of atom coordinates the format is [x1, y1, z1, x2, y2, z3 ....]
    adjMat - adjacency martix of the compound
    atomNames - element names of the atoms in the compound
    featAtoms - list of atoms that are part of the feature
    scale - length of the direction vector
  """
  
  assert len(featAtoms) == 1
  aid = featAtoms[0]
  cpt = array(coords[3*aid:3*aid+3])
  # find the two atoms that are neighbors of this atoms
  nbrs = findNeighbors(aid, adjMat)
  assert len(nbrs) >= 2

  hydrogens = _findHydAtoms(nbrs, atomNames)
  
  if len(nbrs) == 2:
    # there should be no hydrogens in this case
    assert len(hydrogens) == 0
    # in this case the direction is the opposite of the average vector of the two neighbors
    bvec = _findAvgVec(coords, aid, nbrs)
    bvec *= (-1.0*scale)
    bvec += cpt
    return ((cpt, bvec),), 'linear'

  elif len(nbrs) == 3:
    assert len(hydrogens) == 1
    # this is a little more tricky we have to check if the hydrogen is in the plane of the
    # two heavy atoms (i.e. sp2 arrangement) or out of plane (sp3 arrangement)

    # one of the directions will be from hydrogen atom to the heavy atom
    hid = hydrogens[0]
    bvec = array(coords[3*hid:3*hid+3])
    bvec -= cpt
    bvec /= (sqrt(innerproduct(bvec, bvec)))
    bvec *= scale
    bvec += cpt
    if _checkPlanarity(coords, aid, nbrs):
      # only the hydrogen atom direction needs to be used
      return ((cpt, bvec),), 'linear'
    else :
      # we have a non-planar configuration - we will assume sp3 and compute a second direction vector
      ovec = _findAvgVec(coords, aid, nbrs)
      ovec *= (-1.0*scale)
      ovec += cpt
      return ((cpt, bvec), (cpt, ovec),), 'linear'

  elif len(nbrs) >= 4 :
    # in this case we should have two or more hydrogens we will simple use there directions
    res = []
    for hid in hydrogenss:
      bvec = array(coords[3*(hid):3*hid+3])
      bvec -= cpt
      bvec /= (sqrt(innerproduct(bvec, bvec)))
      bvec *= scale
      bvec += cpt
      res.append((cpt, bvec))
    return tuple(res), 'linear'

def GetDonor1FeatVects(coords, adjMat, atomNames, featAtoms, scale=1.5) :
  """
  Get the direction vectors for Donor of type 1

  This is a donor with one heavy atom. It is not clear where we should we should be putting the
  direction vector for this. It should probably be a cone. In this case we will just use the
  direction vector from the donor atom to the heavy atom
    
  ARGUMENTS:
    coords - an list of atom coordinates the format is [x1, y1, z1, x2, y2, z3 ....]
    adjMat - adjacency martix of the compound
    atomNames - element names of the atoms in the compound
    featAtoms - list of atoms that are part of the feature
    scale - length of the direction vector
  """
  assert len(featAtoms) == 1
  aid = featAtoms[0]
  nbrs = findNeighbors(aid, adjMat)

  # find the neighboring heavy atom
  hnbr = -1
  for nid in nbrs:
    if atomNames[nid] != 'H':
      hnbr = nid
      break

  cpt = array(coords[3*aid:3*aid+3])
  v1 = array(coords[3*hnbr:3*hnbr+3])
  v1 -= cpt
  v1 /= (sqrt(innerproduct(v1, v1)))
  v1 *= (-1.0*scale)
  v1 += cpt
  return ((cpt, v1),), 'cone'

def GetAcceptor1FeatVects(coords, adjMat, atomNames, featAtoms, scale=1.5) :
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
    coords - an list of atom coordinates the format is [x1, y1, z1, x2, y2, z3 ....]
    adjMat - adjacency martix of the compound
    atomNames - element names of the atoms in the compound
    featAtoms - list of atoms that are part of the feature
    scale - length of the direction vector
  """
  assert len(featAtoms) == 1
  aid = featAtoms[0]
  nbrs = findNeighbors(aid, adjMat)
  cpt = array(coords[3*aid:3*aid+3])
  # find the adjacent heavy atom
  heavyAt = -1
  for nbr in nbrs:
    if atomNames[nbr] != 'H':
      heavyAt = nbr
      break

  singleBnd = 1
  if adjMat[aid][heavyAt] >= 1.5:
    singleBnd = 0

  # special scale - if the heavy atom is a sulfer (we sould proabably check phosphorous as well
  sulfer = 0
  if atomNames[heavyAt] == 'S':
    sulfer = 1
    
  if singleBnd or sulfer:
    v1 = array(coords[3*heavyAt:3*heavyAt+3])
    v1 -= cpt
    v1 /= (sqrt(innerproduct(v1, v1)))
    v1 *= (-1.0*scale)
    v1 += cpt
    return ((cpt, v1),), 'cone'

  else :
    # ok in this case we will assume that
    # heavy atom is sp2 hybridized and the diretion vectors (two of them)
    # are in the same plane, we will find this plane by looking for one
    # of the neighbors of the heavy atom
    hvNbrs = findNeighbors(heavyAt, adjMat)
    hvNbr = -1
    for nbr in hvNbrs:
      if (nbr != aid) :
        hvNbr = nbr
        break
      
    pt1 = array(coords[3*hvNbr:3*hvNbr+3])
    v1 = array(coords[3*heavyAt:3*heavyAt+3])
    pt1 -= v1
    v1 -= cpt
    rotAxis = cross(v1, pt1)
    rotAxis /= sqrt(innerproduct(rotAxis, rotAxis))
    bv1 = ArbAxisRotation(120, rotAxis, v1)
    bv1 /= sqrt(innerproduct(bv1, bv1))
    bv1 *= scale
    bv1 += cpt
    bv2 = ArbAxisRotation(-120, rotAxis, v1)
    bv2 /= sqrt(innerproduct(bv2, bv2))
    bv2 *= scale
    bv2 += cpt
    return ((cpt, bv1), (cpt, bv2),), 'linear'
  
