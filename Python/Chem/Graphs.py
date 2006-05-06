# $Id$
#
# Copyright (C) 2001-2006 greg landrum and rational discovery llc
#
#   @@ All Rights Reserved  @@
#
""" Python functions for manipulating molecular graphs

In theory much of the functionality in here should be migrating into the
C/C++ codebase.

"""
from Numeric import *
import Chem
import DataStructs
import types

def CharacteristicPolynomial(mol,mat=None):
  """ calculates the characteristic polynomial for a molecular graph

    if mat is not passed in, the molecule's Weighted Adjacency Matrix will
    be used.

    The approach used is the Le Verrier-Faddeev-Frame method described
    in _Chemical Graph Theory, 2nd Edition_ by Nenad Trinajstic (CRC Press,
    1992), pg 76.
    
  """
  nAtoms = mol.GetNumAtoms()
  if mat is None:
    # FIX: complete this:
    #A = mol.GetWeightedAdjacencyMatrix()
    pass
  else:
    A = mat
  I = 1.*identity(nAtoms)
  An = A
  res = zeros(nAtoms+1,Float)
  res[0] = 1.0
  for n in xrange(1,nAtoms+1):
    res[n] = 1./n*trace(An)
    Bn = An - res[n]*I
    An = matrixmultiply(A,Bn)

  res[1:] *= -1
  return res  
      

