#
# Copyright (C) 2001-2022 Greg Landrum and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Python functions for manipulating molecular graphs

In theory much of the functionality in here should be migrating into the
C/C++ codebase.

"""
import types

import numpy

from rdkit import Chem, DataStructs


def CharacteristicPolynomial(mol, mat=None):
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
  I = 1. * numpy.identity(nAtoms)
  An = A
  res = numpy.zeros(nAtoms + 1, float)
  res[0] = 1.0
  for n in range(1, nAtoms + 1):
    res[n] = 1. / n * numpy.trace(An)
    Bn = An - res[n] * I
    An = numpy.dot(A, Bn)

  res[1:] *= -1
  return res
