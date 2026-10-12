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

      The polynomial is expanded from the eigenvalues of the matrix, with the
      roots multiplied in order of increasing magnitude.

      The Le Verrier-Faddeev-Frame method previously used here (described in
      _Chemical Graph Theory, 2nd Edition_ by Nenad Trinajstic, CRC Press, 1992,
      pg 76) is a sequential recursion: coefficient k is produced after k matrix
      products have grown the intermediates to the magnitude of the largest
      coefficient, and is then obtained as a difference of quantities of that
      size. Beyond roughly 80 atoms float64 has no significant digits left, and
      because each coefficient feeds the next matrix the error amplifies. For a
      120-atom chain the true final coefficient is +1 and the recursion returned
      -1.97e26. Nothing overflows, so the failure was silent.

    """
  nAtoms = mol.GetNumAtoms()
  if mat is None:
    # FIX: complete this:
    #A = mol.GetWeightedAdjacencyMatrix()
    pass
  else:
    A = mat
  A = numpy.asarray(A, dtype=float)
  # A molecular adjacency matrix is symmetric, and eigvalsh is backward stable for that case.
  # The test must be EXACT: eigvalsh reads only one triangle, so a merely near-symmetric matrix
  # would silently get the polynomial of its symmetrised version instead of its own. The general
  # branch is kept because this function accepts an arbitrary matrix.
  if A.shape[0] == A.shape[1] and numpy.array_equal(A, A.T):
    roots = numpy.linalg.eigvalsh(A)
    res = numpy.array([1.0])
  else:
    roots = numpy.linalg.eigvals(A)
    res = numpy.array([1.0 + 0j])
  # Ordering the roots by increasing magnitude keeps every partial product near the scale of
  # the final coefficients, so nothing cancels. Multiplying them in arbitrary order (as
  # numpy.poly does) still reaches 4.8e4 relative error at 160 atoms.
  for r in roots[numpy.argsort(numpy.abs(roots))]:
    res = numpy.convolve(res, numpy.array([1.0, -r]))
  return numpy.real(res)
