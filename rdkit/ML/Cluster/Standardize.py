# $Id$
#
# Copyright (C) 2001-2008  greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" contains code for standardization of data matrices for clustering


"""
from rdkit.ML.Data import Stats


def StdDev(mat):
  """ the standard deviation classifier

   This uses _ML.Data.Stats.StandardizeMatrix()_ to do the work

  """
  return Stats.StandardizeMatrix(mat)


methods = [
  ("None", lambda x: x, "No Standardization"),
  ("Standard Deviation", StdDev, "Use the standard deviation"),
]
