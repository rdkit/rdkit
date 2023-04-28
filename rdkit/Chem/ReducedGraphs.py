# $Id$
#
#  Copyright (c) 2013 Greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
# Created by Greg Landrum, July 2013

import numpy

from rdkit.Chem.rdReducedGraphs import *


def TanimotoSimilarity(arr1, arr2):
  numer = arr1.dot(arr2)
  if numer == 0.0:
    return 0.0
  denom = arr1.dot(arr1) + arr2.dot(arr2) - numer
  if denom == 0.0:
    return 0.0
  return numer / denom
