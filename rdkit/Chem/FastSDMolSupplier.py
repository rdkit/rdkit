# $Id$
#
#  Copyright (C) 2007 Greg Landrum
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import sys
import warnings

from rdkit import Chem

warnings.warn(
  "The FastSDMolSupplier class has been deprecated, please use Chem.SDMolSupplier instead",
  DeprecationWarning)


class FastSDMolSupplier(Chem.SDMolSupplier):
  pass
