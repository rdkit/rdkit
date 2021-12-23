# $Id$
#
#  Copyright (C) 2000-2006  greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Matrix operations which may or may not come in handy some day


  **NOTE**: the two functions defined here have been moved to ML.Data.Stats

"""

import sys

from rdkit.ML import files
from rdkit.ML.Data import Stats

FormCovarianceMatrix = Stats.FormCovarianceMatrix
PrincipalComponents = Stats.PrincipalComponents

if __name__ == '__main__':
  fileN = sys.argv[1]
  iV, dV = files.ReadDataFile(fileN)
  eVals, eVects = PrincipalComponents(iV)
  print('eVals: ', eVals)
  print('eVects:', eVects)
