# $Id$
#
#  Copyright (C) 2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
from rdkit.Chem.rdMolChemicalFeatures import *
from rdkit.Chem.rdChemicalFeatures import *

def MCFF_GetFeaturesForMol(self,mol,includeOnly=""):
  res = []
  count = self.GetNumMolFeatures(mol,includeOnly=includeOnly)
  for i in range(count):
    res.append(self.GetMolFeature(mol,i,includeOnly=includeOnly))
  return tuple(res)
MolChemicalFeatureFactory.GetFeaturesForMol = MCFF_GetFeaturesForMol
