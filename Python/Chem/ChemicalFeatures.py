# $Id$
#
#  Copyright (C) 2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
from rdMolChemicalFeatures import *
from rdChemicalFeatures import *

def MCFF_GetFeaturesForMol(self,mol,includeOnly=""):
  res = []
  count = self.GetNumMolFeatures(mol,includeOnly=includeOnly)
  for i in range(count):
    res.append(self.GetMolFeature(mol,i,includeOnly=includeOnly))
  return tuple(res)
MolChemicalFeatureFactory.GetFeaturesForMol = MCFF_GetFeaturesForMol