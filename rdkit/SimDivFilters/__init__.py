#
#  Copyright (C) 2004  Rational Discovery LLC
#    All Rights Reserved
#
from rdkit import rdBase
try:
  import rdSimDivPickers
  from rdkit.Chem.SimDivFilters.rdSimDivPickers import *
except ImportError:
  rdSimDivPickers=None
