#
#  Copyright (C) 2004  Rational Discovery LLC
#    All Rights Reserved
#
from rdkit import rdBase

try:
  from rdkit.SimDivFilters import rdSimDivPickers
  from rdkit.SimDivFilters.rdSimDivPickers import *
except ImportError:
  rdSimDivPickers = None
