""" Actual implementation of the FMCS algorithm
This code should be used by importing rdkit.Chem.MCS

"""
import warnings

import warnings

warnings.warn("the rdkit.Chem.fmcs module is deprecated, use rdkit.Chem.rdFMCS instead.", DeprecationWarning, stacklevel=2)

from rdkit.Chem.fmcs.fmcs import *
