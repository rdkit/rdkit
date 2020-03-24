# -*- coding: utf-8 -*-
"""
MolVS - Molecule Validation and Standardization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MolVS is a python tool built on top of RDKit that performs validation and standardization of chemical structures.

Note that the C++ reimplementation of this is available in the module rdkit.Chem.MolStandardize.rdMolStandardize

:copyright: (c) 2016 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""


import logging

from .standardize import Standardizer, standardize_smiles, enumerate_tautomers_smiles, canonicalize_tautomer_smiles
from .validate import Validator, validate_smiles
from .errors import MolVSError, StandardizeError, ValidateError

from rdkit import Chem
from rdkit.Chem.MolStandardize import tautomer

__title__ = 'MolVS'
__version__ = '0.1.1'
__author__ = 'Matt Swain'
__email__ = 'm.swain@me.com'
__license__ = 'MIT'
__copyright__ = 'Copyright 2016 Matt Swain'


log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())



"""
This module contains tools for reordering the canonical tautomer so generated 
at the first position i.e., among all the tautomers generated the first one 
will be the canonical one after enumerating all the possible tautomers.

"""
def reorderTautomers(self,mol):
    enumerator = self.TautomerEnumerator.enumerate(self,mol)
    canon = enumerator.TautomerCanonicalizer.Canonicalize(self,mol)
    csmi = Chem.MolToSmiles(canon)
    res = [canon]
    tauts = enumerator.enumerate(mol)
    for x in tauts:
        smis = [Chem.MolToSmiles(x)]
    zipped = zip(smis,tauts)
    stpl = sorted((x,y) for x,y in zipped if x!= csmi)
    res += [y for x,y in stpl]
    return res

