"""
Module containing functionality from the Avalon toolkit.

The functions currently exposed are:
  - GetCanonSmiles()   : return the canonical smiles for a molecule
  - GetAvalonFP()      : return the Avalon fingerprint for a molecule as
                         an RDKit ExplicitBitVector
  - GetAvalonCountFP()      : return the Avalon fingerprint for a molecule as
                              an RDKit SparseIntVector
  - Generate2DCoords() : use the Avalon coordinate generator to create
                         a set of 2D coordinates for a molecule
Each function can be called with either an RDKit molecule or some
molecule data as text (e.g. a SMILES or an MDL mol block).

See the individual docstrings for more information.
"""

import enum
from typing import overload

import rdkit.Chem.rdchem
import rdkit.DataStructs.cDataStructs


@overload
def GetCanonSmiles(mol: rdkit.Chem.rdchem.Mol, flags: int = -1) -> str:
    """returns canonical smiles for an RDKit molecule"""

@overload
def GetCanonSmiles(molData: str, isSmiles: bool, flags: int = -1) -> str:
    """
    Returns canonical smiles for some molecule data.
    If the isSmiles argument is true, the data is assumed to be SMILES, otherwise
    MDL mol data is assumed.
    """

@overload
def GetAvalonFP(mol: rdkit.Chem.rdchem.Mol, nBits: int = 512, isQuery: bool = False, resetVect: bool = False, bitFlags: int = 15761407) -> rdkit.DataStructs.cDataStructs.ExplicitBitVect:
    """returns the Avalon fingerprint for an RDKit molecule"""

@overload
def GetAvalonFP(molData: str, isSmiles: bool, nBits: int = 512, isQuery: bool = False, resetVect: bool = False, bitFlags: int = 15761407) -> rdkit.DataStructs.cDataStructs.ExplicitBitVect:
    """
    returns the Avalon fingerprint for some molecule data.
    If the isSmiles argument is true, the data is assumed to be SMILES, otherwise
    MDL mol data is assumed.
    """

@overload
def Generate2DCoords(mol: rdkit.Chem.rdchem.Mol, clearConfs: bool = True) -> int:
    """Generates 2d coordinates for an RDKit molecule"""

@overload
def Generate2DCoords(molData: str, isSmiles: bool) -> str:
    """
    returns an MDL mol block with 2D coordinates for some molecule data.
    If the isSmiles argument is true, the data is assumed to be SMILES, otherwise
    MDL mol data is assumed.
    """

@overload
def GetAvalonFPAsWords(mol: rdkit.Chem.rdchem.Mol, nBits: int = 512, isQuery: bool = False, resetVect: bool = False, bitFlags: int = 15761407) -> list:
    """returns the Avalon fingerprint for an RDKit molecule as a list of ints"""

@overload
def GetAvalonFPAsWords(molData: str, isSmiles: bool, nBits: int = 512, isQuery: bool = False, resetVect: bool = False, bitFlags: int = 15761407) -> list:
    """
    returns the Avalon fingerprint for some molecule data as a list of ints.
    If the isSmiles argument is true, the data is assumed to be SMILES, otherwise
    MDL mol data is assumed.
    """

@overload
def GetAvalonCountFP(mol: rdkit.Chem.rdchem.Mol, nBits: int = 512, isQuery: bool = False, bitFlags: int = 15761407) -> rdkit.DataStructs.cDataStructs.UIntSparseIntVect:
    """returns the Avalon count fingerprint for an RDKit molecule"""

@overload
def GetAvalonCountFP(molData: str, isSmiles: bool, nBits: int = 512, isQuery: bool = False, bitFlags: int = 15761407) -> rdkit.DataStructs.cDataStructs.UIntSparseIntVect:
    """
    returns the Avalon count fingerprint for some molecule data.
    If the isSmiles argument is true, the data is assumed to be SMILES, otherwise
    MDL mol data is assumed.
    """

def InitializeCheckMol(options: str = '') -> int:
    """
    initializes the structure checker.
    The argument should contain option lines separated by embedded newlines.An empty string will be used if the argument is omitted.An non-zero error code is returned in case of failure.
    """

def CloseCheckMolFiles() -> None:
    """close open files used by molecule-checking functions."""

@overload
def CheckMolecule(molstring: str, isSmiles: bool) -> tuple:
    """
    check a molecule passed in as a string.
    If the isSmiles argument is true, the string should represent the SMILES encoding
    of the molecule, otherwise it should be encoded as an MDL molfile.
    The first member of the return tuple contains the bit-encoded corrections made to the molecule.
    If possible, the molecule (corrected when appropriate) is returned as the second member of
    the return tuple. Otherwise, None is returned.
    """

@overload
def CheckMolecule(mol: rdkit.Chem.rdchem.Mol) -> tuple:
    """
    check a molecule passed in as an RDKit molecule.
    The first member of the return tuple contains the bit-encoded corrections made to the molecule.
    If possible, the molecule (corrected when appropriate) is returned as the second member of
    the return tuple. Otherwise, None is returned.
    """

def CheckMoleculeString(molstring: str, isSmiles: bool) -> tuple:
    """
    check a molecule passed in as a string and returns the result as a string.
    If the isSmiles argument is true, the string should represent the SMILES encoding
    of the molecule, otherwise it should be encoded as an MDL molfile.
    The first member of the return tuple contains the bit-encoded corrections made to the molecule.
    If possible, a corrected CTAB for the molecule is returned as the second member of
    the return tuple.
    """

def GetCheckMolLog() -> str:
    """Returns the Struchk log for the last molecules processed."""

avalonSSSBits: int = 32767

avalonSimilarityBits: int = 15761407

class StruChkFlag(enum.IntEnum):
    bad_molecule = 1

    alias_conversion_failed = 2

    transformed = 4

    fragments_found = 8

    either_warning = 16

    stereo_error = 32

    dubious_stereo_removed = 64

    atom_clash = 128

    atom_check_failed = 256

    size_check_failed = 512

    recharged = 1024

    stereo_forced_bad = 2048

    stereo_transformed = 4096

    template_transformed = 8192

class StruChkResult(enum.IntEnum):
    success = 0

    bad_set = 2979

    transformed_set = 29788
