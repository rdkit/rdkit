"""
Module containing functions for interchange of molecules.
Note that this should be considered beta and that the format
and API will very likely change in future releases.
"""

from collections.abc import Iterable

import rdkit.Chem.rdchem


class JSONParseParameters:
    """Parameters controlling the JSON parser"""

    def __init__(self) -> None: ...

    @property
    def setAromaticBonds(self) -> bool:
        """set bond types to aromatic for bonds flagged aromatic"""

    @setAromaticBonds.setter
    def setAromaticBonds(self, arg: bool, /) -> None: ...

    @property
    def strictValenceCheck(self) -> bool:
        """be strict when checking atom valences"""

    @strictValenceCheck.setter
    def strictValenceCheck(self, arg: bool, /) -> None: ...

    @property
    def parseConformers(self) -> bool:
        """parse conformers in the JSON"""

    @parseConformers.setter
    def parseConformers(self, arg: bool, /) -> None: ...

    @property
    def parseProperties(self) -> bool:
        """parse molecular properties in the JSON"""

    @parseProperties.setter
    def parseProperties(self, arg: bool, /) -> None: ...

    @property
    def useHCounts(self) -> bool:
        """
        use atomic H counts from the JSON. You may want to set
        this to False when parsing queries.
        """

    @useHCounts.setter
    def useHCounts(self, arg: bool, /) -> None: ...

class JSONWriteParameters:
    """Parameters controlling the JSON writer"""

    def __init__(self) -> None: ...

    @property
    def useRDKitExtensions(self) -> bool:
        """use RDKit extensions to the commonchem format"""

    @useRDKitExtensions.setter
    def useRDKitExtensions(self, arg: bool, /) -> None: ...

def MolToJSON(mol: rdkit.Chem.rdchem.Mol, params: JSONWriteParameters = ...) -> str:
    """
    Convert a single molecule to JSON

    ARGUMENTS:
      - mol: the molecule to work with
    RETURNS:
      a string
    """

def MolsToJSON(mols: Iterable[rdkit.Chem.rdchem.Mol], params: JSONWriteParameters = ...) -> str:
    """
    Convert a set of molecules to JSON

    ARGUMENTS:
      - mols: the molecules to work with
    RETURNS:
      a string
    """

def JSONToMols(jsonBlock: str, params: JSONParseParameters = ...) -> tuple:
    """
    Convert JSON to a tuple of molecules

    ARGUMENTS:
      - jsonBlock: the molecule to work with
      - params: (optional) JSONParseParameters controlling the JSON parsing
    RETURNS:
      a tuple of Mols
    """
