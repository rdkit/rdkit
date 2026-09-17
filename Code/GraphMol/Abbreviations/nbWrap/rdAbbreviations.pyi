"""Module containing functions for working with molecular abbreviations"""

from collections.abc import Sequence

import rdkit.Chem.rdchem


class AbbreviationDefinition:
    """Abbreviation Definition"""

    def __init__(self) -> None: ...

    @property
    def label(self) -> str:
        """the label"""

    @label.setter
    def label(self, arg: str, /) -> None: ...

    @property
    def displayLabel(self) -> str:
        """the label in a drawing when the bond comes from the right"""

    @displayLabel.setter
    def displayLabel(self, arg: str, /) -> None: ...

    @property
    def displayLabelW(self) -> str:
        """the label in a drawing when the bond comes from the west"""

    @displayLabelW.setter
    def displayLabelW(self, arg: str, /) -> None: ...

    @property
    def mol(self) -> rdkit.Chem.rdchem.Mol | None:
        """
        the query molecule (should have a dummy as the first atom if includesXBonds is true)
        """

    @mol.setter
    def mol(self, arg: rdkit.Chem.rdchem.Mol, /) -> None: ...

    @property
    def includesXBonds(self) -> bool:
        """
        whether or not the abbreviation definition includes bonds to non-abbreviation atoms
        """

    @includesXBonds.setter
    def includesXBonds(self, arg: bool, /) -> None: ...

def GetDefaultAbbreviations() -> list[AbbreviationDefinition]:
    """returns a list of the default abbreviation definitions"""

def GetDefaultLinkers() -> list[AbbreviationDefinition]:
    """returns a list of the default linker definitions"""

def ParseAbbreviations(text: str, removeExtraDummies: bool = False, allowConnectionToDummies: bool = False) -> list[AbbreviationDefinition]:
    """
    Returns a set of abbreviation definitions from a string.
    Format of the text data: A series of lines, each of which contains:

    * label
    * SMARTS
    * displayLabel
    * displayLabelW

    Where 'label' is the label used for the abbreviation,
    'SMARTS' is the SMARTS definition of the abbreviation,
    'displayLabel' is used in drawings to render the abbreviations and
    'displayLabelW' is the display label if a bond comes in from the right.
    The 'displayLabel' and 'displayLabelW' fields are optional.
    Use dummies in the SMARTS to indicate attachment points. The assumption
    is that the first atom is a dummy (one will be added if this is not
    true) and that the second atom is the surrogate for the rest of
    the group.
    """

def ParseLinkers(text: str) -> list[AbbreviationDefinition]:
    """
    Returns a set of linker definitions from a string. Equivalent to calling ParseAbbreviations(text, True True).
    """

def CondenseMolAbbreviations(mol: rdkit.Chem.rdchem.Mol, abbrevs: Sequence[AbbreviationDefinition], maxCoverage: float = 0.4, sanitize: bool = True) -> rdkit.Chem.rdchem.Mol:
    """
    Finds and replaces abbreviations in a molecule. The result is not sanitized.
    """

def LabelMolAbbreviations(mol: rdkit.Chem.rdchem.Mol, abbrevs: Sequence[AbbreviationDefinition], maxCoverage: float = 0.4) -> rdkit.Chem.rdchem.Mol:
    """
    Finds abbreviations and adds to them to a molecule as "SUP" SubstanceGroups.
    """

def CondenseAbbreviationSubstanceGroups(mol: rdkit.Chem.rdchem.Mol) -> rdkit.Chem.rdchem.Mol:
    """
    Finds and replaces abbreviation (i.e. "SUP") substance groups in a molecule. The result is not sanitized.
    """
