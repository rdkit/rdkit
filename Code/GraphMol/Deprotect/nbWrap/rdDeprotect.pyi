"""
Module containing the Deprotect functionality for removing
protecting groups from molecules.
"""

from collections.abc import Sequence

import rdkit.Chem.rdchem


class DeprotectData:
    """
    DeprotectData class, contains a single deprotection reaction and information

     deprotectdata.deprotection_class - functional group being protected
     deprotectdata.reaction_smarts - reaction smarts used for deprotection
     deprotectdata.abbreviation - common abbreviation for the protecting group
     deprotectdata.full_name - full name for the protecting group
    """

    def __init__(self, deprotection_class: str, reaction_smarts: str, abbreviation: str, full_name: str) -> None:
        """
        Construct a new DeprotectData instance.
          >>> reaction_class = "amine"
          >>> reaction_smarts = "[C;R0][C;R0]([C;R0])([O;R0][C;R0](=[O;R0])[NX3;H0,H1:1])C>>[N:1]"
          >>> abbreviation = "Boc"
          >>> full_name = "tert-butyloxycarbonyl"
          >>> data = DeprotectData(reaction_class, reaction_smarts, abbreviation, full_name)
          >>> assert data.isValid()
        """

    @property
    def deprotection_class(self) -> str: ...

    @property
    def full_name(self) -> str: ...

    @property
    def abbreviation(self) -> str: ...

    @property
    def reaction_smarts(self) -> str: ...

    @property
    def example(self) -> str: ...

    def isValid(self) -> bool:
        """Returns True if the DeprotectData has a valid reaction"""

def GetDeprotections() -> list[DeprotectData]:
    """Return the default list of deprotections"""

def Deprotect(mol: rdkit.Chem.rdchem.Mol, deprotections: Sequence[DeprotectData] = ...) -> rdkit.Chem.rdchem.Mol:
    """Return the deprotected version of the molecule."""

def DeprotectInPlace(mol: rdkit.Chem.rdchem.Mol, deprotections: Sequence[DeprotectData] = ...) -> bool:
    """Deprotects the molecule in place."""
