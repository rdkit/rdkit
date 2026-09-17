"""Module containing from chemical feature and functions to generate the"""

from collections.abc import Sequence
import os
from typing import overload

import rdkit.Chem.rdchem
import rdkit.Geometry.rdGeometry


class FeatureFileParseException(ValueError):
    pass

def BuildFeatureFactory(fileName: str | os.PathLike) -> MolChemicalFeatureFactory:
    """Construct a feature factory given a feature definition in a file"""

def BuildFeatureFactoryFromString(fdefString: str) -> MolChemicalFeatureFactory:
    """Construct a feature factory given a feature definition block"""

class MolChemicalFeature:
    """
    Class to represent a chemical feature.
    These chemical features may or may not have been derived from molecule object;
    i.e. it is possible to have a chemical feature that was created just from its type
    and location.
    """

    def GetId(self) -> int:
        """Returns the identifier of the feature"""

    def GetFamily(self) -> str:
        """Get the family to which the feature belongs; donor, acceptor, etc."""

    def GetType(self) -> str:
        """Get the specific type for the feature"""

    @overload
    def GetPos(self, confId: int) -> rdkit.Geometry.rdGeometry.Point3D:
        """Get the location of the chemical feature"""

    @overload
    def GetPos(self) -> rdkit.Geometry.rdGeometry.Point3D:
        """Get the location of the default chemical feature (first position)"""

    def GetAtomIds(self) -> tuple:
        """Get the IDs of the atoms that participate in the feature"""

    def GetMol(self) -> rdkit.Chem.rdchem.Mol:
        """Get the molecule used to derive the features"""

    def GetFactory(self) -> MolChemicalFeatureFactory:
        """Get the factory used to generate this feature"""

    def ClearCache(self) -> None:
        """Clears the cache used to store position information."""

    def SetActiveConformer(self, confId: int) -> None:
        """Sets the conformer to use (must be associated with a molecule)."""

    def GetActiveConformer(self) -> int:
        """Gets the conformer to use."""

class MolChemicalFeatureFactory:
    """Class to featurize a molecule"""

    def GetNumFeatureDefs(self) -> int:
        """Get the number of feature definitions"""

    def GetFeatureFamilies(self) -> tuple[str, ...]:
        """Get a tuple of feature types"""

    def GetFeatureDefs(self) -> dict[str, str]:
        """Get a dictionary with SMARTS definitions for each feature type"""

    def GetNumMolFeatures(self, mol: rdkit.Chem.rdchem.Mol, includeOnly: str = '') -> int:
        """Get the number of features the molecule has"""

    def GetMolFeature(self, mol: rdkit.Chem.rdchem.Mol, idx: int, includeOnly: str = '', recompute: bool = True, confId: int = -1) -> MolChemicalFeature:
        """returns a particular feature (by index)"""

def GetAtomMatch(featMatch: Sequence[MolChemicalFeature], maxAts: int = 1024) -> list:
    """
    Returns an empty list if any of the features passed in share an atom.
    Otherwise a list of lists of atom indices is returned.
    """
