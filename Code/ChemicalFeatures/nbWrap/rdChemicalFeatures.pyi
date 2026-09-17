"""
Module containing free chemical feature functionality
These are features that are not associated with molecules. They are
typically derived from pharmacophores and site-maps.
"""

from typing import overload

import rdkit.Geometry.rdGeometry


class FreeChemicalFeature:
    """
    Class to represent free chemical features.
    These chemical features are not associated with a molecule, though they can be matched
    to molecular features
    """

    @overload
    def __init__(self) -> None:
        """Default Constructor"""

    @overload
    def __init__(self, pickle: str) -> None:
        """Constructor from a pickle string"""

    @overload
    def __init__(self, pickle: bytes) -> None:
        """Constructor from a pickle bytes"""

    @overload
    def __init__(self, family: str, type: str, loc: rdkit.Geometry.rdGeometry.Point3D, id: int = -1) -> None:
        """Constructor with family, type and location specified"""

    @overload
    def __init__(self, family: str, loc: rdkit.Geometry.rdGeometry.Point3D) -> None:
        """constructor with family and location specified, empty type and id"""

    def SetId(self, id: int) -> None:
        """Set the id of the feature"""

    def SetFamily(self, family: str) -> None:
        """Set the family of the feature"""

    def SetType(self, type: str) -> None:
        """Set the specific type for the feature"""

    def GetId(self) -> int:
        """Get the id of the feature"""

    def GetFamily(self) -> str:
        """Get the family of the feature"""

    def GetType(self) -> str:
        """Get the specific type for the feature"""

    def SetPos(self, loc: rdkit.Geometry.rdGeometry.Point3D) -> None:
        """Set the feature position"""

    def GetPos(self) -> rdkit.Geometry.rdGeometry.Point3D:
        """Get the position of the feature"""

    def __getstate__(self) -> tuple[bytes]: ...

    @overload
    def __setstate__(self, arg: tuple[bytes], /) -> None: ...

    @overload
    def __setstate__(self, arg: tuple[str], /) -> None: ...
