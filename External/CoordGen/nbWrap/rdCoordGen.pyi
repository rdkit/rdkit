"""Module containing interface to the CoordGen library."""

import rdkit.Chem.rdchem
import rdkit.Geometry.rdGeometry


class CoordGenParams:
    """Parameters controlling coordinate generation"""

    def __init__(self) -> None: ...

    def SetCoordMap(self, coordMap: dict[int, rdkit.Geometry.rdGeometry.Point2D]) -> None:
        """expects a dictionary of Point2D objects with template coordinates"""

    def SetTemplateMol(self, templ: rdkit.Chem.rdchem.Mol) -> None:
        """sets a molecule to be used as the template"""

    @property
    def coordgenScaling(self) -> float:
        """scaling factor for a single bond"""

    @coordgenScaling.setter
    def coordgenScaling(self, arg: float, /) -> None: ...

    @property
    def dbg_useConstrained(self) -> bool:
        """for debugging use"""

    @dbg_useConstrained.setter
    def dbg_useConstrained(self, arg: bool, /) -> None: ...

    @property
    def dbg_useFixed(self) -> bool:
        """for debugging use"""

    @dbg_useFixed.setter
    def dbg_useFixed(self, arg: bool, /) -> None: ...

    @property
    def templateFileDir(self) -> str:
        """directory containing the templates.mae file"""

    @templateFileDir.setter
    def templateFileDir(self, arg: str, /) -> None: ...

    @property
    def sketcherBestPrecision(self) -> float:
        """highest quality (and slowest) precision setting"""

    @property
    def sketcherStandardPrecision(self) -> float:
        """
        standard quality precision setting, the default for the coordgen project
        """

    @property
    def sketcherQuickPrecision(self) -> float:
        """faster precision setting"""

    @property
    def sketcherCoarsePrecision(self) -> float:
        """
        "coarse" (fastest) precision setting, produces good-quality
        coordinates most of the time, this is the default setting for the RDKit
        """

    @property
    def minimizerPrecision(self) -> float:
        """controls sketcher precision"""

    @minimizerPrecision.setter
    def minimizerPrecision(self, arg: float, /) -> None: ...

    @property
    def treatNonterminalBondsToMetalAsZOBs(self) -> bool: ...

    @treatNonterminalBondsToMetalAsZOBs.setter
    def treatNonterminalBondsToMetalAsZOBs(self, arg: bool, /) -> None: ...

    def __setattr__(self, name: str, value: object | None) -> None: ...

def SetDefaultTemplateFileDir(dir: str) -> None: ...

def AddCoords(mol: rdkit.Chem.rdchem.Mol, params: CoordGenParams | None = None) -> None:
    """
    Add 2D coordinates.
    ARGUMENTS:
       - mol: molecule to modify
       - params: (optional) parameters controlling the coordinate generation
    """
