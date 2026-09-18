"""Module containing a C++ implementation of 2D molecule drawing"""

from collections.abc import Iterable, Sequence
import enum
from typing import Annotated, overload

import numpy
from numpy.typing import NDArray

import rdkit.Chem.rdChemReactions
import rdkit.Chem.rdchem
import rdkit.Geometry.rdGeometry


class MultiColourHighlightStyle(enum.Enum):
    CircleAndLine = 0

    Lasso = 1

CircleAndLine: MultiColourHighlightStyle = MultiColourHighlightStyle.CircleAndLine

Lasso: MultiColourHighlightStyle = MultiColourHighlightStyle.Lasso

class LegendPosition(enum.Enum):
    Bottom = 0

    Top = 1

    Left = 2

    Right = 3

Bottom: LegendPosition = LegendPosition.Bottom

Top: LegendPosition = LegendPosition.Top

Left: LegendPosition = LegendPosition.Left

Right: LegendPosition = LegendPosition.Right

class DrawElement(enum.IntEnum):
    NONE = 0

    PRESHAPES = 1

    BONDS = 2

    ATOMLABELS = 4

    HIGHLIGHTS = 8

    ANNOTATIONS = 16

    RADICALS = 32

    POSTSHAPES = 64

    ALL = 2147483647

NONE: DrawElement = DrawElement.NONE

PRESHAPES: DrawElement = DrawElement.PRESHAPES

BONDS: DrawElement = DrawElement.BONDS

ATOMLABELS: DrawElement = DrawElement.ATOMLABELS

HIGHLIGHTS: DrawElement = DrawElement.HIGHLIGHTS

ANNOTATIONS: DrawElement = DrawElement.ANNOTATIONS

RADICALS: DrawElement = DrawElement.RADICALS

POSTSHAPES: DrawElement = DrawElement.POSTSHAPES

ALL: DrawElement = DrawElement.ALL

class IntStringMap:
    def __getitem__(self, arg: int, /) -> str: ...

    def __setitem__(self, arg0: int, arg1: str, /) -> None: ...

    def __delitem__(self, arg: int, /) -> None: ...

    def __len__(self) -> int: ...

    def __contains__(self, arg: int, /) -> bool: ...

    def __repr__(self) -> str: ...

class MolDrawOptions:
    """Drawing options"""

    def __init__(self) -> None: ...

    @property
    def dummiesAreAttachments(self) -> bool: ...

    @dummiesAreAttachments.setter
    def dummiesAreAttachments(self, arg: bool, /) -> None: ...

    @property
    def circleAtoms(self) -> bool: ...

    @circleAtoms.setter
    def circleAtoms(self, arg: bool, /) -> None: ...

    @property
    def splitBonds(self) -> bool: ...

    @splitBonds.setter
    def splitBonds(self, arg: bool, /) -> None: ...

    @property
    def backgroundColour(self) -> tuple[float, float, float, float]:
        """
        the background colour as an (R,G,B,A) tuple, values should be between 0 and 1
        """

    @backgroundColour.setter
    def backgroundColour(self, arg: tuple[float, ...], /) -> None: ...

    @property
    def queryColour(self) -> tuple[float, float, float, float]:
        """
        the query colour as an (R,G,B,A) tuple, values should be between 0 and 1
        """

    @queryColour.setter
    def queryColour(self, arg: tuple[float, ...], /) -> None: ...

    @property
    def highlightColour(self) -> tuple[float, float, float, float]:
        """
        the highlight colour as an (R,G,B,A) tuple, values should be between 0 and 1
        """

    @highlightColour.setter
    def highlightColour(self, arg: tuple[float, ...], /) -> None: ...

    @property
    def symbolColour(self) -> tuple[float, float, float, float]:
        """
        the symbol colour as an (R,G,B,A) tuple, values should be between 0 and 1
        """

    @symbolColour.setter
    def symbolColour(self, arg: tuple[float, ...], /) -> None: ...

    @property
    def annotationColour(self) -> tuple[float, float, float, float]:
        """
        the annotation colour as an (R,G,B,A) tuple, values should be between 0 and 1
        """

    @annotationColour.setter
    def annotationColour(self, arg: tuple[float, ...], /) -> None: ...

    @property
    def atomNoteColour(self) -> tuple[float, float, float, float]:
        """
        the atom note colour as an (R,G,B,A) tuple, values should be between 0 and 1
        """

    @atomNoteColour.setter
    def atomNoteColour(self, arg: tuple[float, ...], /) -> None: ...

    @property
    def bondNoteColour(self) -> tuple[float, float, float, float]:
        """
        the bond note colour as an (R,G,B,A) tuple, values should be between 0 and 1
        """

    @bondNoteColour.setter
    def bondNoteColour(self, arg: tuple[float, ...], /) -> None: ...

    @property
    def legendColour(self) -> tuple[float, float, float, float]:
        """
        the legend colour as an (R,G,B,A) tuple, values should be between 0 and 1
        """

    @legendColour.setter
    def legendColour(self, arg: tuple[float, ...], /) -> None: ...

    @property
    def variableAttachmentColour(self) -> tuple[float, float, float, float]:
        """
        the variable attachment colour as an (R,G,B,A) tuple, values should be between 0 and 1
        """

    @variableAttachmentColour.setter
    def variableAttachmentColour(self, arg: tuple[float, ...], /) -> None: ...

    def getBackgroundColour(self) -> tuple[float, float, float, float]:
        """method returning the background colour"""

    def getQueryColour(self) -> tuple[float, float, float, float]:
        """method returning the query colour"""

    def getHighlightColour(self) -> tuple[float, float, float, float]:
        """method returning the highlight colour"""

    def setBackgroundColour(self, tpl: tuple[float, ...]) -> None:
        """method for setting the background colour"""

    def setQueryColour(self, tpl: tuple[float, ...]) -> None:
        """method for setting the query colour"""

    def setHighlightColour(self, tpl: tuple[float, ...]) -> None:
        """method for setting the highlight colour"""

    def getSymbolColour(self) -> tuple[float, float, float, float]:
        """method returning the symbol colour"""

    def setSymbolColour(self, tpl: tuple[float, ...]) -> None:
        """method for setting the symbol colour"""

    def getAnnotationColour(self) -> tuple[float, float, float, float]:
        """method returning the annotation colour"""

    def setAnnotationColour(self, tpl: tuple[float, ...]) -> None:
        """method for setting the annotation colour"""

    def setAtomNoteColour(self, tpl: tuple[float, ...]) -> None:
        """method for setting the atom note colour"""

    def getAtomNoteColour(self) -> tuple[float, float, float, float]:
        """method returning the atom note colour"""

    def setBondNoteColour(self, tpl: tuple[float, ...]) -> None:
        """method for setting the bond note colour"""

    def getBondNoteColour(self) -> tuple[float, float, float, float]:
        """method returning the bond note colour"""

    def getLegendColour(self) -> tuple[float, float, float, float]:
        """method returning the legend colour"""

    def setLegendColour(self, tpl: tuple[float, ...]) -> None:
        """method for setting the legend colour"""

    def useDefaultAtomPalette(self) -> None:
        """use the default colour palette for atoms and bonds"""

    def useBWAtomPalette(self) -> None:
        """use a black and white palette for atoms and bonds"""

    def useAvalonAtomPalette(self) -> None:
        """use the Avalon renderer palette for atoms and bonds"""

    def useCDKAtomPalette(self) -> None:
        """use the CDK palette for atoms and bonds"""

    def updateAtomPalette(self, cmap: dict[int, tuple[float, ...]]) -> None:
        """
        updates the palette for atoms and bonds from a dictionary mapping ints to 3-tuples
        """

    def setAtomPalette(self, cmap: dict[int, tuple[float, ...]]) -> None:
        """
        sets the palette for atoms and bonds from a dictionary mapping ints to 3-tuples
        """

    def getAtomPalette(self) -> dict[int, tuple[float, float, float, float]]:
        """
        returns the current atom palette as a dictionary mapping ints to 4-tuples
        """

    @property
    def atomLabels(self) -> IntStringMap:
        """maps indices to atom labels"""

    @atomLabels.setter
    def atomLabels(self, arg: dict[int, str], /) -> None: ...

    @property
    def atomLabelDeuteriumTritium(self) -> bool:
        """labels deuterium as D and tritium as T"""

    @atomLabelDeuteriumTritium.setter
    def atomLabelDeuteriumTritium(self, arg: bool, /) -> None: ...

    @property
    def continuousHighlight(self) -> bool: ...

    @continuousHighlight.setter
    def continuousHighlight(self, arg: bool, /) -> None: ...

    @property
    def fillHighlights(self) -> bool: ...

    @fillHighlights.setter
    def fillHighlights(self, arg: bool, /) -> None: ...

    @property
    def highlightRadius(self) -> float:
        """Default radius for highlight circles."""

    @highlightRadius.setter
    def highlightRadius(self, arg: float, /) -> None: ...

    @property
    def flagCloseContactsDist(self) -> int: ...

    @flagCloseContactsDist.setter
    def flagCloseContactsDist(self, arg: int, /) -> None: ...

    @property
    def atomRegions(self) -> list[list[int]]:
        """regions to outline"""

    @atomRegions.setter
    def atomRegions(self, arg: Sequence[Sequence[int]], /) -> None: ...

    @property
    def includeAtomTags(self) -> bool:
        """include atom tags in output"""

    @includeAtomTags.setter
    def includeAtomTags(self, arg: bool, /) -> None: ...

    @property
    def clearBackground(self) -> bool:
        """clear the background before drawing a molecule"""

    @clearBackground.setter
    def clearBackground(self, arg: bool, /) -> None: ...

    @property
    def legendFontSize(self) -> int:
        """font size in pixels of the legend (if drawn)"""

    @legendFontSize.setter
    def legendFontSize(self, arg: int, /) -> None: ...

    @property
    def legendFraction(self) -> float:
        """fraction of the draw panel to be used for the legend if present"""

    @legendFraction.setter
    def legendFraction(self, arg: float, /) -> None: ...

    @property
    def legendPosition(self) -> LegendPosition:
        """
        legend position enum. Default=Bottom.
        Values: LegendPosition.Bottom, LegendPosition.Top, LegendPosition.Left, LegendPosition.Right.
        """

    @legendPosition.setter
    def legendPosition(self, arg: LegendPosition, /) -> None: ...

    @property
    def legendVerticalText(self) -> bool:
        """when legend is Left or Right, draw text vertically (one char per line)"""

    @legendVerticalText.setter
    def legendVerticalText(self, arg: bool, /) -> None: ...

    @property
    def maxFontSize(self) -> int:
        """maximum font size in pixels. default=40, -1 means no maximum."""

    @maxFontSize.setter
    def maxFontSize(self, arg: int, /) -> None: ...

    @property
    def minFontSize(self) -> int:
        """minimum font size in pixels. default=6, -1 means no minimum."""

    @minFontSize.setter
    def minFontSize(self, arg: int, /) -> None: ...

    @property
    def fixedFontSize(self) -> int:
        """
        font size in pixels. default=-1 means not fixed.  If set,
        always used irrespective of scale, minFontSize and maxFontSize.
        """

    @fixedFontSize.setter
    def fixedFontSize(self, arg: int, /) -> None: ...

    @property
    def baseFontSize(self) -> float:
        """relative size of font.  Defaults to 0.6.  -1 means use default."""

    @baseFontSize.setter
    def baseFontSize(self, arg: float, /) -> None: ...

    @property
    def annotationFontScale(self) -> float:
        """
        Scale of font for atom and bond annotation relative to atom
        label font.  Default=0.75.
        """

    @annotationFontScale.setter
    def annotationFontScale(self, arg: float, /) -> None: ...

    @property
    def fontFile(self) -> str:
        """
        Font file for use with FreeType text drawer.  Can also be
        BuiltinTelexRegular (the default) or BuiltinRobotoRegular.
        """

    @fontFile.setter
    def fontFile(self, arg: str, /) -> None: ...

    @property
    def multipleBondOffset(self) -> float:
        """
        offset for the extra lines in a multiple bond as a fraction of mean bond length
        """

    @multipleBondOffset.setter
    def multipleBondOffset(self, arg: float, /) -> None: ...

    @property
    def padding(self) -> float:
        """Fraction of empty space to leave around molecule.  Default=0.05."""

    @padding.setter
    def padding(self, arg: float, /) -> None: ...

    @property
    def reagentPadding(self) -> float:
        """
        Fraction of empty space to leave around each component
        of a reaction drawing.  Default=0.0.
        """

    @reagentPadding.setter
    def reagentPadding(self, arg: float, /) -> None: ...

    @property
    def bondLineWidth(self) -> float:
        """if positive, this overrides the default line width for bonds"""

    @bondLineWidth.setter
    def bondLineWidth(self, arg: float, /) -> None: ...

    @property
    def scaleBondWidth(self) -> bool:
        """Scales the width of drawn bonds using image scaling."""

    @scaleBondWidth.setter
    def scaleBondWidth(self, arg: bool, /) -> None: ...

    @property
    def scaleHighlightBondWidth(self) -> bool:
        """Scales the width of drawn highlighted bonds using image scaling."""

    @scaleHighlightBondWidth.setter
    def scaleHighlightBondWidth(self, arg: bool, /) -> None: ...

    @property
    def highlightBondWidthMultiplier(self) -> int:
        """
        What to multiply default bond width by for highlighting bonds. Default-8.
        """

    @highlightBondWidthMultiplier.setter
    def highlightBondWidthMultiplier(self, arg: int, /) -> None: ...

    @property
    def prepareMolsBeforeDrawing(self) -> bool:
        """call prepareMolForDrawing() on each molecule passed to DrawMolecules()"""

    @prepareMolsBeforeDrawing.setter
    def prepareMolsBeforeDrawing(self, arg: bool, /) -> None: ...

    @property
    def fixedScale(self) -> float:
        """
        If > 0.0, fixes scale to that fraction of width of
        draw window unless that would make it too big.  Default -1.0 means adjust scale to fit.
        """

    @fixedScale.setter
    def fixedScale(self, arg: float, /) -> None: ...

    @property
    def fixedBondLength(self) -> float:
        """
        If > 0.0, fixes bond length to this number of pixels
        unless that would make it too big.  Default -1.0 means
        no fix.  If both set, fixedScale takes precedence.
        """

    @fixedBondLength.setter
    def fixedBondLength(self, arg: float, /) -> None: ...

    @property
    def rotate(self) -> float:
        """Rotates molecule about centre by this number of degrees,"""

    @rotate.setter
    def rotate(self, arg: float, /) -> None: ...

    @property
    def addStereoAnnotation(self) -> bool:
        """adds R/S and E/Z to drawings. Default False."""

    @addStereoAnnotation.setter
    def addStereoAnnotation(self, arg: bool, /) -> None: ...

    @property
    def showAllCIPCodes(self) -> bool:
        """show all defined CIP codes (no hiding!). Default False."""

    @showAllCIPCodes.setter
    def showAllCIPCodes(self, arg: bool, /) -> None: ...

    @property
    def addAtomIndices(self) -> bool:
        """adds atom indices to drawings. Default False."""

    @addAtomIndices.setter
    def addAtomIndices(self, arg: bool, /) -> None: ...

    @property
    def addBondIndices(self) -> bool:
        """adds bond indices to drawings. Default False."""

    @addBondIndices.setter
    def addBondIndices(self, arg: bool, /) -> None: ...

    @property
    def isotopeLabels(self) -> bool:
        """adds isotope labels on non-dummy atoms. Default True."""

    @isotopeLabels.setter
    def isotopeLabels(self, arg: bool, /) -> None: ...

    @property
    def dummyIsotopeLabels(self) -> bool:
        """adds isotope labels on dummy atoms. Default True."""

    @dummyIsotopeLabels.setter
    def dummyIsotopeLabels(self, arg: bool, /) -> None: ...

    @property
    def atomHighlightsAreCircles(self) -> bool:
        """
        forces atom highlights always to be circles.
        Default (false) is to put ellipses round longer labels.
        """

    @atomHighlightsAreCircles.setter
    def atomHighlightsAreCircles(self, arg: bool, /) -> None: ...

    @property
    def multiColourHighlightStyle(self) -> MultiColourHighlightStyle:
        """
        Either 'CircleAndLine' or 'Lasso', to control style of
        multi-coloured highlighting in DrawMoleculeWithHighlights.
        Default is CircleAndLine.
        """

    @multiColourHighlightStyle.setter
    def multiColourHighlightStyle(self, arg: MultiColourHighlightStyle, /) -> None: ...

    @property
    def centreMoleculesBeforeDrawing(self) -> bool:
        """Moves the centre of the drawn molecule to (0,0). Default False."""

    @centreMoleculesBeforeDrawing.setter
    def centreMoleculesBeforeDrawing(self, arg: bool, /) -> None: ...

    @property
    def additionalAtomLabelPadding(self) -> float:
        """
        additional padding to leave around atom labels.
        Expressed as a fraction of the font size.
        """

    @additionalAtomLabelPadding.setter
    def additionalAtomLabelPadding(self, arg: float, /) -> None: ...

    @property
    def noAtomLabels(self) -> bool:
        """disables inclusion of atom labels in the rendering"""

    @noAtomLabels.setter
    def noAtomLabels(self, arg: bool, /) -> None: ...

    @property
    def explicitMethyl(self) -> bool:
        """Draw terminal methyls explictly.  Default is false."""

    @explicitMethyl.setter
    def explicitMethyl(self, arg: bool, /) -> None: ...

    @property
    def includeMetadata(self) -> bool:
        """
        When possible, include metadata about molecules and reactions to
        allow them to be reconstructed. Default is true.
        """

    @includeMetadata.setter
    def includeMetadata(self, arg: bool, /) -> None: ...

    @property
    def includeRadicals(self) -> bool:
        """
        include radicals in the drawing (it can be useful to turn this off
        for reactions and queries). Default is true.
        """

    @includeRadicals.setter
    def includeRadicals(self, arg: bool, /) -> None: ...

    @property
    def comicMode(self) -> bool:
        """
        simulate hand-drawn lines for bonds. When combined with
        a font like Comic-Sans or Comic-Neue, this gives
        xkcd-like drawings. Default is false.
        """

    @comicMode.setter
    def comicMode(self, arg: bool, /) -> None: ...

    @property
    def variableBondWidthMultiplier(self) -> int:
        """
        what to multiply standard bond width by for variable attachment points.
        """

    @variableBondWidthMultiplier.setter
    def variableBondWidthMultiplier(self, arg: int, /) -> None: ...

    @property
    def variableAtomRadius(self) -> float:
        """radius value to use for atoms involved in variable attachment points."""

    @variableAtomRadius.setter
    def variableAtomRadius(self, arg: float, /) -> None: ...

    @property
    def includeChiralFlagLabel(self) -> bool:
        """
        add a molecule annotation with "ABS" if the chiral
        flag is set. Default is false.
        """

    @includeChiralFlagLabel.setter
    def includeChiralFlagLabel(self, arg: bool, /) -> None: ...

    @property
    def simplifiedStereoGroupLabel(self) -> bool:
        """
        if all specified stereocenters are in a single
        StereoGroup, show a molecule-level annotation instead of
        the individual labels. Default is false.
        """

    @simplifiedStereoGroupLabel.setter
    def simplifiedStereoGroupLabel(self, arg: bool, /) -> None: ...

    @property
    def unspecifiedStereoIsUnknown(self) -> bool:
        """
        if true, double bonds with unspecified stereo are drawn
        crossed, potential stereocenters with unspecified stereo
        are drawn with a wavy bond. Default is false.
        """

    @unspecifiedStereoIsUnknown.setter
    def unspecifiedStereoIsUnknown(self, arg: bool, /) -> None: ...

    @property
    def singleColourWedgeBonds(self) -> bool:
        """
        if true wedged and dashed bonds are drawn using symbolColour
        rather than inheriting their colour from the atoms.
        Default is false.
        """

    @singleColourWedgeBonds.setter
    def singleColourWedgeBonds(self, arg: bool, /) -> None: ...

    @property
    def singleColourBonds(self) -> bool:
        """
        if true all bonds are drawn using symbolColour rather than inheriting their colour from the atoms. Default is false.
        """

    @singleColourBonds.setter
    def singleColourBonds(self, arg: bool, /) -> None: ...

    @property
    def useMolBlockWedging(self) -> bool:
        """
        If the molecule came from a MolBlock, prefer the wedging
        information that provides.  If false, use RDKit rules.
        Default false
        """

    @useMolBlockWedging.setter
    def useMolBlockWedging(self, arg: bool, /) -> None: ...

    @property
    def scalingFactor(self) -> float:
        """
        scaling factor for pixels->angstrom when auto scaling
        being used.  Default is 20.
        """

    @scalingFactor.setter
    def scalingFactor(self, arg: float, /) -> None: ...

    @property
    def drawMolsSameScale(self) -> bool:
        """
        when drawing multiple molecules with DrawMolecules,
        forces them to use the same scale.  Default is true.
        """

    @drawMolsSameScale.setter
    def drawMolsSameScale(self, arg: bool, /) -> None: ...

    @property
    def useComplexQueryAtomSymbols(self) -> bool:
        """
        replace any atom, any hetero, any halo queries
        with complex query symbols A, Q, X, M, optionally followed
        by H if hydrogen is included (except for AH, which stays *).
        Default is true
        """

    @useComplexQueryAtomSymbols.setter
    def useComplexQueryAtomSymbols(self, arg: bool, /) -> None: ...

    @property
    def bracketsAroundAtomLists(self) -> bool:
        """
        Whether to put brackets round atom lists in query atoms.
        Default is true.
        """

    @bracketsAroundAtomLists.setter
    def bracketsAroundAtomLists(self, arg: bool, /) -> None: ...

    @property
    def standardColoursForHighlightedAtoms(self) -> bool:
        """
        If true, highlighted hetero atoms are drawn in standard colours
        rather than black.  Default=False
        """

    @standardColoursForHighlightedAtoms.setter
    def standardColoursForHighlightedAtoms(self, arg: bool, /) -> None: ...

    @property
    def drawingExtentsInclude(self) -> int:
        """
        Drawing extents are computed taking into account only selected
        DrawElement items.  Default=DrawElement.ALL
        """

    @drawingExtentsInclude.setter
    def drawingExtentsInclude(self, arg: int, /) -> None: ...

    def getVariableAttachmentColour(self) -> tuple[float, float, float, float]:
        """method for getting the colour of variable attachment points"""

    def setVariableAttachmentColour(self, tpl: tuple[float, ...]) -> None:
        """method for setting the colour of variable attachment points"""

    @property
    def stereoGroupAndLabel(self) -> str:
        """String to use for enhanced stereo 'AND' groups.  Default='and'."""

    @stereoGroupAndLabel.setter
    def stereoGroupAndLabel(self, arg: str, /) -> None: ...

    @property
    def stereoGroupOrLabel(self) -> str:
        """String to use for enhanced stereo 'OR' groups.  Default='or'."""

    @stereoGroupOrLabel.setter
    def stereoGroupOrLabel(self, arg: str, /) -> None: ...

    @property
    def stereoGroupAbsLabel(self) -> str:
        """String to use for enhanced stereo 'ABS' groups.  Default='abs'."""

    @stereoGroupAbsLabel.setter
    def stereoGroupAbsLabel(self, arg: str, /) -> None: ...

    @property
    def addStereoGroupAnnotation(self) -> bool:
        """Whether to add the enhanced stereo labels.  Default is True."""

    @addStereoGroupAnnotation.setter
    def addStereoGroupAnnotation(self, arg: bool, /) -> None: ...

    def __setattr__(self, name: str, value: object | None) -> None: ...

class MolDraw2D:
    """Drawer abstract base class"""

    def SetFontSize(self, new_size: float) -> None:
        """change the default font size. The units are, roughly, pixels."""

    def FontSize(self) -> float:
        """get the default font size. The units are, roughly, pixels."""

    @overload
    def DrawMolecule(self, mol: rdkit.Chem.rdchem.Mol, highlightAtoms: Iterable[int] | None = None, highlightAtomColors: dict[int, tuple[float, ...]] | None = None, highlightAtomRadii: dict[int, float] | None = None, confId: int = -1, legend: str = '') -> None: ...

    @overload
    def DrawMolecule(self, mol: rdkit.Chem.rdchem.Mol, highlightAtoms: Iterable[int] | None, highlightBonds: Iterable[int] | None, highlightAtomColors: dict[int, tuple[float, ...]] | None = None, highlightBondColors: dict[int, tuple[float, ...]] | None = None, highlightAtomRadii: dict[int, float] | None = None, confId: int = -1, legend: str = '') -> None:
        """renders a molecule"""

    def GetMolSize(self, mol: rdkit.Chem.rdchem.Mol, highlightAtoms: Iterable[int] | None = None, highlightBonds: Iterable[int] | None = None, highlightAtomColors: dict[int, tuple[float, ...]] | None = None, highlightBondColors: dict[int, tuple[float, ...]] | None = None, highlightAtomRadii: dict[int, float] | None = None, confId: int = -1, legend: str = '') -> tuple[int, int]:
        """
        returns the width and height required to draw a molecule at the current size
        """

    def DrawMoleculeWithHighlights(self, mol: rdkit.Chem.rdchem.Mol, legend: str, highlight_atom_map: dict[int, list[tuple[float, ...]]] | None, highlight_bond_map: dict[int, list[tuple[float, ...]]] | None, highlight_radii: dict[int, float] | None, highlight_linewidth_multipliers: dict[int, int] | None, confId: int = -1) -> None:
        """renders a molecule with multiple highlight colours"""

    def DrawMolecules(self, mols: Iterable[rdkit.Chem.rdchem.Mol], highlightAtoms: Sequence[Iterable[int]] | None = None, highlightBonds: Sequence[Iterable[int]] | None = None, highlightAtomColors: Sequence[dict[int, tuple[float, ...]]] | None = None, highlightBondColors: Sequence[dict[int, tuple[float, ...]]] | None = None, highlightAtomRadii: Sequence[dict[int, float]] | None = None, confIds: Iterable[int] | None = None, legends: Iterable[str] | None = None) -> None:
        """renders multiple molecules"""

    def DrawReaction(self, rxn: rdkit.Chem.rdChemReactions.ChemicalReaction, highlightByReactant: bool = False, highlightColorsReactants: list[tuple[float, ...]] | None = None, confIds: Iterable[int] | None = None) -> None:
        """renders a reaction"""

    def Width(self) -> int:
        """get the width of the drawing canvas"""

    def Height(self) -> int:
        """get the height of the drawing canvas"""

    def SetOffset(self, x: int, y: int) -> None:
        """set the offset (in drawing coordinates) for the drawing"""

    def Offset(self) -> rdkit.Geometry.rdGeometry.Point2D:
        """returns the offset (in drawing coordinates) for the drawing"""

    def SetScale(self, width: int, height: int, minv: rdkit.Geometry.rdGeometry.Point2D, maxv: rdkit.Geometry.rdGeometry.Point2D, mol: rdkit.Chem.rdchem.Mol | None = None) -> None:
        """uses the values provided to set the drawing scaling"""

    def FlexiMode(self) -> bool:
        """returns whether or not FlexiMode is being used"""

    def SetFlexiMode(self, mode: bool) -> None:
        """
        when FlexiMode is set, molecules will always been drawn with the default values for bond length, font size, etc.
        """

    def SetLineWidth(self, width: float) -> None:
        """set the line width being used"""

    def SetColour(self, tpl: tuple[float, ...]) -> None:
        """set the color being used fr drawing and filling"""

    def LineWidth(self) -> float:
        """returns the line width being used"""

    def SetFillPolys(self, val: bool) -> None:
        """sets whether or not polygons are filled"""

    def FillPolys(self) -> bool:
        """returns whether or not polygons are being filled"""

    def DrawLine(self, cds1: rdkit.Geometry.rdGeometry.Point2D, cds2: rdkit.Geometry.rdGeometry.Point2D, rawCoords: bool = False) -> None:
        """
        draws a line with the current drawing style. The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels.
        """

    def DrawArrow(self, cds1: rdkit.Geometry.rdGeometry.Point2D, cds2: rdkit.Geometry.rdGeometry.Point2D, asPolygon: bool = False, frac: float = 0.05, angle: float = 0.5235987755982988, color: tuple[float, ...] | None = None, rawCoords: bool = False) -> None:
        """
        draws an arrow with the current drawing style. The coordinates
        are in the molecule frame unless rawCoords is true,
        in which case the coordinates are in pixels.
        If asPolygon is true the head of the
        arrow will be drawn as a triangle, otherwise two lines are used.
        The fraction of the arrow length to use for the head is given by
        frac. The angle of the arrowhead
        (the angle between the main line and each arrowhead line) is given by angle.
        The color is a tuple of 3 floats (0-1) in red, green, blue (RGB) order.
        """

    def DrawTriangle(self, cds1: rdkit.Geometry.rdGeometry.Point2D, cds2: rdkit.Geometry.rdGeometry.Point2D, cds3: rdkit.Geometry.rdGeometry.Point2D, rawCoords: bool = False) -> None:
        """
        draws a triangle with the current drawing style. The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels.
        """

    def DrawPolygon(self, cds: Iterable[rdkit.Geometry.rdGeometry.Point2D], rawCoords: bool = False) -> None:
        """
        draws a polygon with the current drawing style. The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels.
        """

    def DrawEllipse(self, cds1: rdkit.Geometry.rdGeometry.Point2D, cds2: rdkit.Geometry.rdGeometry.Point2D, rawCoords: bool = False) -> None:
        """
        draws a triangle with the current drawing style in the rectangle defined by the two points. The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels.
        """

    def DrawRect(self, cds1: rdkit.Geometry.rdGeometry.Point2D, cds2: rdkit.Geometry.rdGeometry.Point2D, rawCoords: bool = False) -> None:
        """
        draws a rectangle with the current drawing style in the rectangle defined by the two points. The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels.
        """

    def DrawArc(self, center: rdkit.Geometry.rdGeometry.Point2D, radius: float, angle1: float, angle2: float, rawCoords: bool = False) -> None:
        """
        draws an arc with the current drawing style. The coordinates
        are in the molecule frame unless rawCoords is true,
        in which case the coordinates are in pixels.
        The angles are in degrees; angle2 should be > angle1.
        """

    def DrawAttachmentLine(self, cds1: rdkit.Geometry.rdGeometry.Point2D, cds2: rdkit.Geometry.rdGeometry.Point2D, color: tuple[float, ...], len: float = 1.0, nSegments: int = 16, rawCoords: bool = False) -> None:
        """
        draw a line indicating the presence of an attachment point
        (normally a squiggle line perpendicular to a bond).
        The coordinates
        are in the molecule frame unless rawCoords is true,
        in which case the coordinates are in pixels.
        """

    def DrawWavyLine(self, cds1: rdkit.Geometry.rdGeometry.Point2D, cds2: rdkit.Geometry.rdGeometry.Point2D, color1: tuple[float, ...], color2: tuple[float, ...], nSegments: int = 16, vertOffset: float = 0.05, rawCoords: bool = False) -> None:
        """
        draw a line indicating the presence of an attachment point
        (normally a squiggle line perpendicular to a bond).
        The coordinates
        are in the molecule frame unless rawCoords is true,
        in which case the coordinates are in pixels.
        """

    @overload
    def DrawString(self, string: str, pos: rdkit.Geometry.rdGeometry.Point2D, rawCoords: bool = False) -> None:
        """
        add text to the canvas. The coordinates are in the molecule frame unless rawCoords is true, in which case the coordinates are in pixels.
        """

    @overload
    def DrawString(self, string: str, pos: rdkit.Geometry.rdGeometry.Point2D, align: int, rawCoords: bool = False) -> None:
        """
        add aligned text to the canvas. The align argument can be 0
        (=MIDDLE), 1 (=START), or 2 (=END).
        The coordinates
        are in the molecule frame unless rawCoords is true,
        in which case the coordinates are in pixels.
        """

    @overload
    def GetDrawCoords(self, point: rdkit.Geometry.rdGeometry.Point2D) -> rdkit.Geometry.rdGeometry.Point2D:
        """
        get the coordinates in drawing space for a particular point in molecule space
        """

    @overload
    def GetDrawCoords(self, atomIndex: int) -> rdkit.Geometry.rdGeometry.Point2D:
        """get the coordinates in drawing space for a particular atom"""

    def ClearDrawing(self) -> None:
        """clears the drawing by filling it with the background color"""

    def drawOptions(self) -> MolDrawOptions:
        """Returns a modifiable version of the current drawing options"""

    def SetDrawOptions(self, opts: MolDrawOptions) -> None:
        """Copies the drawing options passed in over our drawing options"""

class MolDraw2DSVG(MolDraw2D):
    """SVG molecule drawer"""

    def __init__(self, width: int, height: int, panelWidth: int = -1, panelHeight: int = -1, noFreetype: bool = False) -> None: ...

    def FinishDrawing(self) -> None:
        """add the last bits of SVG to finish the drawing"""

    def AddMoleculeMetadata(self, mol: rdkit.Chem.rdchem.Mol, confId: int = -1) -> None:
        """add RDKit-specific information to the bottom of the drawing"""

    def TagAtoms(self, mol: rdkit.Chem.rdchem.Mol, radius: float = 0.2, events: dict[str, str] | None = None) -> None:
        """allow atom selection in the SVG"""

    def GetDrawingText(self) -> str:
        """return the SVG"""

class MolDraw2DCairo(MolDraw2D):
    """Cairo molecule drawer"""

    def __init__(self, width: int, height: int, panelWidth: int = -1, panelHeight: int = -1, noFreetype: bool = False) -> None: ...

    def FinishDrawing(self) -> None:
        """add the last bits to finish the drawing"""

    def GetDrawingText(self) -> bytes:
        """return the PNG data as a string"""

    def WriteDrawingText(self, fName: str) -> None:
        """write the PNG data to the named file"""

def PrepareMolForDrawing(mol: rdkit.Chem.rdchem.Mol | None, kekulize: bool = True, addChiralHs: bool = True, wedgeBonds: bool = True, forceCoords: bool = False, wavyBonds: bool = False) -> rdkit.Chem.rdchem.Mol:
    """
    Does some cleanup operations on the molecule to prepare it to draw nicely.
    The operations include: kekulization, addition of chiral Hs (so that we can draw
    wedges to them), wedging of bonds at chiral centers, and generation of a 2D
    conformation if the molecule does not already have a conformation

    Returns a modified copy of the molecule.
    """

def PrepareAndDrawMolecule(drawer: MolDraw2D, mol: rdkit.Chem.rdchem.Mol, legend: str = '', highlightAtoms: Iterable[int] | None = None, highlightBonds: Iterable[int] | None = None, highlightAtomColors: dict[int, tuple[float, ...]] | None = None, highlightBondColors: dict[int, tuple[float, ...]] | None = None, highlightAtomRadii: dict[int, float] | None = None, confId: int = -1, kekulize: bool = True) -> None:
    """Preps a molecule for drawing and actually draws it"""

def DrawMoleculeACS1996(drawer: MolDraw2D, mol: rdkit.Chem.rdchem.Mol, legend: str = '', highlightAtoms: Iterable[int] | None = None, highlightBonds: Iterable[int] | None = None, highlightAtomColors: dict[int, tuple[float, ...]] | None = None, highlightBondColors: dict[int, tuple[float, ...]] | None = None, highlightAtomRadii: dict[int, float] | None = None, confId: int = -1) -> None:
    """Draws molecule in ACS 1996 mode."""

class ContourParams:
    """Parameters for drawing contours"""

    def __init__(self) -> None: ...

    @property
    def setScale(self) -> bool:
        """
        set the scale of the drawing object (useful if you draw the grid/contours first)
        """

    @setScale.setter
    def setScale(self, arg: bool, /) -> None: ...

    @property
    def dashNegative(self) -> bool:
        """use a dashed line for negative contours"""

    @dashNegative.setter
    def dashNegative(self, arg: bool, /) -> None: ...

    @property
    def fillGrid(self) -> bool:
        """colors the grid in addition to drawing contours"""

    @fillGrid.setter
    def fillGrid(self, arg: bool, /) -> None: ...

    @property
    def gridResolution(self) -> float:
        """set the resolution of the grid"""

    @gridResolution.setter
    def gridResolution(self, arg: float, /) -> None: ...

    @property
    def contourWidth(self) -> float:
        """line width of the contours"""

    @contourWidth.setter
    def contourWidth(self, arg: float, /) -> None: ...

    @property
    def extraGridPadding(self) -> float:
        """extra space (in molecule coords) around the grid"""

    @extraGridPadding.setter
    def extraGridPadding(self, arg: float, /) -> None: ...

    @property
    def drawAsLines(self) -> bool:
        """draw the contours as continuous lines isntead of line segments"""

    @drawAsLines.setter
    def drawAsLines(self, arg: bool, /) -> None: ...

    @property
    def coordScaleForQuantization(self) -> float:
        """
        scaling factor used to convert coordinates to ints when forming the continuous lines
        """

    @coordScaleForQuantization.setter
    def coordScaleForQuantization(self, arg: float, /) -> None: ...

    @property
    def isovalScaleForQuantization(self) -> float:
        """
        scaling factor used to convert isovalues to ints when forming the continuous lines
        """

    @isovalScaleForQuantization.setter
    def isovalScaleForQuantization(self, arg: float, /) -> None: ...

    @property
    def useFillThreshold(self) -> bool:
        """use a magnitude threshold to determine if a grid point is filled"""

    @useFillThreshold.setter
    def useFillThreshold(self, arg: bool, /) -> None: ...

    @property
    def fillThreshold(self) -> float:
        """magnitude threshold to determine if a grid point is filled"""

    @fillThreshold.setter
    def fillThreshold(self, arg: float, /) -> None: ...

    @property
    def fillThresholdIsFraction(self) -> bool:
        """if true, fillThreshold is a fraction of the range of the data"""

    @fillThresholdIsFraction.setter
    def fillThresholdIsFraction(self, arg: bool, /) -> None: ...

    @property
    def colourMap(self) -> tuple[tuple[float, float, float, float], ...]:
        """the color map to use when filling the grid"""

    @colourMap.setter
    def colourMap(self, arg: Sequence[tuple[float, ...]], /) -> None: ...

    @property
    def contourColour(self) -> tuple[float, float, float, float]:
        """the color to use for drawing the contours"""

    @contourColour.setter
    def contourColour(self, arg: tuple[float, ...], /) -> None: ...

    def setContourColour(self, colour: tuple[float, ...]) -> None: ...

    def setColourMap(self, colours: Sequence[tuple[float, ...]]) -> None: ...

    def __setattr__(self, name: str, value: object | None) -> None: ...

def ContourAndDrawGaussians(drawer: MolDraw2D, locs: Iterable[rdkit.Geometry.rdGeometry.Point2D], heights: Iterable[float], widths: Iterable[float], nContours: int = 10, levels: Iterable[float] | None = None, params: ContourParams = ..., mol: rdkit.Chem.rdchem.Mol | None = None) -> None:
    """
    Generates and draws contours for a set of gaussians

    - drawer: the MolDraw2D object to use
    - locs: locations of the gaussians
    - heights: the heights (or weights) of the gaussians
    - widths: the standard deviations of the gaussians
    - nContours: the number of contours to draw
    - levels: the contours to use
    - ps: additional parameters controlling the contouring.
    - mol: molecule used to help set scale.

    The values are calculated on a grid with spacing params.gridResolution.
    If params.setScale  is set, the grid size will be calculated based on the
    locations of the gaussians and params.extraGridPadding. Otherwise the current
    size of the viewport will be used.

    If the levels argument is empty, the contour levels will be determined
    automatically from the max and min values on the grid and levels will
    be updated to include the contour levels.

    If params.fillGrid is set, the data on the grid will also be drawn using
    the color scheme in params.colourMap

    If mol is not 0, uses the molecule to help set the scale, assuming that
    it will be drawn over the plot, so needs to fit on it.
    """

def ContourAndDrawGrid(drawer: MolDraw2D, data: Annotated[NDArray[numpy.float64], dict(shape=(None, None), order='C')], xcoords: Iterable[float], ycoords: Iterable[float], nContours: int = 10, levels: Iterable[float] | None = None, params: ContourParams = ..., mol: rdkit.Chem.rdchem.Mol | None = None) -> None:
    """
    Generates and draws contours for data on a grid

    - drawer: the MolDraw2D object to use
    - data: numpy array with the data to be contoured
    - xcoords: the x coordinates of the grid
    - ycoords: the y coordinates of the grid
    - nContours: the number of contours to draw
    - levels: the contours to use
    - ps: additional parameters controlling the contouring
    - mol: molecule used to help set scale.

    The values are calculated on a grid with spacing params.gridResolution.
    If params.setScale  is set, the grid size will be calculated based on the
    locations of the gaussians and params.extraGridPadding. Otherwise the current
    size of the viewport will be used.

    If the levels argument is empty, the contour levels will be determined
    automatically from the max and min values on the grid and levels will
    be updated to include the contour levels.

    If params.fillGrid is set, the data on the grid will also be drawn using
    the color scheme in params.colourMap

    If mol is not 0, uses the molecule to help set the scale, assuming that
    it will be drawn over the plot, so needs to fit on it.
    """

def UpdateMolDrawOptionsFromJSON(opts: MolDrawOptions, json: str | bytes) -> None: ...

def UpdateDrawerParamsFromJSON(drawer: MolDraw2D, json: str | bytes) -> None: ...

def MolToSVG(mol: rdkit.Chem.rdchem.Mol, width: int = 300, height: int = 300, highlightAtoms: Iterable[int] | None = None, kekulize: bool = True, lineWidthMult: int = 1, includeAtomCircles: bool = True, confId: int = -1) -> str:
    """Returns svg for a molecule"""

def MolToACS1996SVG(mol: rdkit.Chem.rdchem.Mol, legend: str = '', highlightAtoms: Iterable[int] | None = None, highlightBonds: Iterable[int] | None = None, highlightAtomColors: dict[int, tuple[float, ...]] | None = None, highlightBondColors: dict[int, tuple[float, ...]] | None = None, highlightAtomRadii: dict[int, float] | None = None, confId: int = -1) -> str:
    """Returns ACS 1996 mode svg for a molecule"""

def SetACS1996Mode(drawOptions: MolDrawOptions, meanBondLength: float) -> None:
    """
    Set the draw options to produce something as close as possible to
    the ACS 1996 guidelines as described at
    https://en.wikipedia.org/wiki/Wikipedia:Manual_of_Style/Chemistry/Structure_drawing

    - MolDrawOptions opt - the options what will be changed
    - float meanBondLength - mean bond length of the molecule

    Works best if the MolDraw2D object is created with width and height -1 (a
    flexiCanvas).
    The mean bond length may be calculated with MeanBondLength.
    It is used to calculate the offset for the lines in multiple bonds.

    Options changed are:
      bondLineWidth = 0.6
      scaleBondWidth = false
      scalingFactor = 14.4 / meanBondLen
      multipleBondOffset = 0.18
      highlightBondWidthMultiplier = 32
      setMonochromeMode - black and white
      fixedFontSize = 10
      additionalAtomLabelPadding = 0.066
      fontFile - if it isn't set already, then if RDBASE is set and the file
                 exists, uses $RDBASE/Data/Fonts/FreeSans.ttf.  Otherwise uses
                 BuiltinRobotoRegular.
    """

def MeanBondLength(mol: rdkit.Chem.rdchem.Mol, confId: int = -1) -> float:
    """Calculate the mean bond length for the molecule."""

@overload
def SetDarkMode(d2d: MolDrawOptions) -> None:
    """set dark mode for a MolDrawOptions object"""

@overload
def SetDarkMode(d2d: MolDraw2D) -> None:
    """set dark mode for a MolDraw2D object"""

@overload
def SetMonochromeMode(options: MolDrawOptions, fgColour: tuple[float, ...], bgColour: tuple[float, ...]) -> None:
    """set monochrome mode for a MolDrawOptions object"""

@overload
def SetMonochromeMode(drawer: MolDraw2D, fgColour: tuple[float, ...], bgColour: tuple[float, ...]) -> None:
    """set monochrome mode for a MolDraw2D object"""
