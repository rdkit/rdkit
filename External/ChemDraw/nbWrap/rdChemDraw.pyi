"""
Module containing classes and functions for working with ChemDraw files.
"""

import enum
import os

import rdkit.Chem.rdchem


class CDXFormat(enum.Enum):
    CDX = 1

    CDXML = 2

class NeedsCleanPolicy(enum.Enum):
    TrustSource = 0

    TrustExplicitHydrogens = 1

def MolsFromChemDrawFile(filename: str | os.PathLike, sanitize: bool = True, removeHs: bool = True, needsCleanPolicy: NeedsCleanPolicy = NeedsCleanPolicy.TrustSource) -> tuple:
    """
    Extract all molecules from a ChemDraw file.

    Note that the ChemDraw format is large and complex, the RDKit doesn't support
    full functionality, just the base ones required for molecule and
    reaction parsing.

    ARGUMENTS:

      - filename: the chemdraw filename (.cdx/.cdxml)

      - sanitize: if True, sanitize the molecules [default True]

      - removeHs: if True, convert explicit Hs into implicit Hs. [default True]

      - needsCleanPolicy: how to handle `NeedsClean` hydrogen metadata.
      `TrustSource` honors `NeedsClean` by allowing sanitization to
      recompute hydrogens. `TrustExplicitHydrogens` preserves the literal
      source metadata when sanitize is True. [default TrustSource]

    RETURNS:
      a tuple of parsed Mol objects.
    """

def MolsFromChemDrawBlock(block: str | bytes, sanitize: bool = True, removeHs: bool = True, needsCleanPolicy: NeedsCleanPolicy = NeedsCleanPolicy.TrustSource) -> tuple:
    """
    Extract all molecules from a ChemDraw block.

    Note that the ChemDraw format is large and complex, the RDKit doesn't support
    full functionality, just the base ones required for molecule and
    reaction parsing.

    ARGUMENTS:

      - block: the CDX/CDXML block

      - sanitize: if True, sanitize the molecules [default True]

      - removeHs: if True, convert explicit Hs into implicit Hs. [default True]

      - needsCleanPolicy: how to handle `NeedsClean` hydrogen metadata.
      `TrustSource` honors `NeedsClean` by allowing sanitization to
      recompute hydrogens. `TrustExplicitHydrogens` preserves the literal
      source metadata when sanitize is True. [default TrustSource]

    RETURNS:
      a tuple of parsed Mol objects.
    """

def ReactionsFromChemDrawFile(filename: str | os.PathLike, sanitize: bool = False, removeHs: bool = False) -> tuple:
    """
    Extract all reactions from a ChemDraw file.

    Note that the ChemDraw format is large and complex, the RDKit doesn't support
    full functionality, just the base ones required for molecule and
    reaction parsing.

    ARGUMENTS:

      - filename: the chemdraw filename (.cdx/.cdxml)

      - sanitize: if True, sanitize the molecules [default True]

      - removeHs: if True, convert explicit Hs into implicit Hs. [default True]

    RETURNS:
      a tuple of parsed ChemicalReaction objects.
    """

def ReactionsFromChemDrawBlock(rxnblock: str | bytes, sanitize: bool = False, removeHs: bool = False) -> tuple:
    """
    Extract all reactions from a ChemDraw text block.

    Note that the ChemDraw format is large and complex, the RDKit doesn't support
    full functionality, just the base ones required for molecule and
    reaction parsing.

    ARGUMENTS:

      - rxnblock: the ChemDraw text block

      - sanitize: if True, sanitize the molecules [default True]

      - removeHs: if True, convert explicit Hs into implicit Hs. [default True]

    RETURNS:
      a tuple of parsed ChemicalReaction objects.
    """

def MolToChemDrawBlock(mol: rdkit.Chem.rdchem.Mol, format: CDXFormat = CDXFormat.CDXML) -> str:
    """
    Convert a molecule into a chemdraw string using the specified format

    ARGUMENTS:

      - mol: the molecule to convert

      - format: The ChemDraw format to use, CDXML/CDX [default CDXML]

    RETURNS:
      the ChemDraw string.
    """
