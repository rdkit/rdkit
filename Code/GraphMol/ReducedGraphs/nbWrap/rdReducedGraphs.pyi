"""Module containing functions to generate and work with reduced graphs"""

from typing import Annotated

import numpy
from numpy.typing import NDArray

import rdkit.Chem.rdchem


def GenerateMolExtendedReducedGraph(mol: rdkit.Chem.rdchem.Mol, atomTypes: object | None = None) -> rdkit.Chem.rdchem.Mol:
    """Returns the reduced graph for a molecule"""

def GenerateErGFingerprintForReducedGraph(mol: rdkit.Chem.rdchem.Mol, atomTypes: object | None = None, fuzzIncrement: float = 0.3, minPath: int = 1, maxPath: int = 15) -> Annotated[NDArray[numpy.float64], dict(shape=(None,))]:
    """Returns the ErG fingerprint vector for a reduced graph"""

def GetErGFingerprint(mol: rdkit.Chem.rdchem.Mol, atomTypes: object | None = None, fuzzIncrement: float = 0.3, minPath: int = 1, maxPath: int = 15) -> Annotated[NDArray[numpy.float64], dict(shape=(None,))]:
    """Returns the ErG fingerprint vector for a molecule"""
