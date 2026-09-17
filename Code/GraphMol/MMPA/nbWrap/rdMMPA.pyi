"""Module containing a C++ implementation of code for doing MMPA"""

from collections.abc import Sequence
from typing import overload

import rdkit.Chem.rdchem


@overload
def FragmentMol(mol: rdkit.Chem.rdchem.Mol, maxCuts: int = 3, maxCutBonds: int = 20, pattern: str = '[#6+0;!$(*=,#[!#6])]!@!=!#[*]', resultsAsMols: bool = True) -> tuple[tuple, ...]: ...

@overload
def FragmentMol(mol: rdkit.Chem.rdchem.Mol, minCuts: int, maxCuts: int, maxCutBonds: int, pattern: str = '[#6+0;!$(*=,#[!#6])]!@!=!#[*]', resultsAsMols: bool = True) -> tuple[tuple, ...]: ...

@overload
def FragmentMol(mol: rdkit.Chem.rdchem.Mol, bondsToCut: Sequence[int], minCuts: int = 1, maxCuts: int = 3, resultsAsMols: bool = True) -> tuple[tuple, ...]:
    """Does the fragmentation necessary for an MMPA analysis"""
