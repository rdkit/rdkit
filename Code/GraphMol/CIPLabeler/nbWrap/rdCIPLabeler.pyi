"""
Module containing a function to assign stereochemical labels"
based on an accurate CIP rules implementation. This algoritm is a port of
https://github.com/SiMolecule/centres, which was originally written by John
Mayfield. The original algorithm is described in:

Hanson, R. M., Musacchio, S., Mayfield, J. W., Vainio, M. J., Yerin, A., Redkin, D.
Algorithmic Analysis of Cahn--Ingold--Prelog Rules of Stereochemistry:
Proposals for Revised Rules and a Guide for Machine Implementation.
J. Chem. Inf. Model. 2018, 58, 1755-1765.
"""

from collections.abc import Iterable

import rdkit.Chem.rdchem


class MaxIterationsExceeded(RuntimeError):
    pass

def AssignCIPLabels(mol: rdkit.Chem.rdchem.Mol, atomsToLabel: Iterable[int] | None = None, bondsToLabel: Iterable[int] | None = None, maxRecursiveIterations: int = 0) -> None:
    """
    New implementation of Stereo assignment using a true CIP ranking.
    On return:  The molecule to contains CIP flags
    Errors:  when maxRecursiveIterations is exceeded, throws a MaxIterationsExceeded error
    ARGUMENTS:

     - mol: the molecule
     - atomsToLabel: (optional) list of atoms to label
     - bondsToLabel: (optional) list of bonds to label
     - maxRecursiveIterations: (optional) protects against pseudo-infinite
    recursion for highly symmetrical structures.
     A value of 1,250,000 take about 1 second.  Most structures requires less than 10,000iterations.
     A peptide with MW~3000 took about 100 iterations, and a 20,000 mw protein took about 600 iterations
    (0 = default - no limit)
    """
