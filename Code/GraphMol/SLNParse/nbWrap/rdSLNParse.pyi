"""
Module containing classes and functions for working with Sybyl line notation (SLN).
"""

import rdkit.Chem.rdchem


class SLNParseException(ValueError):
    pass

def MolFromSLN(SLN: str, sanitize: bool = True, debugParser: bool = False) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a molecule from an SLN string.

        ARGUMENTS:

        - SLN: the SLN string

        - sanitize: (optional) toggles sanitization of the molecule.
          Defaults to True.

      RETURNS:

        a Mol object, None on failure.

      NOTE: the SLN should not contain query information or properties. To build a
        query from SLN, use MolFromQuerySLN.
    """

def MolFromQuerySLN(SLN: str, mergeHs: bool = True, debugParser: bool = False) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a query molecule from an SLN string.

      ARGUMENTS:

        - SLN: the SLN string

        - mergeHs: (optional) toggles the merging of explicit Hs in the query into the attached
          heavy atoms. Defaults to False.

      RETURNS:

        a Mol object suitable for using in substructure queries, None on failure.
    """
