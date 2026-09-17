"""
Module containing functions to set partial charges - currently Gasteiger Charges
"""

import rdkit.Chem.rdchem


def ComputeGasteigerCharges(mol: rdkit.Chem.rdchem.Mol, nIter: int = 12, throwOnParamFailure: bool = False) -> None:
    """
    Compute Gasteiger partial charges for molecule

    The charges are computed using an iterative procedure presented in

    Ref : J.Gasteiger, M. Marseli, Iterative Equalization of Oribital Electronegatiity
    A Rapid Access to Atomic Charges, Tetrahedron Vol 36 p3219 1980

    The computed charges are stored on each atom are stored a computed property (under the name
    _GasteigerCharge). In addition, each atom also stored the total charge for the implicit hydrogens
    on the atom (under the property name _GasteigerHCharge)

    ARGUMENTS:

       - mol : the molecule of interrest
       - nIter : number of iteration (defaults to 12)
       - throwOnParamFailure : toggles whether or not an exception should be raised if parameters
         for an atom cannot be found.  If this is false (the default), all parameters for unknown
         atoms will be set to zero.  This has the effect of removing that atom from the iteration.
    """
