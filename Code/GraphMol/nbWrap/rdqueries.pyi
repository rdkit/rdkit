"""Module containing RDKit functionality for querying molecules."""

import rdkit.Chem.rdchem
import rdkit.DataStructs.cDataStructs


def AtomNumEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where AtomNum is equal to the target value.
    """

def AtomNumLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where AtomNum is less than the target value.                      \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def AtomNumGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where AtomNum is equal to the target value.                       \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def ExplicitValenceEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where ExplicitValence is equal to the target value.
    """

def ExplicitValenceLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where ExplicitValence is less than the target value.                      \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def ExplicitValenceGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where ExplicitValence is equal to the target value.                       \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def TotalValenceEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where TotalValence is equal to the target value.
    """

def TotalValenceLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where TotalValence is less than the target value.                      \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def TotalValenceGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where TotalValence is equal to the target value.                       \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def ExplicitDegreeEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where ExplicitDegree is equal to the target value.
    """

def ExplicitDegreeLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where ExplicitDegree is less than the target value.                      \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def ExplicitDegreeGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where ExplicitDegree is equal to the target value.                       \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def TotalDegreeEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where TotalDegree is equal to the target value.
    """

def TotalDegreeLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where TotalDegree is less than the target value.                      \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def TotalDegreeGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where TotalDegree is equal to the target value.                       \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def HCountEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where HCount is equal to the target value.
    """

def HCountLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where HCount is less than the target value.                      \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def HCountGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where HCount is equal to the target value.                       \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def MassEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where Mass is equal to the target value.
    """

def MassLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where Mass is less than the target value.                      \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def MassGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where Mass is equal to the target value.                       \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def IsotopeEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where Isotope is equal to the target value.
    """

def IsotopeLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where Isotope is less than the target value.                      \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def IsotopeGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where Isotope is equal to the target value.                       \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def FormalChargeEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where FormalCharge is equal to the target value.
    """

def FormalChargeLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where FormalCharge is less than the target value.                      \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def FormalChargeGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where FormalCharge is equal to the target value.                       \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def HybridizationEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where Hybridization is equal to the target value.
    """

def HybridizationLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where Hybridization is less than the target value.                      \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def HybridizationGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where Hybridization is equal to the target value.                       \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def InNRingsEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where InNRings is equal to the target value.
    """

def InNRingsLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where InNRings is less than the target value.                      \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def InNRingsGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where InNRings is equal to the target value.                       \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def MinRingSizeEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where MinRingSize is equal to the target value.
    """

def MinRingSizeLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where MinRingSize is less than the target value.                      \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def MinRingSizeGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where MinRingSize is equal to the target value.                       \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def RingBondCountEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where RingBondCount is equal to the target value.
    """

def RingBondCountLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where RingBondCount is less than the target value.                      \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def RingBondCountGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where RingBondCount is equal to the target value.                       \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def NumRadicalElectronsEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where NumRadicalElectrons is equal to the target value.
    """

def NumRadicalElectronsLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where NumRadicalElectrons is less than the target value.                      \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def NumRadicalElectronsGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where NumRadicalElectrons is equal to the target value.                       \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def NumHeteroatomNeighborsEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where NumHeteroatomNeighbors is equal to the target value.
    """

def NumHeteroatomNeighborsLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where NumHeteroatomNeighbors is less than the target value.                      \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def NumHeteroatomNeighborsGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where NumHeteroatomNeighbors is equal to the target value.                       \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def NumAliphaticHeteroatomNeighborsEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where NumAliphaticHeteroatomNeighbors is equal to the target value.
    """

def NumAliphaticHeteroatomNeighborsLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where NumAliphaticHeteroatomNeighbors is less than the target value.                      \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def NumAliphaticHeteroatomNeighborsGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where NumAliphaticHeteroatomNeighbors is equal to the target value.                       \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def NonHydrogenDegreeEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where NonHydrogenDegree is equal to the target value.
    """

def NonHydrogenDegreeLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where NonHydrogenDegree is less than the target value.                      \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def NonHydrogenDegreeGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches atoms where NonHydrogenDegree is equal to the target value.                       \\
    NOTE: the direction of comparison is reversed relative to the C++ API
    """

def IsUnsaturatedQueryAtom(negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """Returns a QueryAtom that matches atoms when IsUnsaturated is True."""

def IsAromaticQueryAtom(negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """Returns a QueryAtom that matches atoms when IsAromatic is True."""

def IsAliphaticQueryAtom(negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """Returns a QueryAtom that matches atoms when IsAliphatic is True."""

def IsInRingQueryAtom(negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """Returns a QueryAtom that matches atoms when IsInRing is True."""

def HasChiralTagQueryAtom(negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """Returns a QueryAtom that matches atoms when HasChiralTag is True."""

def MissingChiralTagQueryAtom(negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """Returns a QueryAtom that matches atoms when MissingChiralTag is True."""

def IsBridgeheadQueryAtom(negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """Returns a QueryAtom that matches atoms when IsBridgehead is True."""

def AAtomQueryAtom(negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """Returns a QueryAtom that matches atoms when AAtom is True."""

def AHAtomQueryAtom(negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """Returns a QueryAtom that matches atoms when AHAtom is True."""

def XAtomQueryAtom(negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """Returns a QueryAtom that matches atoms when XAtom is True."""

def XHAtomQueryAtom(negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """Returns a QueryAtom that matches atoms when XHAtom is True."""

def QAtomQueryAtom(negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """Returns a QueryAtom that matches atoms when QAtom is True."""

def QHAtomQueryAtom(negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """Returns a QueryAtom that matches atoms when QHAtom is True."""

def MAtomQueryAtom(negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """Returns a QueryAtom that matches atoms when MAtom is True."""

def MHAtomQueryAtom(negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """Returns a QueryAtom that matches atoms when MHAtom is True."""

def HasPropQueryAtom(propname: str, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches when the property 'propname'
    exists in the atom.
    """

def HasIntPropWithValueQueryAtom(propname: str, val: int, negate: bool = False, tolerance: int = 0) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches when the property 'propname'
    has the specified int value.
    """

def HasBoolPropWithValueQueryAtom(propname: str, val: bool, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches when the property 'propname'
    has the specified boolean value.
    """

def HasStringPropWithValueQueryAtom(propname: str, val: str, negate: bool = False) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches when the property 'propname'
    has the specified string value.
    """

def HasDoublePropWithValueQueryAtom(propname: str, val: float, negate: bool = False, tolerance: float = 0.0) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches when the property 'propname'
    has the specified value +- tolerance
    """

def HasBitVectPropWithValueQueryAtom(propname: str, val: rdkit.DataStructs.cDataStructs.ExplicitBitVect, negate: bool = False, tolerance: float = 0) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Returns a QueryAtom that matches when the property 'propname'
    has the specified explicit bit vector value.  The Tolerance is the allowed
    Tanimoto difference
    """

def HasPropQueryBond(propname: str, negate: bool = False) -> rdkit.Chem.rdchem.QueryBond:
    """
    Returns a QueryBond that matches when the property 'propname'
    exists in the bond.
    """

def HasIntPropWithValueQueryBond(propname: str, val: int, negate: bool = False, tolerance: int = 0) -> rdkit.Chem.rdchem.QueryBond:
    """
    Returns a QueryBond that matches when the property 'propname'
    has the specified int value.
    """

def HasBoolPropWithValueQueryBond(propname: str, val: bool, negate: bool = False) -> rdkit.Chem.rdchem.QueryBond:
    """
    Returns a QueryBond that matches when the property 'propname'
    has the specified boolean value.
    """

def HasStringPropWithValueQueryBond(propname: str, val: str, negate: bool = False) -> rdkit.Chem.rdchem.QueryBond:
    """
    Returns a QueryBond that matches when the property 'propname'
    has the specified string value.
    """

def HasDoublePropWithValueQueryBond(propname: str, val: float, negate: bool = False, tolerance: float = 0.0) -> rdkit.Chem.rdchem.QueryBond:
    """
    Returns a QueryBond that matches when the property 'propname'
    has the specified value +- tolerance
    """

def ReplaceAtomWithQueryAtom(mol: rdkit.Chem.rdchem.Mol, atom: rdkit.Chem.rdchem.Atom) -> rdkit.Chem.rdchem.QueryAtom:
    """
    Changes the given atom in the molecule to a query atom and returns
    the atom which can then be modified, for example with additional query
    constraints added.  The new atom is otherwise a copy of the old.
    If the atom already has a query, nothing will be changed.
    """
