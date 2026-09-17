"""
Module containing classes and functions for working with chemical reactions.
"""

from collections.abc import Iterable, Mapping, Sequence
import enum
from typing import overload

import rdkit.Chem.rdchem
import rdkit.Chem.rdmolfiles
import rdkit.Chem.rdmolops
import rdkit.DataStructs.cDataStructs


class ChemicalReactionParserException(ValueError):
    pass

class ChemicalReactionException(ValueError):
    pass

class FingerprintType(enum.Enum):
    AtomPairFP = 1

    TopologicalTorsion = 2

    MorganFP = 3

    RDKitFP = 4

    PatternFP = 5

class ReactionFingerprintParams:
    """
    A class for storing parameters to manipulate the calculation of
    fingerprints of chemical reactions.
    """

    @overload
    def __init__(self) -> None:
        """Constructor, takes no arguments"""

    @overload
    def __init__(self, includeAgents: bool, bitRatioAgents: float, nonAgentWeight: int, agentWeight: int, fpSize: int, fpType: FingerprintType) -> None: ...

    @property
    def fpSize(self) -> int: ...

    @fpSize.setter
    def fpSize(self, arg: int, /) -> None: ...

    @property
    def fpType(self) -> FingerprintType: ...

    @fpType.setter
    def fpType(self, arg: FingerprintType, /) -> None: ...

    @property
    def bitRatioAgents(self) -> float: ...

    @bitRatioAgents.setter
    def bitRatioAgents(self, arg: float, /) -> None: ...

    @property
    def nonAgentWeight(self) -> int: ...

    @nonAgentWeight.setter
    def nonAgentWeight(self, arg: int, /) -> None: ...

    @property
    def agentWeight(self) -> int: ...

    @agentWeight.setter
    def agentWeight(self, arg: int, /) -> None: ...

    @property
    def includeAgents(self) -> bool: ...

    @includeAgents.setter
    def includeAgents(self, arg: bool, /) -> None: ...

    def __setattr__(self, name: str, value: object | None) -> None: ...

class ChemicalReaction:
    """
    A class for storing and applying chemical reactions.

    Sample Usage:
      >>> from rdkit import Chem
      >>> from rdkit.Chem import rdChemReactions
      >>> rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]')
      >>> reacts = (Chem.MolFromSmiles('C(=O)O'),Chem.MolFromSmiles('CNC'))
      >>> products = rxn.RunReactants(reacts)
      >>> len(products)
      1
      >>> len(products[0])
      1
      >>> Chem.MolToSmiles(products[0][0])
      'CN(C)C=O'
    """

    @overload
    def __init__(self) -> None:
        """Constructor, takes no arguments"""

    @overload
    def __init__(self, binStr: bytes) -> None: ...

    @overload
    def __init__(self, binStr: str) -> None: ...

    @overload
    def __init__(self, other: ChemicalReaction) -> None: ...

    def GetNumReactantTemplates(self) -> int:
        """returns the number of reactants this reaction expects"""

    def GetNumProductTemplates(self) -> int:
        """returns the number of products this reaction generates"""

    def GetNumAgentTemplates(self) -> int:
        """returns the number of agents this reaction expects"""

    def AddReactantTemplate(self, mol: rdkit.Chem.rdchem.Mol) -> int:
        """adds a reactant (a Molecule) to the reaction"""

    def AddProductTemplate(self, mol: rdkit.Chem.rdchem.Mol) -> int:
        """adds a product (a Molecule)"""

    def AddAgentTemplate(self, mol: rdkit.Chem.rdchem.Mol) -> int:
        """adds a agent (a Molecule)"""

    def RemoveUnmappedReactantTemplates(self, thresholdUnmappedAtoms: float = 0.2, moveToAgentTemplates: bool = True, targetList: list[rdkit.Chem.rdchem.Mol] | None = None) -> None:
        """
        Removes molecules with an atom mapping ratio below
        thresholdUnmappedAtoms from reactant templates to the agent
        templates or to a given targetList
        """

    def RemoveUnmappedProductTemplates(self, thresholdUnmappedAtoms: float = 0.2, moveToAgentTemplates: bool = True, targetList: list[rdkit.Chem.rdchem.Mol] | None = None) -> None:
        """
        Removes molecules with an atom mapping ratio below
        thresholdUnmappedAtoms from product templates to the agent
        templates or to a given targetList
        """

    def RemoveAgentTemplates(self, targetList: list[rdkit.Chem.rdchem.Mol] | None = None) -> None:
        """
        Removes agents from reaction. If targetList is provide the agents
        will be transferred to that list.
        """

    def RunReactants(self, reactants: Sequence[rdkit.Chem.rdchem.Mol], maxProducts: int = 1000) -> tuple[tuple[rdkit.Chem.rdchem.Mol, ...], ...]:
        """
        apply the reaction to a sequence of reactant molecules and return
        the products as a tuple of tuples.  If maxProducts is not zero,
         stop the reaction when maxProducts have been generated [default=1000]
        """

    def RunReactant(self, reactant: rdkit.Chem.rdchem.Mol, reactionIdx: int) -> tuple[tuple[rdkit.Chem.rdchem.Mol, ...], ...]:
        """apply the reaction to a single reactant"""

    def RunReactantInPlace(self, reactant: rdkit.Chem.rdchem.Mol, removeUnmatchedAtoms: bool = True) -> bool:
        """
        apply the reaction to a single reactant in place. The reactant
        itself is modified. This can only be used for single reactant -
        single product reactions.
        """

    def Initialize(self, silent: bool = False) -> None:
        """initializes the reaction so that it can be used"""

    def IsInitialized(self) -> bool:
        """checks if the reaction is ready for use"""

    def Validate(self, silent: bool = False) -> tuple[int, int]:
        """
        checks the reaction for potential problems, returns (numWarnings,numErrors)
        """

    def GetProductTemplate(self, which: int) -> rdkit.Chem.rdchem.Mol:
        """returns one of our product templates"""

    def GetReactantTemplate(self, which: int) -> rdkit.Chem.rdchem.Mol:
        """returns one of our reactant templates"""

    def GetAgentTemplate(self, which: int) -> rdkit.Chem.rdchem.Mol:
        """returns one of our agent templates"""

    @overload
    def ToBinary(self) -> bytes: ...

    @overload
    def ToBinary(self, propertyFlags: object) -> bytes:
        """Returns a binary string representation of the reaction."""

    def IsMoleculeReactant(self, mol: rdkit.Chem.rdchem.Mol) -> bool:
        """
        returns whether or not the molecule has a substructure match to one of the reactants.
        """

    def IsMoleculeProduct(self, mol: rdkit.Chem.rdchem.Mol) -> bool:
        """
        returns whether or not the molecule has a substructure match to one of the products.
        """

    def IsMoleculeAgent(self, mol: rdkit.Chem.rdchem.Mol) -> bool:
        """
        returns whether or not the molecule has a substructure match to one of the agents.
        """

    def GetReactingAtoms(self, mappedAtomsOnly: bool = False) -> tuple[tuple[int, ...], ...]:
        """
        returns a sequence of sequences with the atoms that change in the reaction
        """

    def AddRecursiveQueriesToReaction(self, queries: dict[str, rdkit.Chem.rdchem.Mol] = {}, propName: str = 'molFileValue', getLabels: bool = False) -> object:
        """adds recursive queries and returns reactant labels"""

    def GetReactants(self) -> list[rdkit.Chem.rdchem.Mol]:
        """get the reactant templates"""

    def GetProducts(self) -> list[rdkit.Chem.rdchem.Mol]:
        """get the product templates"""

    def GetAgents(self) -> list[rdkit.Chem.rdchem.Mol]:
        """get the agent templates"""

    def GetSubstructParams(self) -> rdkit.Chem.rdchem.SubstructMatchParameters:
        """get the parameter object controlling the substructure matching"""

    def SetProp(self, key: str, val: str, computed: bool = False) -> None:
        """
        Sets a molecular property

          ARGUMENTS:
            - key: the name of the property to be set (a string).
            - value: the property value (a string).
            - computed: (optional) marks the property as being computed.
                        Defaults to False.
        """

    def SetDoubleProp(self, key: str, val: float, computed: bool = False) -> None:
        """
        Sets a double valued molecular property

          ARGUMENTS:
            - key: the name of the property to be set (a string).
            - value: the property value as a double.
            - computed: (optional) marks the property as being computed.
                        Defaults to 0.
        """

    def SetIntProp(self, key: str, val: int, computed: bool = False) -> None:
        """
        Sets an integer valued molecular property

          ARGUMENTS:
            - key: the name of the property to be set (an unsigned number).
            - value: the property value as an integer.
            - computed: (optional) marks the property as being computed.
                        Defaults to False.
        """

    def SetUnsignedProp(self, key: str, val: int, computed: bool = False) -> None:
        """
        Sets an unsigned integer valued molecular property

          ARGUMENTS:
            - key: the name of the property to be set (a string).
            - value: the property value as an unsigned integer.
            - computed: (optional) marks the property as being computed.
                        Defaults to False.
        """

    def SetBoolProp(self, key: str, val: bool, computed: bool = False) -> None:
        """
        Sets a boolean valued molecular property

          ARGUMENTS:
            - key: the name of the property to be set (a string).
            - value: the property value as a bool.
            - computed: (optional) marks the property as being computed.
                        Defaults to False.
        """

    def HasProp(self, key: str) -> int:
        """
        Queries a molecule to see if a particular property has been assigned.

          ARGUMENTS:
            - key: the name of the property to check for (a string).
        """

    def GetProp(self, key: str) -> object:
        """
        Returns the value of the property.

          ARGUMENTS:
            - key: the name of the property to return (a string).

          RETURNS: a string

          NOTE:
            - If the property has not been set, a KeyError exception will be raised.
        """

    def GetDoubleProp(self, key: str) -> object:
        """
        Returns the double value of the property if possible.

          ARGUMENTS:
            - key: the name of the property to return (a string).

          RETURNS: a double

          NOTE:
            - If the property has not been set, a KeyError exception will be raised.
        """

    def GetIntProp(self, key: str) -> object:
        """
        Returns the integer value of the property if possible.

          ARGUMENTS:
            - key: the name of the property to return (a string).

          RETURNS: an integer

          NOTE:
            - If the property has not been set, a KeyError exception will be raised.
        """

    def GetUnsignedProp(self, key: str) -> object:
        """
        Returns the unsigned int value of the property if possible.

          ARGUMENTS:
            - key: the name of the property to return (a string).

          RETURNS: an unsigned integer

          NOTE:
            - If the property has not been set, a KeyError exception will be raised.
        """

    def GetBoolProp(self, key: str) -> object:
        """
        Returns the Bool value of the property if possible.

          ARGUMENTS:
            - key: the name of the property to return (a string).

          RETURNS: a bool

          NOTE:
            - If the property has not been set, a KeyError exception will be raised.
        """

    def ClearProp(self, key: str) -> None:
        """
        Removes a property from the reaction.

          ARGUMENTS:
            - key: the name of the property to clear (a string).
        """

    def ClearComputedProps(self) -> None:
        """Removes all computed properties from the reaction."""

    def GetPropNames(self, includePrivate: bool = False, includeComputed: bool = False) -> list[str]:
        """
        Returns a tuple with all property names for this reaction.

          ARGUMENTS:
            - includePrivate: (optional) toggles inclusion of private properties in the result set.
                              Defaults to 0.
            - includeComputed: (optional) toggles inclusion of computed properties in the result set.
                              Defaults to 0.

          RETURNS: a tuple of strings
        """

    def GetPropsAsDict(self, includePrivate: bool = False, includeComputed: bool = False, autoConvertStrings: bool = True) -> dict:
        """
        Returns a dictionary populated with the reaction's properties.
         n.b. Some properties are not able to be converted to python types.

          ARGUMENTS:
            - includePrivate: (optional) toggles inclusion of private properties in the result set.
                              Defaults to False.
            - includeComputed: (optional) toggles inclusion of computed properties in the result set.
                              Defaults to False.

          RETURNS: a dictionary
        """

    def __setstate__(self, arg: object, /) -> None: ...

    def __getstate__(self) -> tuple[bytes, dict]: ...

def ReactionFromSmarts(SMARTS: str, replacements: Mapping[str, str] = {}, useSmiles: bool = False) -> ChemicalReaction:
    """
    construct a ChemicalReaction from a reaction SMARTS string.
    see the documentation for rdkit.Chem.MolFromSmiles for an explanation
    of the replacements argument.
    """

@overload
def ReactionToSmarts(reaction: ChemicalReaction) -> str: ...

@overload
def ReactionToSmarts(reaction: ChemicalReaction, params: rdkit.Chem.rdmolfiles.SmilesWriteParams) -> str:
    """construct a reaction SMARTS string for a ChemicalReaction"""

def ReactionFromSmiles(SMILES: str, replacements: Mapping[str, str] = {}) -> ChemicalReaction:
    """
    construct a ChemicalReaction from a reaction SMILES string.
    see the documentation for rdkit.Chem.MolFromSmiles for an explanation
    of the replacements argument.
    """

@overload
def ReactionToSmiles(reaction: ChemicalReaction, canonical: bool = True) -> str: ...

@overload
def ReactionToSmiles(reaction: ChemicalReaction, params: rdkit.Chem.rdmolfiles.SmilesWriteParams) -> str:
    """construct a reaction SMILES string for a ChemicalReaction"""

@overload
def ReactionToCXSmarts(reaction: ChemicalReaction) -> str:
    """construct a reaction SMARTS string for a ChemicalReaction"""

@overload
def ReactionToCXSmarts(reaction: ChemicalReaction, params: rdkit.Chem.rdmolfiles.SmilesWriteParams, flags: int = CXSmilesFields.CX_ALL) -> str:
    """construct a reaction CXSMARTS string for a ChemicalReaction"""

@overload
def ReactionToCXSmiles(reaction: ChemicalReaction, canonical: bool = True) -> str:
    """construct a reaction SMILES string for a ChemicalReaction"""

@overload
def ReactionToCXSmiles(reaction: ChemicalReaction, params: rdkit.Chem.rdmolfiles.SmilesWriteParams, flags: int = CXSmilesFields.CX_ALL) -> str:
    """construct a reaction CXSMILES string for a ChemicalReaction"""

def ReactionFromRxnFile(filename: str, sanitize: bool = False, removeHs: bool = False, strictParsing: bool = True) -> ChemicalReaction | None:
    """construct a ChemicalReaction from an MDL rxn file"""

def ReactionFromRxnBlock(rxnblock: str, sanitize: bool = False, removeHs: bool = False, strictParsing: bool = True) -> ChemicalReaction:
    """construct a ChemicalReaction from a string in MDL rxn format"""

def ReactionFromMrvFile(filename: str, sanitize: bool = False, removeHs: bool = False) -> ChemicalReaction | None:
    """construct a ChemicalReaction from an Marvin (mrv) rxn file"""

def ReactionFromMrvBlock(rxnblock: str | bytes, sanitize: bool = False, removeHs: bool = False) -> ChemicalReaction | None:
    """construct a ChemicalReaction from a string in Marvin (mrv) format"""

def MrvFileIsReaction(filename: str) -> bool:
    """returns whether or not an MRV file contains reaction data"""

def MrvBlockIsReaction(mrvData: str) -> bool:
    """returns whether or not an MRV block contains reaction data"""

def ReactionsFromCDXMLFile(filename: str, sanitize: bool = False, removeHs: bool = False) -> tuple[ChemicalReaction, ...]:
    """construct a tuple of ChemicalReactions from a CDXML rxn file"""

def ReactionsFromCDXMLBlock(rxnblock: str | bytes, sanitize: bool = False, removeHs: bool = False) -> tuple[ChemicalReaction, ...]:
    """construct a tuple of ChemicalReactions from a string in CDXML format"""

def ReactionToRxnBlock(reaction: ChemicalReaction, separateAgents: bool = False, forceV3000: bool = False) -> str:
    """construct a string in MDL rxn format for a ChemicalReaction"""

def ReactionToMrvBlock(reaction: ChemicalReaction, prettyPrint: bool = False) -> str:
    """construct a string in Marvin (MRV) rxn format for a ChemicalReaction"""

def ReactionToMrvFile(reaction: ChemicalReaction, filename: str, prettyPrint: bool = False) -> None:
    """write a Marvin (MRV) rxn file for a ChemicalReaction"""

def ReactionToV3KRxnBlock(reaction: ChemicalReaction, separateAgents: bool = False) -> str:
    """construct a string in MDL v3000 rxn format for a ChemicalReaction"""

def ReactionFromPNGFile(fname: str) -> ChemicalReaction:
    """construct a ChemicalReaction from metadata in a PNG file"""

def ReactionFromPNGString(data: bytes) -> ChemicalReaction:
    """construct a ChemicalReaction from an string with PNG data"""

def ReactionMetadataToPNGFile(mol: ChemicalReaction, filename: object, includePkl: bool = True, includeSmiles: bool = True, includeSmarts: bool = False, includeMol: bool = False) -> bytes:
    """
    Reads the contents of a PNG file and adds metadata about a reaction to
    it. The modified file contents are returned.
    """

def ReactionMetadataToPNGString(mol: ChemicalReaction, pngdata: bytes, includePkl: bool = True, includeSmiles: bool = True, includeSmarts: bool = False, includeRxn: bool = False) -> bytes:
    """
    Adds metadata about a reaction to the PNG string passed in.
    The modified string is returned.
    """

def ReactionFromMolecule(mol: rdkit.Chem.rdchem.Mol) -> ChemicalReaction:
    """
    construct a ChemicalReaction from an molecule if the RXN role property of the molecule is set
    """

def ReactionToMolecule(reaction: ChemicalReaction) -> rdkit.Chem.rdchem.Mol:
    """construct a molecule for a ChemicalReaction with RXN role property set"""

def Compute2DCoordsForReaction(reaction: ChemicalReaction, spacing: float = 1.0, updateProps: bool = True, canonOrient: bool = True, nFlipsPerSample: int = 0, nSample: int = 0, sampleSeed: int = 0, permuteDeg4Nodes: bool = False, bondLength: float = -1.0) -> None:
    """
    Compute 2D coordinates for a reaction.
      ARGUMENTS:
         - reaction - the reaction of interest
         - spacing - the amount of space left between components of the reaction
         - canonOrient - orient the reactants and products in a canonical way
         - updateProps - if set, properties such as conjugation and
            hybridization will be calculated for the reactant and product
            templates before generating coordinates. This should result in
            better depictions, but can lead to errors in some cases.
         - nFlipsPerSample - number of rotatable bonds that are
                    flipped at random at a time.
         - nSample - Number of random samplings of rotatable bonds.
         - sampleSeed - seed for the random sampling process.
         - permuteDeg4Nodes - allow permutation of bonds at a degree 4
                     node during the sampling process
         - bondLength - change the default bond length for depiction
    """

def CreateDifferenceFingerprintForReaction(reaction: ChemicalReaction, ReactionFingerPrintParams: ReactionFingerprintParams = ...) -> rdkit.DataStructs.cDataStructs.UIntSparseIntVect:
    """
    construct a difference fingerprint for a ChemicalReaction by
    subtracting the reactant fingerprint from the product fingerprint
    """

def CreateStructuralFingerprintForReaction(reaction: ChemicalReaction, ReactionFingerPrintParams: ReactionFingerprintParams = ...) -> rdkit.DataStructs.cDataStructs.ExplicitBitVect:
    """
    construct a structural fingerprint for a ChemicalReaction by
    concatenating the reactant fingerprint and the product fingerprint
    """

def IsReactionTemplateMoleculeAgent(molecule: rdkit.Chem.rdchem.Mol, agentThreshold: float) -> bool:
    """
    tests if a molecule can be classified as an agent depending on the ratio of mapped atoms and a give threshold
    """

def HasReactionAtomMapping(rxn: ChemicalReaction) -> bool:
    """tests if a reaction obtains any atom mapping"""

def HasReactionSubstructMatch(reaction: ChemicalReaction, queryReaction: ChemicalReaction, includeAgents: bool = False) -> bool:
    """tests if the queryReaction is a substructure of a reaction"""

def HasAgentTemplateSubstructMatch(reaction: ChemicalReaction, queryReaction: ChemicalReaction) -> bool:
    """
    tests if the agents of a queryReaction are the same as those of a reaction
    """

def HasProductTemplateSubstructMatch(reaction: ChemicalReaction, queryReaction: ChemicalReaction) -> bool:
    """
    tests if the products of a queryReaction are substructures of the products of a reaction
    """

def HasReactantTemplateSubstructMatch(reaction: ChemicalReaction, queryReaction: ChemicalReaction) -> bool:
    """
    tests if the reactants of a queryReaction are substructures of the reactants of a reaction
    """

def UpdateProductsStereochemistry(reaction: ChemicalReaction) -> None:
    """
    Caution: This is an expert-user function which will change a property (molInversionFlag) of your products.
              This function is called by default using the RXN or SMARTS parser for reactions and should really only be called if reactions have been constructed some other way.
              The function updates the stereochemistry of the product by considering 4 different cases: inversion, retention, removal, and introduction
    """

def ReduceProductToSideChains(product: rdkit.Chem.rdchem.Mol, addDummyAtoms: bool = True) -> rdkit.Chem.rdchem.Mol:
    """
    reduce the product of a reaction to the side chains added by the reaction.
                  The output is a molecule with attached wildcards indicating where the product was attached.
                  The dummy atom has the same reaction-map number as the product atom (if available).
    """

def RemoveMappingNumbersFromReactions(reaction: ChemicalReaction) -> None:
    """Removes the mapping numbers from the molecules of a reaction"""

def PreprocessReaction(reaction: ChemicalReaction, queries: dict[str, rdkit.Chem.rdchem.Mol] = {}, propName: str = 'molFileValue') -> tuple[int, int, int, int, tuple[tuple[tuple[int, str], ...], ...]]:
    """
    A function for preprocessing reactions with more specific queries.
    Queries are indicated by labels on atoms (molFileAlias property by default)
    When these labels are found, more specific queries are placed on the atoms.
    By default, the available quieries come from
      FilterCatalog.GetFlattenedFunctionalGroupHierarchy(True)n
    Sample Usage:
      >>> from rdkit import Chem, RDConfig
      >>> from rdkit.Chem import MolFromSmiles, AllChem
      >>> from rdkit.Chem.rdChemReactions import PreprocessReaction
      >>> import os
      >>> testFile = os.path.join(RDConfig.RDCodeDir,'Chem','SimpleEnum','test_data','boronic1.rxn')
      >>> rxn = AllChem.ReactionFromRxnFile(testFile)
      >>> rxn.Initialize()
      >>> nWarn,nError,nReacts,nProds,reactantLabels = PreprocessReaction(rxn)
      >>> nWarn
      0
      >>> nError
      0
      >>> nReacts
      2
      >>> nProds
      1
      >>> reactantLabels
      (((0, 'halogen.bromine.aromatic'),), ((1, 'boronicacid'),))

    If there are functional group labels in the input reaction (via atoms with molFileValue properties),
    the corresponding atoms will have queries added to them so that they only match such things. We can
    see this here:
      >>> rxn = AllChem.ReactionFromRxnFile(testFile)
      >>> rxn.Initialize()
      >>> r1 = rxn.GetReactantTemplate(0)
      >>> m1 = Chem.MolFromSmiles('CCBr')
      >>> m2 = Chem.MolFromSmiles('c1ccccc1Br')

    These both match because the reaction file itself just has R1-Br:
      >>> m1.HasSubstructMatch(r1)
      True
      >>> m2.HasSubstructMatch(r1)
      True

    After preprocessing, we only match the aromatic Br:
      >>> d = PreprocessReaction(rxn)
      >>> m1.HasSubstructMatch(r1)
      False
      >>> m2.HasSubstructMatch(r1)
      True

    We also support or queries in the values field (separated by commas):
      >>> testFile = os.path.join(RDConfig.RDCodeDir,'Chem','SimpleEnum','test_data','azide_reaction.rxn')
      >>> rxn = AllChem.ReactionFromRxnFile(testFile)
      >>> rxn.Initialize()
      >>> reactantLabels = PreprocessReaction(rxn)[-1]
      >>> reactantLabels
      (((1, 'azide'),), ((1, 'carboxylicacid,acidchloride'),))
      >>> m1 = Chem.MolFromSmiles('CC(=O)O')
      >>> m2 = Chem.MolFromSmiles('CC(=O)Cl')
      >>> m3 = Chem.MolFromSmiles('CC(=O)N')
      >>> r2 = rxn.GetReactantTemplate(1)
      >>> m1.HasSubstructMatch(r2)
      True
      >>> m2.HasSubstructMatch(r2)
      True
      >>> m3.HasSubstructMatch(r2)
      False

    unrecognized final group types are returned as None:
      >>> testFile = os.path.join(RDConfig.RDCodeDir,'Chem','SimpleEnum','test_data','bad_value1.rxn')
      >>> rxn = AllChem.ReactionFromRxnFile(testFile)
      >>> rxn.Initialize()
      >>> nWarn,nError,nReacts,nProds,reactantLabels = PreprocessReaction(rxn)
      Traceback (most recent call last):
        ...
      KeyError: 'boromicacid'

    One unrecognized group type in a comma-separated list makes the whole thing fail:
      >>> testFile = os.path.join(RDConfig.RDCodeDir,'Chem','SimpleEnum','test_data','bad_value2.rxn')
      >>> rxn = AllChem.ReactionFromRxnFile(testFile)
      >>> rxn.Initialize()
      >>> nWarn,nError,nReacts,nProds,reactantLabels = PreprocessReaction(rxn)
      Traceback (most recent call last):
        ...
      KeyError: 'carboxylicacid,acidchlroide'
      >>> testFile = os.path.join(RDConfig.RDCodeDir,'Chem','SimpleEnum','test_data','bad_value3.rxn')
      >>> rxn = AllChem.ReactionFromRxnFile(testFile)
      >>> rxn.Initialize()
      >>> nWarn,nError,nReacts,nProds,reactantLabels = PreprocessReaction(rxn)
      Traceback (most recent call last):
        ...
      KeyError: 'carboxyliccaid,acidchloride'
      >>> rxn = rdChemReactions.ChemicalReaction()
      >>> rxn.Initialize()
      >>> nWarn,nError,nReacts,nProds,reactantLabels = PreprocessReaction(rxn)
      >>> reactantLabels
      ()
      >>> reactantLabels == ()
      True
    """

class SanitizeFlags(enum.Enum):
    SANITIZE_NONE = 0

    SANITIZE_ATOM_MAPS = 2

    SANITIZE_RGROUP_NAMES = 1

    SANITIZE_ADJUST_REACTANTS = 4

    SANITIZE_MERGEHS = 8

    SANITIZE_ALL = 4294967295

SANITIZE_NONE: SanitizeFlags = SanitizeFlags.SANITIZE_NONE

SANITIZE_ATOM_MAPS: SanitizeFlags = SanitizeFlags.SANITIZE_ATOM_MAPS

SANITIZE_RGROUP_NAMES: SanitizeFlags = SanitizeFlags.SANITIZE_RGROUP_NAMES

SANITIZE_ADJUST_REACTANTS: SanitizeFlags = SanitizeFlags.SANITIZE_ADJUST_REACTANTS

SANITIZE_MERGEHS: SanitizeFlags = SanitizeFlags.SANITIZE_MERGEHS

SANITIZE_ALL: SanitizeFlags = SanitizeFlags.SANITIZE_ALL

def GetDefaultAdjustParams() -> rdkit.Chem.rdmolops.AdjustQueryParameters:
    """Returns the default adjustment parameters for reactant templates"""

def GetChemDrawRxnAdjustParams() -> rdkit.Chem.rdmolops.AdjustQueryParameters:
    """
    (deprecated, see MatchOnlyAtRgroupsAdjustParams)
    	Returns the chemdraw style adjustment parameters for reactant templates
    """

def MatchOnlyAtRgroupsAdjustParams() -> rdkit.Chem.rdmolops.AdjustQueryParameters:
    """Only match at the specified rgroup locations in the reactant templates"""

def SanitizeRxn(rxn: ChemicalReaction, sanitizeOps: int = 4294967295, params: rdkit.Chem.rdmolops.AdjustQueryParameters | None = None, catchErrors: bool = False) -> SanitizeFlags:
    """
    Does some sanitization of the reactant and product templates of a reaction.

        - The reaction is modified in place.
        - If sanitization fails, an exception will be thrown unless catchErrors is set

      ARGUMENTS:

        - rxn: the reaction to be modified
        - sanitizeOps: (optional) reaction sanitization operations to be carried out
          these should be constructed by or'ing together the
          operations in rdkit.Chem.rdChemReactions.SanitizeFlags
        - optional adjustment parameters for changing the meaning of the substructure
          matching done in the templates.  The default is
          rdkit.Chem.rdChemReactions.DefaultRxnAdjustParams which aromatizes
          kekule structures if possible.
        - catchErrors: (optional) if provided, instead of raising an exception
          when sanitization fails (the default behavior), the
          first operation that failed (as defined in rdkit.Chem.rdChemReactions.SanitizeFlags)
          is returned. Zero is returned on success.

      The operations carried out by default are:
        1) fixRGroups(): sets R group labels on mapped dummy atoms when possible
        2) fixAtomMaps(): attempts to set atom maps on unmapped R groups
        3) adjustTemplate(): calls adjustQueryProperties() on all reactant templates
        4) fixHs(): merges explicit Hs in the reactant templates that don't map to heavy atoms
    """

def SanitizeRxnAsMols(rxn: ChemicalReaction, sanitizeOps: int = 268435455) -> None:
    """
    Does the usual molecular sanitization on each reactant, agent, and product of the reaction
    """

class EnumerateLibraryBase:
    def __bool__(self) -> bool: ...

    def __iter__(self) -> object: ...

    def next(self) -> tuple:
        """Return the next molecule from the enumeration."""

    def __next__(self) -> tuple:
        """Return the next molecule from the enumeration."""

    def nextSmiles(self) -> list[list[str]]:
        """Return the next smiles string from the enumeration."""

    def Serialize(self) -> bytes:
        """
        Serialize the library to a binary string.
        Note that the position in the library is serialized as well.  Care should
        be taken when serializing.  See GetState/SetState for position manipulation.
        """

    @overload
    def InitFromString(self, data: str) -> None: ...

    @overload
    def InitFromString(self, data: bytes) -> None:
        """Initialize the library from a binary string"""

    def GetPosition(self) -> list[int]:
        """
        Returns the current enumeration position into the reagent vectors, as
        returned by GetReagents().  They do not necessarily refer to
        the input reagent sets as they only refer to reagents compatible
        with the reaction.
        """

    def GetState(self) -> str:
        """
        Returns the current enumeration state (position) of the library.
        This position can be used to restart the library from a known position
        """

    def SetState(self, state: str) -> None:
        """Sets the enumeration state (position) of the library."""

    def ResetState(self) -> None:
        """
        Returns the current enumeration state (position) of the library to the start.
        """

    def GetReaction(self) -> ChemicalReaction:
        """Returns the chemical reaction for this library"""

    def GetEnumerator(self) -> EnumerationStrategyBase:
        """Returns the enumation strategy for the current library"""

class EnumerationParams:
    """
    EnumerationParams
    Controls some aspects of how the enumeration is performed.
    Options:
      reagentMaxMatchCount [ default Infinite ]
        This specifies how many times the reactant template can match a reagent.

      sanePartialProducts [default false]
        If true, forces all products of the reagent plus the product templates
         pass chemical sanitization.  Note that if the product template itself
         does not pass sanitization, then none of the products will.
    """

    def __init__(self) -> None: ...

    @property
    def reagentMaxMatchCount(self) -> int: ...

    @reagentMaxMatchCount.setter
    def reagentMaxMatchCount(self, arg: int, /) -> None: ...

    @property
    def sanePartialProducts(self) -> bool: ...

    @sanePartialProducts.setter
    def sanePartialProducts(self, arg: bool, /) -> None: ...

class EnumerateLibrary(EnumerateLibraryBase):
    """
    EnumerateLibrary
    This class allows easy enumeration of reactions.  Simply provide a reaction
    and a set of reagents and you are off the races.

    Note that this functionality should be considered beta and that the API may
    change in a future release.

    EnumerateLibrary follows the python enumerator protocol, for example:

    library = EnumerateLibrary(rxn, bbs)
    for products in library:
       ... do something with the product

    It is useful to sanitize reactions before hand:

    SanitizeRxn(rxn)
    library = EnumerateLibrary(rxn, bbs)

    If ChemDraw style reaction semantics are prefereed, you can apply
    the ChemDraw parameters:

    SanitizeRxn(rxn, params=GetChemDrawRxnAdjustParams())

    For one, this enforces only matching RGroups and assumes all atoms
    have fully satisfied valences.

    Each product has the same output as applying a set of reagents to
    the libraries reaction.

    This can be a bit confusing as each product can have multiple molecules
    generated.  The returned data structure is as follows:

       [ [products1], [products2],... ]
    Where products1 are the molecule products for the reactions first product
    template and products2 are the molecule products for the second product
    template.  Since each reactant can match more than once, there may be
    multiple product molecules for each template.

    for products in library:
        for results_for_product_template in products:
            for mol in results_for_product_template:
                Chem.MolToSmiles(mol) # finally have a molecule!

    For sufficiently large libraries, using this iteration strategy is not
    recommended as the library may contain more products than atoms in the
    universe.  To help with this, you can supply an enumeration strategy.
    The default strategy is a CartesianProductStrategy which enumerates
    everything.  RandomSampleStrategy randomly samples the products but
    this strategy never terminates, however, python supplies itertools:

    import itertools
    library = EnumerateLibrary(rxn, bbs, rdChemReactions.RandomSampleStrategy())
    for result in itertools.islice(library, 1000):
        # do something with the first 1000 samples

    for result in itertools.islice(library, 1000):
        # do something with the next 1000 samples

    Libraries are also serializable, including their current state:

    s = library.Serialize()
    library2 = EnumerateLibrary()
    library2.InitFromString(s)
    for result in itertools.islice(libary2, 1000):
        # do something with the next 1000 samples
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, rxn: ChemicalReaction, reagents: Iterable[Iterable[rdkit.Chem.rdchem.Mol]], params: EnumerationParams = ...) -> None: ...

    @overload
    def __init__(self, rxn: ChemicalReaction, reagents: Iterable[Iterable[rdkit.Chem.rdchem.Mol]], enumerator: EnumerationStrategyBase, params: EnumerationParams = ...) -> None: ...

    def GetReagents(self) -> list:
        """
        Return the reagents used in this library.  These are the subset
        of the input reagents that are compatible with the reaction so may
        be smaller than the input reagent sets.
        """

class EnumerationStrategyBase:
    def __bool__(self) -> bool: ...

    def Type(self) -> str:
        """Returns the enumeration strategy type as a string."""

    def Skip(self, skipCount: int) -> bool:
        """
        Skip the next Nth results. note: this may be an expensive operation
        depending on the enumeration strategy used. It is recommended to use
        the enumerator state to advance to a known position
        """

    def __copy__(self) -> EnumerationStrategyBase: ...

    def GetNumPermutations(self) -> int:
        """
        Returns the total number of results for this enumeration strategy.
        Note that some strategies are effectively infinite.
        """

    def GetPosition(self) -> list[int]:
        """
        Return the current indices into the arrays of reagents, as
        returned by GetReagents().  They do not necessarily refer to
        the input reagent sets as they only refer to reagents compatible
        with the reaction.
        """

    def next(self) -> list[int]:
        """Return the next indices into the arrays of reagents"""

    def __next__(self) -> list[int]:
        """Return the next indices into the arrays of reagents"""

    def Initialize(self, rxn: ChemicalReaction, ob: Iterable[Iterable[rdkit.Chem.rdchem.Mol]]) -> None: ...

class CartesianProductStrategy(EnumerationStrategyBase):
    """
    CartesianProductStrategy produces a standard walk through all possible
    reagent combinations:

    (0,0,0), (1,0,0), (2,0,0) ...
    """

    def __init__(self) -> None: ...

    def __copy__(self) -> EnumerationStrategyBase: ...

class RandomSampleStrategy(EnumerationStrategyBase):
    """
    RandomSampleStrategy simply randomly samples from the reagent sets.
    Note that this strategy never halts and can produce duplicates.
    """

    def __init__(self) -> None: ...

    def __copy__(self) -> EnumerationStrategyBase: ...

class RandomSampleAllBBsStrategy(EnumerationStrategyBase):
    """
    RandomSampleAllBBsStrategy randomly samples from the reagent sets
    with the constraint that all building blocks are samples as early as possible.
    Note that this strategy never halts and can produce duplicates.
    """

    def __init__(self) -> None: ...

    def __copy__(self) -> EnumerationStrategyBase: ...

class EvenSamplePairsStrategy(EnumerationStrategyBase):
    """
    Randomly sample Pairs evenly from a collection of building blocks
    This is a good strategy for choosing a relatively small selection
    of building blocks from a larger set.  As the amount of work needed
    to retrieve the next evenly sample building block grows with the
    number of samples, this method performs progressively worse as the
    number of samples gets larger.
    See EnumerationStrategyBase for more details.
    """

    def __init__(self) -> None: ...

    def __copy__(self) -> EnumerationStrategyBase: ...

    def Stats(self) -> str:
        """
        Return the statistics log of the pairs used in the current enumeration.
        """

def EnumerateLibraryCanSerialize() -> bool:
    """
    Returns True if the EnumerateLibrary is serializable (requires boost serialization)
    """
