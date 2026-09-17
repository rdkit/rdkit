"""
Module containing RDKit functionality for working with molecular file formats.
"""

from collections.abc import Iterable, Mapping, Sequence
import enum
import os
from typing import overload

import rdkit.Chem.rdchem


class BadFileException(OSError):
    pass

class FileParseException(RuntimeError):
    pass

def MolFromTPLFile(fileName: str | os.PathLike, sanitize: bool = True, skipFirstConf: bool = False) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a molecule from a TPL file.

        ARGUMENTS:

          - fileName: name of the file to read

          - sanitize: (optional) toggles sanitization of the molecule.
            Defaults to True.

          - skipFirstConf: (optional) skips reading the first conformer.
            Defaults to False.
            This should be set to True when reading TPLs written by 
            the CombiCode.

        RETURNS:

          a Mol object, None on failure.
    """

def MolFromTPLBlock(tplBlock: str | bytes, sanitize: bool = True, skipFirstConf: bool = False) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a molecule from a TPL block.

        ARGUMENTS:

          - tplBlock: string containing the TPL block

          - sanitize: (optional) toggles sanitization of the molecule.
            Defaults to True.

          - skipFirstConf: (optional) skips reading the first conformer.
            Defaults to False.
            This should be set to True when reading TPLs written by 
            the CombiCode.

        RETURNS:

          a Mol object, None on failure.
    """

def MolFromMrvFile(molFileName: str | os.PathLike, sanitize: bool = True, removeHs: bool = True) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a molecule from a Marvin (Mrv) file.

        ARGUMENTS:

          - fileName: name of the file to read

          - sanitize: (optional) toggles sanitization of the molecule.
            Defaults to true.

          - removeHs: (optional) toggles removing hydrogens from the molecule.
            This only make sense when sanitization is done.
            Defaults to true.

        RETURNS:

          a Mol object, None on failure.
    """

def MolFromMrvBlock(mrvBlock: str | bytes, sanitize: bool = True, removeHs: bool = True) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a molecule from a Marvin (mrv) block.

        ARGUMENTS:

          - molBlock: string containing the Marvin block

          - sanitize: (optional) toggles sanitization of the molecule.
            Defaults to True.

          - removeHs: (optional) toggles removing hydrogens from the molecule.
            This only make sense when sanitization is done.
            Defaults to true.

        RETURNS:

          a Mol object, None on failure.
    """

def MolFromXYZFile(xyzFileName: str) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a molecule from an XYZ file.
        ARGUMENTS:

          - xyzFileName: name of the file to read

        RETURNS:

          a Mol object, None on failure.
    """

def MolFromXYZBlock(xyzBlock: str | bytes) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a molecule from an XYZ string.
        ARGUMENTS:

          - xyzBlock: the XYZ data to read

        RETURNS:

          a Mol object, None on failure.
    """

def MolFromRDKitSVG(svg: str | bytes, sanitize: bool = True, removeHs: bool = True) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a molecule from an RDKit-generate SVG string.
        ARGUMENTS:

          - svg: string containing the SVG data (must include molecule
          metadata)

          - sanitize: (optional) toggles sanitization of the molecule.
            Defaults to True.

          - removeHs: (optional) toggles removing hydrogens from the molecule.
            This only make sense when sanitization is done.
            Defaults to true.

        RETURNS:

          a Mol object, None on failure.

        NOTE: this functionality should be considered beta.
    """

def MolFromMol2File(mol2FileName: str | os.PathLike, sanitize: bool = True, removeHs: bool = True, cleanupSubstructures: bool = True) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a molecule from a Tripos Mol2 file.
        NOTE:
          The parser expects the atom-typing scheme used by Corina.
          Atom types from Tripos' dbtranslate are less supported.
          Other atom typing schemes are unlikely to work.

        ARGUMENTS:                                  \\

          - mol2FileName: name of the file to read

          - sanitize: (optional) toggles sanitization of the molecule.
            Defaults to true.

          - removeHs: (optional) toggles removing hydrogens from the molecule.
            This only make sense when sanitization is done.
            Defaults to true.

          - cleanupSubstructures: (optional) toggles standardizing some 
            substructures found in mol2 files.
            Defaults to true.

        RETURNS:

          a Mol object, None on failure.
    """

def MolFromMol2Block(mol2Block: str, sanitize: bool = True, removeHs: bool = True, cleanupSubstructures: bool = True) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a molecule from a Tripos Mol2 block.
        NOTE:
          The parser expects the atom-typing scheme used by Corina.
          Atom types from Tripos' dbtranslate are less supported.
          Other atom typing schemes are unlikely to work.

        ARGUMENTS:

          - mol2Block: string containing the Mol2 block

          - sanitize: (optional) toggles sanitization of the molecule.
            Defaults to True.

          - removeHs: (optional) toggles removing hydrogens from the molecule.
            This only make sense when sanitization is done.
            Defaults to true.

          - cleanupSubstructures: (optional) toggles standardizing some
            substructures found in mol2 files.
            Defaults to true.

        RETURNS:

          a Mol object, None on failure.
    """

def MolFromMolFile(molFileName: str | os.PathLike, sanitize: bool = True, removeHs: bool = True, strictParsing: bool = True) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a molecule from a Mol file.
        ARGUMENTS:

          - molFileName: name of the file to read

          - sanitize: (optional) toggles sanitization of the molecule.
            Defaults to true.

          - removeHs: (optional) toggles removing hydrogens from the molecule.
            This only make sense when sanitization is done.
            Defaults to true.

          - strictParsing: (optional) if this is false, the parser is more lax
          about.
            correctness of the content.
            Defaults to true.

        RETURNS:

          a Mol object, None on failure.
    """

def MolFromMolBlock(molBlock: str | bytes, sanitize: bool = True, removeHs: bool = True, strictParsing: bool = True) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a molecule from a Mol block.
        ARGUMENTS:

          - molBlock: string containing the Mol block

          - sanitize: (optional) toggles sanitization of the molecule.
            Defaults to True.

          - removeHs: (optional) toggles removing hydrogens from the molecule.
            This only make sense when sanitization is done.
            Defaults to true.

          - strictParsing: (optional) if this is false, the parser is more lax
          about.
            correctness of the content.
            Defaults to true.

        RETURNS:

          a Mol object, None on failure.
    """

class MolWriterParams:
    """Parameters controlling Mol writing"""

    @property
    def includeStereo(self) -> bool:
        """toggles inclusion of stereochemistry information (default=True)"""

    @includeStereo.setter
    def includeStereo(self, arg: bool, /) -> None: ...

    @property
    def kekulize(self) -> bool:
        """
        triggers kekulization of the molecule before it is written (default=True)
        """

    @kekulize.setter
    def kekulize(self, arg: bool, /) -> None: ...

    @property
    def forceV3000(self) -> bool:
        """
        force generation a V3000 mol block (happens automatically with more than 999 atoms or bonds)(default=False)
        """

    @forceV3000.setter
    def forceV3000(self, arg: bool, /) -> None: ...

    @property
    def precision(self) -> int:
        """precision of coordinates (only available in V3000)(default=false)"""

    @precision.setter
    def precision(self, arg: int, /) -> None: ...

    def __setattr__(self, name: str, value: object | None) -> None: ...

class MolFromSCSRParams:
    """Parameters controlling conversion of an SCSRMol to a Mol"""

    def __init__(self) -> None: ...

    @property
    def includeLeavingGroups(self) -> bool:
        """include leaving groups atoms if not substited at that position"""

    @includeLeavingGroups.setter
    def includeLeavingGroups(self, arg: bool, /) -> None: ...

    @property
    def scsrTemplateNames(self) -> SCSRTemplateNames:
        """
        If True, the first template name in the Sgroup is used as the Sgroup label
        """

    @scsrTemplateNames.setter
    def scsrTemplateNames(self, arg: SCSRTemplateNames, /) -> None: ...

    @property
    def scsrBaseHbondOptions(self) -> SCSRBaseHbondOptions:
        """One of Ignore, UseSapAll(default) , UseSapOne, Auto"""

    @scsrBaseHbondOptions.setter
    def scsrBaseHbondOptions(self, arg: SCSRBaseHbondOptions, /) -> None: ...

    def __setattr__(self, name: str, value: object | None) -> None: ...

def MolFromSCSRBlock(molBlock: str, sanitize: bool = True, removeHs: bool = True, molFromSCSRParams: MolFromSCSRParams | None = None) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a molecule from an SCSR Mol block.
            ARGUMENTS:

              - molBlock: string containing the SCSR Mol block

              - sanitize: (optional) toggles sanitization of the molecule.
                Defaults to True.

              - removeHs: (optional) toggles removing hydrogens from the
              molecule.
                This only make sense when sanitization is done.
                Defaults to true.

              - molFromSCSRParams : MolFromSCSRParams to control conversion
           RETURNS :
           a Mol object, None on failure.
    """

def MolFromSCSRFile(filename: str | os.PathLike, sanitize: bool = True, removeHs: bool = True, molFromSCSRParams: MolFromSCSRParams | None = None) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a molecule from an SCSR Mol block.
            ARGUMENTS:

              - filename: string containing the SCSR filename

              - sanitize: (optional) toggles sanitization of the molecule.
                Defaults to True.

              - removeHs: (optional) toggles removing hydrogens from the
              molecule.
                This only make sense when sanitization is done.
                Defaults to true.

              - molFromSCSRParams : MolFromSCSRParams to control conversion
           RETURNS :
           a Mol object, None on failure.
    """

@overload
def MolToMolBlock(mol: rdkit.Chem.rdchem.Mol, params: MolWriterParams, confId: int = -1) -> str:
    """
    Returns a Mol block for a molecule
        Arguments:
          - mol: the molecule
          - params: the MolWriterParams
          - confId: (optional) selects which conformation to output (-1 =
          default)

        RETURNS:

          a string
    """

@overload
def MolToMolBlock(mol: rdkit.Chem.rdchem.Mol, includeStereo: bool = True, confId: int = -1, kekulize: bool = True, forceV3000: bool = False) -> str:
    """
    Returns a Mol block for a molecule
        ARGUMENTS:

          - mol: the molecule
          - includeStereo: (optional) toggles inclusion of stereochemical
            information in the output
          - confId: (optional) selects which conformation to output (-1 =
          default)
          - kekulize: (optional) triggers kekulization of the molecule before
          it's written,
            as suggested by the MDL spec.
          - forceV3000 (optional) force generation a V3000 mol block (happens
          automatically with
            more than 999 atoms or bonds)

        RETURNS:

          a string
    """

@overload
def MolToV3KMolBlock(mol: rdkit.Chem.rdchem.Mol, params: MolWriterParams, confId: int = -1) -> str:
    """
    Returns a V3000 Mol block for a molecule
         ARGUMENTS:
       \\
           - mol: the molecule
           - params: the MolWriterParams
           - confId: (optional) selects which conformation to output (-1 =
           default)
       \\
         RETURNS:
       \\
           a string
    """

@overload
def MolToV3KMolBlock(mol: rdkit.Chem.rdchem.Mol, includeStereo: bool = True, confId: int = -1, kekulize: bool = True) -> str:
    """
    Returns a V3000 Mol block for a molecule
         ARGUMENTS:
       \\
           - mol: the molecule
           - includeStereo: (optional) toggles inclusion of stereochemical
             information in the output
           - confId: (optional) selects which conformation to output (-1 =
           default)
           - kekulize: (optional) triggers kekulization of the molecule before
           it's written,
             as suggested by the MDL spec.
       \\
         RETURNS:
       \\
           a string
    """

def MolToV2KMolBlock(mol: rdkit.Chem.rdchem.Mol, params: MolWriterParams | None = None, confId: int = -1) -> str:
    """
    Returns a V2000 Mol block for a molecule
         ARGUMENTS:

           - mol: the molecule
           - params: the MolWriterParams
           - confId: (optional) selects which conformation to output (-1 =
           default)

         RETURNS:

           a string

         NOTE: this function throws a ValueError if the molecule has more than
         999 atoms, bonds, or SGroups
    """

@overload
def MolToMolFile(mol: rdkit.Chem.rdchem.Mol, filename: str, params: MolWriterParams, confId: int = -1) -> None:
    """
    Writes a Mol file for a molecule
        ARGUMENTS:

          - mol: the molecule
          - filename: the file to write to
          - params: the MolWriterParams
          - confId: (optional) selects which conformation to output (-1 =
          default)
    """

@overload
def MolToMolFile(mol: rdkit.Chem.rdchem.Mol, filename: str, includeStereo: bool = True, confId: int = -1, kekulize: bool = True, forceV3000: bool = False) -> None:
    """
    Writes a Mol file for a molecule
        ARGUMENTS:

          - mol: the molecule
          - filename: the file to write to
          - includeStereo: (optional) toggles inclusion of stereochemical
            information in the output
          - confId: (optional) selects which conformation to output (-1 =
          default)
          - kekulize: (optional) triggers kekulization of the molecule before
          it's written,
            as suggested by the MDL spec.
          - forceV3000 (optional) force generation a V3000 mol block (happens
          automatically with
            more than 999 atoms or bonds)
    """

@overload
def MolToV3KMolFile(mol: rdkit.Chem.rdchem.Mol, filename: str, params: MolWriterParams = True, confId: int = -1) -> None:
    """
    Writes a V3000 Mol file for a molecule
        ARGUMENTS:

          - mol: the molecule
          - filename: the file to write to
          - params: the MolWriterParams
          - confId: (optional) selects which conformation to output (-1 =
          default)
    """

@overload
def MolToV3KMolFile(mol: rdkit.Chem.rdchem.Mol, filename: str, includeStereo: bool = True, confId: int = -1, kekulize: bool = True) -> None:
    """
    Writes a V3000 Mol file for a molecule
        ARGUMENTS:

          - mol: the molecule
          - filename: the file to write to
          - includeStereo: (optional) toggles inclusion of stereochemical
            information in the output
          - confId: (optional) selects which conformation to output (-1 =
          default)
          - kekulize: (optional) triggers kekulization of the molecule before
          it's written,
            as suggested by the MDL spec.
    """

@overload
def MolToMrvBlock(mol: rdkit.Chem.rdchem.Mol, includeStereo: bool = True, confId: int = -1, kekulize: bool = True, prettyPrint: bool = False) -> str:
    """
    Returns a Marvin (Mrv) Mol block for a molecule
        ARGUMENTS:

          - mol: the molecule
          - includeStereo: (optional) toggles inclusion of stereochemical
            information in the output
          - confId: (optional) selects which conformation to output (-1 =
          default)
          - kekulize: (optional) triggers kekulization of the molecule before
          it's written.
          - prettyPrint: (optional) makes the output more human readable.

        RETURNS:

          a string
    """

@overload
def MolToMrvBlock(mol: rdkit.Chem.rdchem.Mol, params: MrvWriterParams, confId: int = -1) -> str:
    """
    Returns a Marvin (Mrv) Mol block for a molecule
        ARGUMENTS:

          - mol: the molecule
          - params: marvin write params
          - confId: (optional) selects which conformation to output (-1 =
          default)

        RETURNS:

          a string
    """

class MrvWriterParams:
    """Parameters controlling Marvin (Mrv) writing"""

    def __init__(self) -> None: ...

    @property
    def includeStereo(self) -> bool:
        """include stereochemical information"""

    @includeStereo.setter
    def includeStereo(self, arg: bool, /) -> None: ...

    @property
    def kekulize(self) -> bool:
        """kekulize the molecule before it is written"""

    @kekulize.setter
    def kekulize(self, arg: bool, /) -> None: ...

    @property
    def prettyPrint(self) -> bool:
        """make the output more human readable"""

    @prettyPrint.setter
    def prettyPrint(self, arg: bool, /) -> None: ...

    @property
    def precision(self) -> int:
        """number of significant digits in the coordinates"""

    @precision.setter
    def precision(self, arg: int, /) -> None: ...

@overload
def MolToMrvFile(mol: rdkit.Chem.rdchem.Mol, filename: str, includeStereo: bool = True, confId: int = -1, kekulize: bool = True, prettyPrint: bool = False) -> None:
    """
    Writes a Marvin (MRV) file for a molecule
         ARGUMENTS:

           - mol: the molecule
           - filename: the file to write to
           - includeStereo: (optional) toggles inclusion of stereochemical
             information in the output
           - confId: (optional) selects which conformation to output (-1 =
           default)
           - kekulize: (optional) triggers kekulization of the molecule before
           it's written.
           - prettyPrint: (optional) makes the output more human readable.
    """

@overload
def MolToMrvFile(mol: rdkit.Chem.rdchem.Mol, filename: str, params: MrvWriterParams, confId: int = -1) -> None:
    """
    Writes a Marvin (MRV) file for a molecule
        ARGUMENTS:

          - mol: the molecule
          - filename: the file to write to
          - params: marvin write params
          - confId: (optional) selects which conformation to output (-1 =
          default)
    """

def MolToCMLBlock(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, kekulize: bool = True) -> str:
    """
    Writes a CML block for a molecule
        ARGUMENTS:

          - mol: the molecule
          - confId: (optional) selects which conformation to output
          - kekulize: (optional) triggers kekulization of the molecule before
          it's written
    """

def MolToCMLFile(mol: rdkit.Chem.rdchem.Mol, filename: str, confId: int = -1, kekulize: bool = True) -> None:
    """
    Writes a CML file for a molecule
        ARGUMENTS:

          - mol: the molecule
          - filename: the file to write to
          - confId: (optional) selects which conformation to output
          - kekulize: (optional) triggers kekulization of the molecule before
          it's written
    """

def MolToXYZBlock(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, precision: int = 6) -> str:
    """
    Returns a XYZ block for a molecule
        ARGUMENTS:

          - mol: the molecule
          - confId: (optional) selects which conformation to output (-1 =
          default)
          - precision: precision of the coordinates

        RETURNS:

          a string
    """

def MolToXYZFile(mol: rdkit.Chem.rdchem.Mol, filename: str, confId: int = -1, precision: int = 6) -> None:
    """
    Writes a XYZ file for a molecule
        ARGUMENTS:

          - mol: the molecule
          - filename: the file to write to
          - confId: (optional) selects which conformation to output (-1 =
          default)
          - precision: precision of the coordinates
    """

class SmilesParserParams:
    """Parameters controlling SMILES parsing"""

    def __init__(self) -> None: ...

    @property
    def debugParse(self) -> int:
        """controls the amount of debugging information produced"""

    @debugParse.setter
    def debugParse(self, arg: int, /) -> None: ...

    @property
    def parseName(self) -> bool:
        """controls whether or not the molecule name is also parsed"""

    @parseName.setter
    def parseName(self, arg: bool, /) -> None: ...

    @property
    def allowCXSMILES(self) -> bool:
        """controls whether or not the CXSMILES extensions are parsed"""

    @allowCXSMILES.setter
    def allowCXSMILES(self, arg: bool, /) -> None: ...

    @property
    def strictCXSMILES(self) -> bool:
        """
        controls whether or not problems in CXSMILES parsing causes molecule parsing to fail
        """

    @strictCXSMILES.setter
    def strictCXSMILES(self, arg: bool, /) -> None: ...

    @property
    def sanitize(self) -> bool:
        """
        controls whether or not the molecule is sanitized before being returned
        """

    @sanitize.setter
    def sanitize(self, arg: bool, /) -> None: ...

    @property
    def removeHs(self) -> bool:
        """controls whether or not Hs are removed before the molecule is returned"""

    @removeHs.setter
    def removeHs(self, arg: bool, /) -> None: ...

    def __setattr__(self, name: str, value: object | None) -> None: ...

class SmartsParserParams:
    """Parameters controlling SMARTS parsing"""

    def __init__(self) -> None: ...

    @property
    def debugParse(self) -> int:
        """controls the amount of debugging information produced"""

    @debugParse.setter
    def debugParse(self, arg: int, /) -> None: ...

    @property
    def parseName(self) -> bool:
        """controls whether or not the molecule name is also parsed"""

    @parseName.setter
    def parseName(self, arg: bool, /) -> None: ...

    @property
    def allowCXSMILES(self) -> bool:
        """controls whether or not the CXSMILES extensions are parsed"""

    @allowCXSMILES.setter
    def allowCXSMILES(self, arg: bool, /) -> None: ...

    @property
    def strictCXSMILES(self) -> bool:
        """
        controls whether or not problems in CXSMILES parsing causes molecule parsing to fail
        """

    @strictCXSMILES.setter
    def strictCXSMILES(self, arg: bool, /) -> None: ...

    @property
    def mergeHs(self) -> bool:
        """toggles merging H atoms in the SMARTS into neighboring atoms"""

    @mergeHs.setter
    def mergeHs(self, arg: bool, /) -> None: ...

    def __setattr__(self, name: str, value: object | None) -> None: ...

@overload
def MolFromSmiles(SMILES: str | bytes, params: SmilesParserParams) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a molecule from a SMILES string.
           ARGUMENTS:

             - SMILES: the smiles string

             - params: used to provide optional parameters for the SMILES
             parsing

           RETURNS:

             a Mol object, None on failure.
    """

@overload
def MolFromSmiles(SMILES: str | bytes, sanitize: bool = True, replacements: Mapping[str, str] = {}) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a molecule from a SMILES string.
        ARGUMENTS:

          - SMILES: the smiles string

          - sanitize: (optional) toggles sanitization of the molecule.
            Defaults to True.

          - replacements: (optional) a dictionary of replacement strings (see
          below)
            Defaults to {}.

        RETURNS:

          a Mol object, None on failure.

         The optional replacements dict can be used to do string substitution
         of abbreviations in the input SMILES. The set of substitutions is
         repeatedly looped through until the string no longer changes. It is
         the responsibility of the caller to make sure that substitutions
         results in legal and sensible SMILES.

         Examples of replacements:

           CC{Q}C with {'{Q}':'OCCO'} -> CCOCCOC
           C{A}C{Q}C with {'{Q}':'OCCO', '{A}':'C1(CC1)'} -> CC1(CC1)COCCOC
           C{A}C{Q}C with {'{Q}':'{X}CC{X}', '{A}':'C1CC1', '{X}':'N'} ->
           CC1CC1CNCCNC
    """

def AtomFromSmiles(SMILES: str) -> rdkit.Chem.rdchem.Atom | None:
    """Construct an atom from a SMILES string"""

def BondFromSmiles(SMILES: str) -> rdkit.Chem.rdchem.Bond | None:
    """Construct a bond from a SMILES string"""

@overload
def MolFromSmarts(SMARTS: str | bytes, mergeHs: bool = False, replacements: Mapping[str, str] = {}) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a molecule from a SMARTS string.
        ARGUMENTS:

          - SMARTS: the smarts string

          - mergeHs: (optional) toggles the merging of explicit Hs in the query
          into the attached
            atoms.  So, for example, 'C[H]' becomes '[C;!H0]'.
            Defaults to 0.

          - replacements: (optional) a dictionary of replacement strings (see
          below)
            Defaults to {}. See the documentation for MolFromSmiles for an
            explanation.

        RETURNS:

          a Mol object, None on failure.
    """

@overload
def MolFromSmarts(SMARTS: str | bytes, params: SmartsParserParams) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a molecule from a SMARTS string.
           ARGUMENTS:

             - SMARTS: the smarts string

             - params: used to provide optional parameters for the SMARTS
             parsing

           RETURNS:

             a Mol object, None on failure.
    """

def AtomFromSmarts(SMARTS: str) -> rdkit.Chem.rdchem.Atom | None:
    """Construct an atom from a SMARTS string"""

def BondFromSmarts(SMARTS: str) -> rdkit.Chem.rdchem.Bond | None:
    """Construct a bond from a SMARTS string"""

class SmilesWriteParams:
    """Parameters controlling SMILES writing"""

    def __init__(self) -> None: ...

    @property
    def doIsomericSmiles(self) -> bool:
        """include stereochemistry and isotope information"""

    @doIsomericSmiles.setter
    def doIsomericSmiles(self, arg: bool, /) -> None: ...

    @property
    def doKekule(self) -> bool:
        """
        kekulize the molecule before generating the SMILES and output single/double bonds. NOTE that the output is not canonical and that this will thrown an exception if the molecule cannot be kekulized
        """

    @doKekule.setter
    def doKekule(self, arg: bool, /) -> None: ...

    @property
    def canonical(self) -> bool:
        """generate canonical SMILES"""

    @canonical.setter
    def canonical(self, arg: bool, /) -> None: ...

    @property
    def cleanStereo(self) -> bool:
        """chiral centers are removed if they have duplicate sidechains"""

    @cleanStereo.setter
    def cleanStereo(self, arg: bool, /) -> None: ...

    @property
    def allBondsExplicit(self) -> bool:
        """include symbols for all bonds"""

    @allBondsExplicit.setter
    def allBondsExplicit(self, arg: bool, /) -> None: ...

    @property
    def allHsExplicit(self) -> bool:
        """provide hydrogen counts for every atom"""

    @allHsExplicit.setter
    def allHsExplicit(self, arg: bool, /) -> None: ...

    @property
    def doRandom(self) -> bool:
        """randomize the output order. The resulting SMILES is not canonical"""

    @doRandom.setter
    def doRandom(self, arg: bool, /) -> None: ...

    @property
    def rootedAtAtom(self) -> int:
        """
        make sure the SMILES starts at the specified atom. The resulting SMILES is not canonical
        """

    @rootedAtAtom.setter
    def rootedAtAtom(self, arg: int, /) -> None: ...

    @property
    def includeDativeBonds(self) -> bool:
        """
        include the RDKit extension for dative bonds. Otherwise dative bonds will be written as single bonds
        """

    @includeDativeBonds.setter
    def includeDativeBonds(self, arg: bool, /) -> None: ...

    @property
    def ignoreAtomMapNumbers(self) -> bool:
        """ignore atom map numbers when canonicalizing the molecule"""

    @ignoreAtomMapNumbers.setter
    def ignoreAtomMapNumbers(self, arg: bool, /) -> None: ...

    def __setattr__(self, name: str, value: object | None) -> None: ...

@overload
def MolToSmiles(mol: rdkit.Chem.rdchem.Mol, params: SmilesWriteParams) -> str:
    """Returns the canonical SMILES string for a molecule"""

@overload
def MolToSmiles(mol: rdkit.Chem.rdchem.Mol, isomericSmiles: bool = True, kekuleSmiles: bool = False, rootedAtAtom: int = -1, canonical: bool = True, allBondsExplicit: bool = False, allHsExplicit: bool = False, doRandom: bool = False, ignoreAtomMapNumbers: bool = False) -> str:
    """
    Returns the canonical SMILES string for a molecule
        ARGUMENTS:

          - mol: the molecule
          - isomericSmiles: (optional) include information about
          stereochemistry in
            the SMILES.  Defaults to true.
          - kekuleSmiles: (optional) use the Kekule form (no aromatic bonds)
          in
            the SMILES.  Defaults to false.
          - rootedAtAtom: (optional) if non-negative, this forces the SMILES
            to start at a particular atom. Defaults to -1.  If not -1,
            overrides
            canonical setting.
          - canonical: (optional) if false no attempt will be made to
          canonicalize
            the molecule. Defaults to true.
          - allBondsExplicit: (optional) if true, all bond orders will be
          explicitly indicated
            in the output SMILES. Defaults to false.
          - allHsExplicit: (optional) if true, all H counts will be explicitly
          indicated
            in the output SMILES. Defaults to false.
          - doRandom: (optional) if true, randomize the traversal of the
          molecule graph,
            so we can generate random smiles. Defaults to false.  If true,
            overrides
            canonical setting.
          - ignoreAtomMapNumbers (optional) if true, ignores any atom map
          numbers when
            canonicalizing the molecule

        RETURNS:

          a string
    """

@overload
def MolFragmentToSmiles(mol: rdkit.Chem.rdchem.Mol, params: SmilesWriteParams, atomsToUse: Iterable[int], bondsToUse: Iterable[int] | None = None, atomSymbols: Iterable[str] | None = None, bondSymbols: Iterable[str] | None = None) -> str:
    """
    Returns the canonical SMILES string for a fragment of a
            molecule
        ARGUMENTS:

          - mol: the molecule
          - params: the SmilesWriteParams
          - atomsToUse : a list of atoms to include in the fragment
          - bondsToUse : (optional) a list of bonds to include in the fragment
            if not provided, all bonds between the atoms provided
            will be included.
          - atomSymbols : (optional) a list with the symbols to use for the
          atoms
            in the SMILES. This should have be mol.GetNumAtoms() long.
          - bondSymbols : (optional) a list with the symbols to use for the
          bonds
            in the SMILES. This should have be mol.GetNumBonds() long.

        RETURNS:

          a string
    """

@overload
def MolFragmentToSmiles(mol: rdkit.Chem.rdchem.Mol, atomsToUse: Iterable[int], bondsToUse: Iterable[int] | None = None, atomSymbols: Iterable[str] | None = None, bondSymbols: Iterable[str] | None = None, isomericSmiles: bool = True, kekuleSmiles: bool = False, rootedAtAtom: int = -1, canonical: bool = True, allBondsExplicit: bool = False, allHsExplicit: bool = False) -> str:
    """
    Returns the canonical SMILES string for a fragment of a
            molecule
        ARGUMENTS:

          - mol: the molecule
          - atomsToUse : a list of atoms to include in the fragment
          - bondsToUse : (optional) a list of bonds to include in the fragment
            if not provided, all bonds between the atoms provided
            will be included.
          - atomSymbols : (optional) a list with the symbols to use for the
          atoms
            in the SMILES. This should have be mol.GetNumAtoms() long.
          - bondSymbols : (optional) a list with the symbols to use for the
          bonds
            in the SMILES. This should have be mol.GetNumBonds() long.
          - isomericSmiles: (optional) include information about
          stereochemistry in
            the SMILES.  Defaults to true.
          - kekuleSmiles: (optional) use the Kekule form (no aromatic bonds)
          in
            the SMILES.  Defaults to false.
          - rootedAtAtom: (optional) if non-negative, this forces the SMILES
            to start at a particular atom. Defaults to -1.  If not -1,
            over-rides
            setting for canonical.
          - canonical: (optional) if false no attempt will be made to
          canonicalize
            the molecule. Defaults to true.
          - allBondsExplicit: (optional) if true, all bond orders will be
          explicitly indicated
            in the output SMILES. Defaults to false.
          - allHsExplicit: (optional) if true, all H counts will be explicitly
          indicated
            in the output SMILES. Defaults to false.

        RETURNS:

          a string
    """

class CXSmilesFields(enum.IntEnum):
    CX_NONE = 0

    CX_ATOM_LABELS = 1

    CX_MOLFILE_VALUES = 2

    CX_COORDS = 4

    CX_RADICALS = 8

    CX_ATOM_PROPS = 16

    CX_LINKNODES = 32

    CX_ENHANCEDSTEREO = 64

    CX_SGROUPS = 128

    CX_POLYMER = 256

    CX_BOND_CFG = 512

    CX_BOND_ATROPISOMER = 1024

    CX_COORDINATE_BONDS = 2048

    CX_ZERO_BONDS = 8192

    CX_ALL = 2147483647

    CX_ALL_BUT_COORDS = 2147483643

class SCSRBaseHbondOptions(enum.Enum):
    Ignore = 0

    UseSapAll = 1

    UseSapOne = 2

    Auto = 3

class SCSRTemplateNames(enum.Enum):
    UseFirstName = 1

    UseSecondName = 2

    AsEntered = 0

class RestoreBondDirOption(enum.Enum):
    RestoreBondDirOptionClear = 1

    RestoreBondDirOptionTrue = 0

@overload
def MolToCXSmiles(mol: rdkit.Chem.rdchem.Mol, params: SmilesWriteParams, flags: int = CXSmilesFields.CX_ALL, restoreBondDirs: RestoreBondDirOption = RestoreBondDirOption.RestoreBondDirOptionClear) -> str:
    """Returns the CXSMILES string for a molecule"""

@overload
def MolToCXSmiles(mol: rdkit.Chem.rdchem.Mol, isomericSmiles: bool = True, kekuleSmiles: bool = False, rootedAtAtom: int = -1, canonical: bool = True, allBondsExplicit: bool = False, allHsExplicit: bool = False, doRandom: bool = False) -> str:
    """
    Returns the CXSMILES string for a molecule
        ARGUMENTS:

          - mol: the molecule
          - isomericSmiles: (optional) include information about
          stereochemistry in
            the SMILES.  Defaults to true.
          - kekuleSmiles: (optional) use the Kekule form (no aromatic bonds)
          in
            the SMILES.  Defaults to false.
          - rootedAtAtom: (optional) if non-negative, this forces the SMILES
            to start at a particular atom. Defaults to -1.
          - canonical: (optional) if false no attempt will be made to
          canonicalize
            the molecule. Defaults to true.
          - allBondsExplicit: (optional) if true, all bond orders will be
          explicitly indicated
            in the output SMILES. Defaults to false.
          - allHsExplicit: (optional) if true, all H counts will be explicitly
          indicated
            in the output SMILES. Defaults to false.
          - doRandom: (optional) if true, randomizes the traversal of the
          molecule graph,
            so we can generate random smiles. Defaults to false.  If true,
            overrides
            canonical setting.

        RETURNS:

          a string
    """

@overload
def MolFragmentToCXSmiles(mol: rdkit.Chem.rdchem.Mol, params: SmilesWriteParams, atomsToUse: Iterable[int], bondsToUse: Iterable[int] | None = None, atomSymbols: Iterable[str] | None = None, bondSymbols: Iterable[str] | None = None) -> str:
    """
    Returns the CXSMILES string for a fragment of a molecule
        ARGUMENTS:

          - mol: the molecule
          - params: the SmilesWriteParams
          - atomsToUse : a list of atoms to include in the fragment
          - bondsToUse : (optional) a list of bonds to include in the fragment
            if not provided, all bonds between the atoms provided
            will be included.
          - atomSymbols : (optional) a list with the symbols to use for the
          atoms
            in the SMILES. This should have be mol.GetNumAtoms() long.
          - bondSymbols : (optional) a list with the symbols to use for the
          bonds
            in the SMILES. This should have be mol.GetNumBonds() long.

        RETURNS:

          a string
    """

@overload
def MolFragmentToCXSmiles(mol: rdkit.Chem.rdchem.Mol, atomsToUse: Iterable[int], bondsToUse: Iterable[int] | None = None, atomSymbols: Iterable[str] | None = None, bondSymbols: Iterable[str] | None = None, isomericSmiles: bool = True, kekuleSmiles: bool = False, rootedAtAtom: int = -1, canonical: bool = True, allBondsExplicit: bool = False, allHsExplicit: bool = False) -> str:
    """
    Returns the CXSMILES string for a fragment of a molecule
        ARGUMENTS:

          - mol: the molecule
          - atomsToUse : a list of atoms to include in the fragment
          - bondsToUse : (optional) a list of bonds to include in the fragment
            if not provided, all bonds between the atoms provided
            will be included.
          - atomSymbols : (optional) a list with the symbols to use for the
          atoms
            in the SMILES. This should have be mol.GetNumAtoms() long.
          - bondSymbols : (optional) a list with the symbols to use for the
          bonds
            in the SMILES. This should have be mol.GetNumBonds() long.
          - isomericSmiles: (optional) include information about
          stereochemistry in
            the SMILES.  Defaults to true.
          - kekuleSmiles: (optional) use the Kekule form (no aromatic bonds)
          in
            the SMILES.  Defaults to false.
          - rootedAtAtom: (optional) if non-negative, this forces the SMILES
            to start at a particular atom. Defaults to -1.  If not -1,
            overrides
            canonical setting.
          - canonical: (optional) if false no attempt will be made to
          canonicalize
            the molecule. Defaults to true.
          - allBondsExplicit: (optional) if true, all bond orders will be
          explicitly indicated
            in the output SMILES. Defaults to false.
          - allHsExplicit: (optional) if true, all H counts will be explicitly
          indicated
            in the output SMILES. Defaults to false.

        RETURNS:

          a string
    """

@overload
def MolToSmarts(mol: rdkit.Chem.rdchem.Mol, isomericSmiles: bool = True, rootedAtAtom: int = -1) -> str:
    """
    Returns a SMARTS string for a molecule
        ARGUMENTS:

          - mol: the molecule
          - isomericSmiles: (optional) include information about
          stereochemistry in
            the SMARTS.  Defaults to true.
          - rootedAtomAtom: (optional) the atom index to start the SMARTS
          from.

        RETURNS:

          a string
    """

@overload
def MolToSmarts(mol: rdkit.Chem.rdchem.Mol, params: SmilesWriteParams) -> str:
    """
    Returns a SMARTS string for a molecule
        ARGUMENTS:

          - mol: the molecule
          - params: SmilesWriteParams controlling the SMARTS generation

        RETURNS:

          a string
    """

def MolFragmentToSmarts(mol: rdkit.Chem.rdchem.Mol, atomsToUse: Iterable[int], bondsToUse: Iterable[int] | None = None, isomericSmarts: bool = True) -> str:
    """
    Returns a SMARTS string for a fragment of a molecule
        ARGUMENTS:

          - mol: the molecule
          - atomsToUse: indices of atoms to include in the SMARTS string
          - bondsToUse: indices of bonds to include in the SMARTS string
          (optional)
          - isomericSmarts: (optional) include information about
          stereochemistry in
            the SMARTS.  Defaults to true.

        RETURNS:

          a string
    """

def MolToCXSmarts(mol: rdkit.Chem.rdchem.Mol, isomericSmiles: bool = True) -> str:
    """
    Returns a SMARTS string for a molecule
        ARGUMENTS:

          - mol: the molecule
          - isomericSmiles: (optional) include information about
          stereochemistry in
            the SMARTS.  Defaults to true.

        RETURNS:

          a string
    """

def MolFragmentToCXSmarts(mol: rdkit.Chem.rdchem.Mol, atomsToUse: Iterable[int], bondsToUse: Iterable[int] | None = None, isomericSmarts: bool = True) -> str:
    """
    Returns a SMARTS string for a fragment of a molecule
        ARGUMENTS:

          - mol: the molecule
          - atomsToUse: indices of atoms to include in the SMARTS string
          - bondsToUse: indices of bonds to include in the SMARTS string
          (optional)
          - isomericSmarts: (optional) include information about
          stereochemistry in
            the SMARTS.  Defaults to true.

        RETURNS:

          a string
    """

def MolToTPLFile(mol: rdkit.Chem.rdchem.Mol, fileName: str, partialChargeProp: str = '_GasteigerCharge', writeFirstConfTwice: bool = False) -> None:
    """
    Writes a molecule to a TPL file.
        ARGUMENTS:

          - mol: the molecule
          - fileName: name of the file to write
          - partialChargeProp: name of the property to use for partial charges
            Defaults to '_GasteigerCharge'.
          - writeFirstConfTwice: Defaults to False.
            This should be set to True when writing TPLs to be read by
            the CombiCode.
    """

def MolToTPLBlock(mol: rdkit.Chem.rdchem.Mol, partialChargeProp: str = '_GasteigerCharge', writeFirstConfTwice: bool = False) -> str:
    """
    Returns the Tpl block for a molecule.
        ARGUMENTS:

          - mol: the molecule
          - partialChargeProp: name of the property to use for partial charges
            Defaults to '_GasteigerCharge'.
          - writeFirstConfTwice: Defaults to False.
            This should be set to True when writing TPLs to be read by
            the CombiCode.

        RETURNS:

          a string
    """

def MolFromPDBFile(pdbFileName: str | os.PathLike, sanitize: bool = True, removeHs: bool = True, flavor: int = 0, proximityBonding: bool = True) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a molecule from a PDB file.
        ARGUMENTS:

          - pdbFileName: name of the file to read

          - sanitize: (optional) toggles sanitization of the molecule.
            Defaults to true.

          - removeHs: (optional) toggles removing hydrogens from the molecule.
            This only make sense when sanitization is done.
            Defaults to true.

          - flavor: (optional)

          - proximityBonding: (optional) toggles automatic proximity bonding

        RETURNS:

          a Mol object, None on failure.
    """

def MolFromPDBBlock(molBlock: str | bytes, sanitize: bool = True, removeHs: bool = True, flavor: int = 0, proximityBonding: bool = True) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a molecule from a PDB block.
        ARGUMENTS:

          - molBlock: string containing the PDB block

          - sanitize: (optional) toggles sanitization of the molecule.
            Defaults to True.

          - removeHs: (optional) toggles removing hydrogens from the molecule.
            This only make sense when sanitization is done.
            Defaults to true.

          - flavor: (optional)

          - proximityBonding: (optional) toggles automatic proximity bonding

        RETURNS:

          a Mol object, None on failure.
    """

def MolToPDBBlock(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, flavor: int = 0) -> str:
    """
    Returns a PDB block for a molecule
        ARGUMENTS:

          - mol: the molecule
          - confId: (optional) selects which conformation to output (-1 =
          default)
          - flavor: (optional)
                  - flavor & 1 : Write MODEL/ENDMDL lines around each record
                  - flavor & 2 : Don't write any CONECT records
                  - flavor & 4 : Write CONECT records in both directions
                  - flavor & 8 : Don't use multiple CONECTs to encode bond
                  order

                  - flavor & 16 : Write MASTER record
                  - flavor & 32 : Write TER record

        RETURNS:

          a string
    """

def MolToPDBFile(mol: rdkit.Chem.rdchem.Mol, filename: str, confId: int = -1, flavor: int = 0) -> None:
    """
    Writes a PDB file for a molecule
        ARGUMENTS:

          - mol: the molecule
          - filename: name of the file to write
          - confId: (optional) selects which conformation to output (-1 =
          default)
          - flavor: (optional)
                  - flavor & 1 : Write MODEL/ENDMDL lines around each record
                  - flavor & 2 : Don't write any CONECT records
                  - flavor & 4 : Write CONECT records in both directions
                  - flavor & 8 : Don't use multiple CONECTs to encode bond
                  order

                  - flavor & 16 : Write MASTER record
                  - flavor & 32 : Write TER record
    """

def MolFromSequence(text: str | bytes, sanitize: bool = True, flavor: int = 0) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a molecule from a sequence string (currently
            supports standard amino acids, DNA and RNA bases).
        ARGUMENTS:

          - text: string containing the sequence

          - sanitize: (optional) toggles sanitization of the molecule.
            Defaults to True.

          - flavor: (optional)
              - 0 Protein, L amino acids (default)
              - 1 Protein, D amino acids
              - 2 RNA, no cap
              - 3 RNA, 5' cap
              - 4 RNA, 3' cap
              - 5 RNA, both caps
              - 6 DNA, no cap
              - 7 DNA, 5' cap
              - 8 DNA, 3' cap
              - 9 DNA, both caps

        RETURNS:

          a Mol object, None on failure.
    """

def MolToSequence(mol: rdkit.Chem.rdchem.Mol) -> str:
    """
    Returns the sequence string for a molecule
        ARGUMENTS:

          - mol: the molecule

        NOTE: the molecule should contain monomer information in
        AtomMonomerInfo structures

        RETURNS:

          a string
    """

def MolFromFASTA(text: str | bytes, sanitize: bool = True, flavor: int = 0) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a molecule from a FASTA string (currently supports
            standard amino acids, DNA and RNA bases).
        ARGUMENTS:

          - text: string containing the FASTA

          - sanitize: (optional) toggles sanitization of the molecule.
            Defaults to True.

      - flavor: (optional)
          - 0 Protein, L amino acids (default)
          - 1 Protein, D amino acids
          - 2 RNA, no cap
          - 3 RNA, 5' cap
          - 4 RNA, 3' cap
          - 5 RNA, both caps
          - 6 DNA, no cap
          - 7 DNA, 5' cap
          - 8 DNA, 3' cap
          - 9 DNA, both caps
        RETURNS:

          a Mol object, None on failure.
    """

def MolToFASTA(mol: rdkit.Chem.rdchem.Mol) -> str:
    """
    Returns the FASTA string for a molecule
        ARGUMENTS:

          - mol: the molecule

        NOTE: the molecule should contain monomer information in
        AtomMonomerInfo structures

        RETURNS:

          a string
    """

def MolFromHELM(text: str | bytes, sanitize: bool = True) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a molecule from a HELM string (currently supports
            standard
            amino acids, DNA and RNA bases).
        ARGUMENTS:

          - text: string containing the HELM

          - sanitize: (optional) toggles sanitization of the molecule.
            Defaults to true.

        RETURNS:

          a Mol object, None on failure.
    """

def MolToHELM(mol: rdkit.Chem.rdchem.Mol) -> str:
    """
    Returns the HELM string for a molecule
        ARGUMENTS:

          - mol: the molecule

        NOTE: the molecule should contain monomer information in
        AtomMonomerInfo structures

        RETURNS:

          a string
    """

def CanonicalRankAtoms(mol: rdkit.Chem.rdchem.Mol, breakTies: bool = True, includeChirality: bool = True, includeIsotopes: bool = True, includeAtomMaps: bool = True, includeChiralPresence: bool = False) -> list[int]:
    """
    Returns the canonical atom ranking for each atom of a
            molecule fragment.
        If breakTies is False, this returns the symmetry class for each atom.
        The symmetry class is used by the canonicalization routines to type
        each atom based on the whole chemistry of the molecular graph.  Any
        atom with the same rank (symmetry class) is indistinguishable.  For
        example:

          >>> mol = MolFromSmiles('C1NCN1')
          >>> list(CanonicalRankAtoms(mol, breakTies=False))
          [0,1,0,1]

        In this case the carbons have the same symmetry class and the nitrogens
        have the same
        symmetry class.  From the perspective of the Molecular Graph, they are
        identical.

        ARGUMENTS:

          - mol: the molecule
          - breakTies: (optional) force breaking of ranked ties [default=True]
          - includeChirality: (optional) use chiral information when computing
          rank [default=True]
          - includeIsotopes: (optional) use isotope information when computing
          rank [default=True]
          - includeAtomMaps: (optional) use atom map information when computing
          rank [default=True]
          - includeChiralPresence: (optional) use information about whether or
          not chirality is specified when computing rank [default=False]

        RETURNS:

          a string
    """

def CanonicalRankAtomsInFragment(mol: rdkit.Chem.rdchem.Mol, atomsToUse: Iterable[int], bondsToUse: Iterable[int] | None = None, atomSymbols: Iterable[str] | None = None, breakTies: bool = True, includeChirality: bool = True, includeIsotopes: bool = True, includeAtomMaps: bool = True, includeChiralPresence: bool = False) -> list[int]:
    """
    Returns the canonical atom ranking for each atom of a
            molecule fragment
        See help(CanonicalRankAtoms) for more information.

         >>> mol = MolFromSmiles('C1NCN1.C1NCN1')
         >>> list(CanonicalRankAtomsInFragment(mol, atomsToUse=range(0,4),
         breakTies=False))
         [4,6,4,6,-1,-1,-1,-1]
         >>> list(CanonicalRankAtomsInFragment(mol, atomsToUse=range(4,8),
         breakTies=False))
         [-1,-1,-1,-1,4,6,4,6]

        ARGUMENTS:

          - mol: the molecule
          - atomsToUse : a list of atoms to include in the fragment
          - bondsToUse : (optional) a list of bonds to include in the fragment
            if not provided, no bonds will be used
          - atomSymbols : (optional) a list with the symbols to use for the
          atoms
            in the SMILES. This should have be mol.GetNumAtoms() long.
          - breakTies: (optional) force breaking of ranked ties
          - includeChirality: (optional) use chiral information when computing
          rank [default=True]
          - includeIsotopes: (optional) use isotope information when computing
          rank [default=True]
          - includeAtomMaps: (optional) use atom map information when computing
          rank [default=True]
          - includeChiralPresence: (optional) use information about whether or
          not chirality is specified when computing rank [default=False]

        RETURNS:

          a string
    """

def CanonicalizeEnhancedStereo(mol: rdkit.Chem.rdchem.Mol) -> None: ...

def CreateAtomIntPropertyList(mol: rdkit.Chem.rdchem.Mol, propName: str, missingValueMarker: str = '', lineSize: int = 190) -> None:
    """
    creates a list property on the molecule from individual atom property values
    """

def CreateAtomDoublePropertyList(mol: rdkit.Chem.rdchem.Mol, propName: str, missingValueMarker: str = '', lineSize: int = 190) -> None:
    """
    creates a list property on the molecule from individual atom property values
    """

def CreateAtomBoolPropertyList(mol: rdkit.Chem.rdchem.Mol, propName: str, missingValueMarker: str = '', lineSize: int = 190) -> None:
    """
    creates a list property on the molecule from individual atom property values
    """

def CreateAtomStringPropertyList(mol: rdkit.Chem.rdchem.Mol, propName: str, missingValueMarker: str = '', lineSize: int = 190) -> None:
    """
    creates a list property on the molecule from individual atom property values
    """

def CreateBondIntPropertyList(mol: rdkit.Chem.rdchem.Mol, propName: str, missingValueMarker: str = '', lineSize: int = 190) -> None:
    """
    creates a list property on the molecule from individual bond property values
    """

def CreateBondDoublePropertyList(mol: rdkit.Chem.rdchem.Mol, propName: str, missingValueMarker: str = '', lineSize: int = 190) -> None:
    """
    creates a list property on the molecule from individual bond property values
    """

def CreateBondBoolPropertyList(mol: rdkit.Chem.rdchem.Mol, propName: str, missingValueMarker: str = '', lineSize: int = 190) -> None:
    """
    creates a list property on the molecule from individual bond property values
    """

def CreateBondStringPropertyList(mol: rdkit.Chem.rdchem.Mol, propName: str, missingValueMarker: str = '', lineSize: int = 190) -> None:
    """
    creates a list property on the molecule from individual bond property values
    """

def MolToRandomSmilesVect(mol: rdkit.Chem.rdchem.Mol, numSmiles: int, randomSeed: int = 0, isomericSmiles: bool = True, kekuleSmiles: bool = False, allBondsExplicit: bool = False, allHsExplicit: bool = False) -> list[str]:
    """returns a list of SMILES generated using the randomSmiles algorithm"""

class PNGMetadataParams:
    """Parameters controlling metadata included in PNG images"""

    def __init__(self) -> None: ...

    @property
    def includePkl(self) -> bool:
        """toggles inclusion of molecule pickle (default=True)"""

    @includePkl.setter
    def includePkl(self, arg: bool, /) -> None: ...

    @property
    def includeSmiles(self) -> bool:
        """toggles inclusion of molecule CXSMILES (default=True)"""

    @includeSmiles.setter
    def includeSmiles(self, arg: bool, /) -> None: ...

    @property
    def includeMol(self) -> bool:
        """toggles inclusion of molecule molblock (default=False)"""

    @includeMol.setter
    def includeMol(self, arg: bool, /) -> None: ...

    @property
    def propertyFlags(self) -> int:
        """
        choose properties to be included in the pickle (default=rdkit.Chem.rdchem.PropertyPickleOptions.NoProps)
        """

    @propertyFlags.setter
    def propertyFlags(self, arg: int, /) -> None: ...

    @property
    def smilesWriteParams(self) -> SmilesWriteParams:
        """
        choose SmilesWriteParams for the CXSMILES string (default=rdkit.Chem.rdmolfiles.SmilesWriteParams())
        """

    @smilesWriteParams.setter
    def smilesWriteParams(self, arg: SmilesWriteParams, /) -> None: ...

    @property
    def cxSmilesFlags(self) -> int:
        """
        choose CXSMILES fields to be included in the CXSMILES string (default=rdkit.Chem.rdmolfiles.CXSmilesFields.CX_ALL)
        """

    @cxSmilesFlags.setter
    def cxSmilesFlags(self, arg: int, /) -> None: ...

    @property
    def restoreBondDirs(self) -> RestoreBondDirOption:
        """
        choose what to do with bond dirs in the CXSMILES string (default=rdkit.Chem.rdmolfiles.RestoreBondDirOption.RestoreBondDirOptionClear)
        """

    @restoreBondDirs.setter
    def restoreBondDirs(self, arg: RestoreBondDirOption, /) -> None: ...

    def __setattr__(self, name: str, value: object | None) -> None: ...

def MolFromPNGString(png: bytes, params: SmilesParserParams | None = None) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a molecule from metadata in a PNG string.

           ARGUMENTS:

             - png: the PNG string

             - params: used to provide optional parameters for the metadata
             parsing

           RETURNS:
             a Mol object, None on failure.
    """

def MolFromPNGFile(filename: str | os.PathLike, params: SmilesParserParams | None = None) -> rdkit.Chem.rdchem.Mol | None:
    """
    Construct a molecule from metadata in a PNG file.

           ARGUMENTS:

             - filename: the PNG filename

             - params: used to provide optional parameters for the metadata
             parsing

           RETURNS:
             a Mol object, None on failure.
    """

def MolsFromPNGString(png: bytes, tag: str = 'rdkitPKL', params: SmilesParserParams | None = None) -> tuple[rdkit.Chem.rdchem.Mol, ...]:
    """returns a tuple of molecules constructed from the PNG string"""

def MolsFromPNGFile(filename: str | os.PathLike, tag: str = 'rdkitPKL', params: SmilesParserParams | None = None) -> tuple[rdkit.Chem.rdchem.Mol, ...]:
    """returns a tuple of molecules constructed from the PNG file"""

@overload
def MolsFromCDXMLFile(filename: str | os.PathLike, sanitize: bool = True, removeHs: bool = True) -> tuple[rdkit.Chem.rdchem.Mol, ...]:
    """
    Construct a molecule from a cdxml file.

           Note that the CDXML format is large and complex, the RDKit doesn't
           support full functionality, just the base ones required for molecule
           and reaction parsing.

           ARGUMENTS:

             - filename: the cdxml filename

             - sanitize: if True, sanitize the molecules [default True]

             - removeHs: if True, convert explicit Hs into implicit Hs.
             [default True]

           RETURNS:
             an iterator of parsed Mol objects.
    """

@overload
def MolsFromCDXMLFile(filename: str | os.PathLike, params: CDXMLParserParams) -> tuple[rdkit.Chem.rdchem.Mol, ...]:
    """
    Construct a molecule from a cdxml file.

           Note: that the CDXML format is large and complex, the RDKit doesn't
           support full functionality, just the base ones required for molecule
           and reaction parsing.

           Note: If the ChemDraw extensions are available,
              CDXMLFormat::Auto attempts to see if the input string is CDXML or
              CDX,
           If not, it defaults to CDXML

           ARGUMENTS:

             - filename: the cdxml filename

             - pyParams: CDXParserParams, see CDXParserParams for usage

           RETURNS:
             a tuple  of parsed Mol objects.
    """

@overload
def MolsFromCDXML(cdxml: str | bytes, sanitize: bool = True, removeHs: bool = True) -> tuple[rdkit.Chem.rdchem.Mol, ...]:
    """
    Construct a molecule from a cdxml string.

           Note that the CDXML format is large and complex, the RDKit doesn't
           support full functionality, just the base ones required for molecule
           and reaction parsing.

           ARGUMENTS:

             - cdxml: the cdxml string

             - sanitize: if True, sanitize the molecules [default True]

             - removeHs: if True, convert explicit Hs into implicit Hs.
             [default True]

           RETURNS:
             an iterator of parsed Mol objects.
    """

@overload
def MolsFromCDXML(cdxml: str | bytes, params: CDXMLParserParams) -> tuple[rdkit.Chem.rdchem.Mol, ...]:
    """
    Construct a molecule from a cdxml string.

           Note that the CDXML format is large and complex, the RDKit doesn't
           support full functionality, just the base ones required for molecule
           and reaction parsing.

           Note: in this function CDXMLFormat::Auto currently defaults to CDXML

           ARGUMENTS:

             - cdxml: the cdxml string

             - pyParams: CDXParserParams, see CDXParserParams for usage

           RETURNS:
             a tuple of parsed Mol objects.
    """

class CDXMLFormat(enum.Enum):
    CDXML = 0

    CDX = 1

    Auto = 2

class CDXMLParserParams:
    """Parameters controlling conversion of a CDXML document to molecules"""

    @overload
    def __init__(self) -> None:
        """Construct a default CDXMLFormat"""

    @overload
    def __init__(self, sanitize: bool, removeHs: bool, format: CDXMLFormat) -> None: ...

    @property
    def sanitize(self) -> bool:
        """
        controls whether or not the molecule is sanitized before being returned
        """

    @sanitize.setter
    def sanitize(self, arg: bool, /) -> None: ...

    @property
    def removeHs(self) -> bool:
        """controls whether or not Hs are removed before the molecule is returned"""

    @removeHs.setter
    def removeHs(self, arg: bool, /) -> None: ...

    @property
    def format(self) -> CDXMLFormat:
        """
        ChemDraw format One of Auto, CDXML, CDX.  For data streams, Auto defaults to CDXML
        """

    @format.setter
    def format(self, arg: CDXMLFormat, /) -> None: ...

    def __setattr__(self, name: str, value: object | None) -> None: ...

def MolToCDXMLBlock(mol: rdkit.Chem.rdchem.Mol, format: CDXMLFormat = CDXMLFormat.CDXML) -> object:
    """
    brief write a CDX or CDXML block from a molecule

          The RDKit is optionally built with the Revvity ChemDraw parser
          If this is available, CDX and CDXML can be written
            Note that the CDXML format is large and complex, the RDKit doesn't
            support full functionality, just the base ones required for molecule and
            reaction parsing.

          Note: If the ChemDraw extensions are unavailable, an exception will be thrown
           please use the support function HasChemDrawCDXSupport() to check
           whether ChemDraw writing support is enabled.

          Note: For CDXML this returns a UTF-8 string <str>
                For CDX this returns a byte sting <bytes>

          ARGUMENTS:

            - mol: the molecule to write

            - format: CDXMLFormat [default CDXML]

          RETURNS:
            the CDXML or CDX block
    """

def HasChemDrawCDXSupport() -> bool:
    """
    Returns true if the RDKit is built with ChemDraw CDX
            support
    """

@overload
def MolMetadataToPNGFile(mol: rdkit.Chem.rdchem.Mol, filename: str | os.PathLike, includePkl: bool = True, includeSmiles: bool = True, includeMol: bool = False) -> bytes:
    """
    Adds molecular metadata to PNG data read from a file.

         ARGUMENTS:

           - mol: the molecule

           - filename: the PNG filename

           - includePkl: include the RDKit's internal binary format in the output

           - includeSmiles: include CXSmiles in the output

           - includeMol: include CTAB (Mol) in the output

         RETURNS:
           the updated PNG data
    """

@overload
def MolMetadataToPNGFile(mol: rdkit.Chem.rdchem.Mol, filename: str | os.PathLike, params: PNGMetadataParams) -> bytes:
    """
    Adds molecular metadata to PNG data read from a file.

         ARGUMENTS:

           - mol: the molecule

           - filename: the PNG filename

           - params: an instance of PNGMetadataParams

         RETURNS:
           the updated PNG data
    """

@overload
def MolMetadataToPNGString(mol: rdkit.Chem.rdchem.Mol, png: bytes, includePkl: bool = True, includeSmiles: bool = True, includeMol: bool = False) -> bytes:
    """
    Adds molecular metadata to a PNG string.

         ARGUMENTS:

           - mol: the molecule

           - png: the PNG string

           - includePkl: include the RDKit's internal binary format in the output

           - includeSmiles: include CXSmiles in the output

           - includeMol: include CTAB (Mol) in the output

         RETURNS:
           the updated PNG data
    """

@overload
def MolMetadataToPNGString(mol: rdkit.Chem.rdchem.Mol, png: bytes, params: PNGMetadataParams) -> bytes:
    """
    Adds molecular metadata to a PNG string.

         ARGUMENTS:

           - mol: the molecule

           - png: the PNG string

           - params: an instance of PNGMetadataParams

         RETURNS:
           the updated PNG data
    """

def AddMetadataToPNGFile(metadata: dict[str, str], filename: str | os.PathLike) -> bytes:
    """
    Adds metadata to PNG data read from a file.

         ARGUMENTS:

           - metadata: dict with the metadata to be written
                       (keys and values should be strings)

           - filename: the PNG filename

         RETURNS:
           the updated PNG data
    """

def AddMetadataToPNGString(metadata: dict[str, str], png: bytes) -> bytes:
    """
    Adds metadata to a PNG string.

         ARGUMENTS:

           - metadata: dict with the metadata to be written
                       (keys and values should be strings)

           - png: the PNG string

         RETURNS:
           the updated PNG data
    """

def MetadataFromPNGFile(filename: object, asList: bool = False) -> object:
    """
    Returns a dict with all metadata from the PNG file. Keys are strings, values are bytes. If asList is True, a list of (key, value) tuples is returned; this enables retrieving multiple values sharing the same key.
    """

def MetadataFromPNGString(png: bytes, asList: bool = False) -> object:
    """
    Returns a dict with all metadata from the PNG string. Keys are strings, values are bytes. If asList is True, a list of (key, value) tuples is returned; this enables retrieving multiple values sharing the same key.
    """

class SDMolSupplier:
    """
    A class which supplies molecules from an SD file.

         Usage examples:

         1) Lazy evaluation: the molecules are not constructed until we ask for them:

              >>> suppl = SDMolSupplier('in.sdf')
              >>> for mol in suppl:
              ...    mol.GetNumAtoms()

         2) Lazy evaluation 2:

              >>> suppl = SDMolSupplier('in.sdf')
              >>> mol1 = next(suppl)
              >>> mol2 = next(suppl)
              >>> suppl.reset()
              >>> mol3 = next(suppl)
              # mol3 and mol1 are the same:
              >>> MolToSmiles(mol3)==MolToSmiles(mol1)

         3) Random Access:

              >>> suppl = SDMolSupplier('in.sdf')
              >>> mol1 = suppl[0]
              >>> mol2 = suppl[1]
              # NOTE: this will generate an IndexError if the supplier doesn't have that many
              molecules.

         4) Random Access 2: looping over all molecules

              >>> suppl = SDMolSupplier('in.sdf')
              >>> nMols = len(suppl)
              >>> for i in range(nMols):
              ...   suppl[i].GetNumAtoms()

      Properties in the SD file are used to set properties on each molecule.
      The properties are accessible using the mol.GetProp(propName) method.
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, fileName: str | os.PathLike, sanitize: bool = True, removeHs: bool = True, strictParsing: bool = True) -> None: ...

    def __enter__(self) -> SDMolSupplier: ...

    def __exit__(self, excType: object | None = None, excValue: object | None = None, traceback: object | None = None) -> bool: ...

    def __iter__(self) -> SDMolSupplier: ...

    def __next__(self) -> rdkit.Chem.rdchem.Mol | None:
        """Returns the next molecule in the file. Raises _StopIteration_ on EOF."""

    def __getitem__(self, idx: int) -> rdkit.Chem.rdchem.Mol | None: ...

    def reset(self) -> None:
        """Resets our position in the file to the beginning."""

    def __len__(self) -> int: ...

    @overload
    def SetData(self, data: str, sanitize: bool = True, removeHs: bool = True, strictParsing: bool = True) -> None: ...

    @overload
    def SetData(self, data: bytes, sanitize: bool = True, removeHs: bool = True, strictParsing: bool = True) -> None:
        """Sets the text to be parsed."""

    def GetItemText(self, index: int) -> str:
        """Returns the text for an item."""

    def atEnd(self) -> bool:
        """Returns whether or not we have hit EOF."""

    def GetProcessPropertyLists(self) -> bool:
        """
        Returns whether or not any property lists that are present will be processed when reading molecules.
        """

    def SetProcessPropertyLists(self, val: bool) -> None:
        """
        Sets whether or not any property lists that are present will be processed when reading molecules.
        """

class ForwardSDMolSupplier:
    """
    A class which supplies molecules from a file-like object containing SD data.

      Usage examples:

        1) Lazy evaluation: the molecules are not constructed until we ask for them:

           >>> suppl = ForwardSDMolSupplier(file('in.sdf'))
           >>> for mol in suppl:
           ...    if mol is not None: mol.GetNumAtoms()

        2) we can also read from compressed files:

           >>> import gzip
           >>> suppl = ForwardSDMolSupplier(gzip.open('in.sdf.gz'))
           >>> for mol in suppl:
           ...   if mol is not None: print mol.GetNumAtoms()

      Properties in the SD file are used to set properties on each molecule.
      The properties are accessible using the mol.GetProp(propName) method.
    """

    @overload
    def __init__(self, filename: str | os.PathLike, sanitize: bool = True, removeHs: bool = True, strictParsing: bool = True) -> None: ...

    @overload
    def __init__(self, fileobj: object, sanitize: bool = True, removeHs: bool = True, strictParsing: bool = True) -> None: ...

    def __enter__(self) -> ForwardSDMolSupplier: ...

    def __exit__(self, excType: object | None = None, excValue: object | None = None, traceback: object | None = None) -> bool: ...

    def __next__(self) -> rdkit.Chem.rdchem.Mol | None:
        """Returns the next molecule in the file. Raises _StopIteration_ on EOF."""

    def atEnd(self) -> bool:
        """Returns whether or not we have hit EOF."""

    def GetEOFHitOnRead(self) -> bool:
        """Returns whether EOF was hit while parsing the previous entry."""

    def __iter__(self) -> ForwardSDMolSupplier: ...

    def GetProcessPropertyLists(self) -> bool:
        """
        Returns whether or not any property lists that are present will be processed when reading molecules.
        """

    def SetProcessPropertyLists(self, val: bool) -> None:
        """
        Sets whether or not any property lists that are present will be processed when reading molecules.
        """

class TDTMolSupplier:
    """
    A class which supplies molecules from a TDT file.

      Usage examples:

        1) Lazy evaluation: the molecules are not constructed until we ask for them:

           >>> suppl = TDTMolSupplier('in.smi')
           >>> for mol in suppl:
           ...    mol.GetNumAtoms()

        2) Lazy evaluation 2:

           >>> suppl = TDTMolSupplier('in.smi')
           >>> mol1 = next(suppl)
           >>> mol2 = next(suppl)
           >>> suppl.reset()
           >>> mol3 = next(suppl)
           # mol3 and mol1 are the same:
           >>> MolToSmiles(mol3)==MolToSmiles(mol1)

        3) Random Access:  all molecules are constructed as soon as we ask for the
           length:

           >>> suppl = TDTMolSupplier('in.smi')
           >>> nMols = len(suppl)
           >>> for i in range(nMols):
           ...   suppl[i].GetNumAtoms()

      Properties in the file are used to set properties on each molecule.
      The properties are accessible using the mol.GetProp(propName) method.
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, fileName: str | os.PathLike, nameRecord: str = '', confId2D: int = -1, confId3D: int = -1, sanitize: bool = True) -> None: ...

    def __enter__(self) -> TDTMolSupplier: ...

    def __exit__(self, excType: object | None = None, excValue: object | None = None, traceback: object | None = None) -> bool: ...

    def __iter__(self) -> TDTMolSupplier: ...

    def __next__(self) -> rdkit.Chem.rdchem.Mol | None:
        """Returns the next molecule in the file. Raises _StopIteration_ on EOF."""

    def __getitem__(self, idx: int) -> rdkit.Chem.rdchem.Mol | None: ...

    def reset(self) -> None:
        """Resets our position in the file to the beginning."""

    def __len__(self) -> int: ...

    def SetData(self, data: str, nameRecord: str = '', confId2D: int = -1, confId3D: int = -1, sanitize: bool = True) -> None:
        """Sets the text to be parsed."""

    def GetItemText(self, index: int) -> str:
        """Returns the text for an item."""

class SmilesMolSupplier:
    """
    A class which supplies molecules from a text file.

      Usage examples:

        1) Lazy evaluation: the molecules are not constructed until we ask for them:

           >>> suppl = SmilesMolSupplier('in.smi')
           >>> for mol in suppl:
           ...    mol.GetNumAtoms()

        2) Lazy evaluation 2:

           >>> suppl = SmilesMolSupplier('in.smi')
           >>> mol1 = next(suppl)
           >>> mol2 = next(suppl)
           >>> suppl.reset()
           >>> mol3 = next(suppl)
           # mol3 and mol1 are the same:
           >>> MolToSmiles(mol3)==MolToSmiles(mol1)

        3) Random Access: all molecules are constructed as soon as we ask for the
           length:

           >>> suppl = SmilesMolSupplier('in.smi')
           >>> nMols = len(suppl)
           >>> for i in range(nMols):
           ...   suppl[i].GetNumAtoms()

      If the input file has a title line and more than two columns (smiles and id), the
      additional columns will be used to set properties on each molecule. The properties
      are accessible using the mol.GetProp(propName) method.
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, data: str | os.PathLike, delimiter: str = ' ', smilesColumn: int = 0, nameColumn: int = 1, titleLine: bool = True, sanitize: bool = True) -> None:
        """
        Constructor

          ARGUMENTS:

            - fileName: name of the file to be read

            - delimiter: (optional) text delimiter (a string). Defauts to ' '.

            - smilesColumn: (optional) index of the column containing the SMILES
              data. Defaults to 0.

            - nameColumn: (optional) index of the column containing molecule names.
              Defaults to 1.

            - titleLine: (optional) set this toggle if the file contains a title line.
              Defaults to 1.

            - sanitize: (optional) toggles sanitization of molecules as they are read.
              Defaults to 1.
        """

    def __enter__(self) -> SmilesMolSupplier: ...

    def __exit__(self, excType: object | None = None, excValue: object | None = None, traceback: object | None = None) -> bool: ...

    def __iter__(self) -> SmilesMolSupplier: ...

    def __next__(self) -> rdkit.Chem.rdchem.Mol | None:
        """Returns the next molecule in the file. Raises _StopIteration_ on EOF."""

    def __getitem__(self, idx: int) -> rdkit.Chem.rdchem.Mol | None: ...

    def reset(self) -> None:
        """Resets our position in the file to the beginning."""

    def __len__(self) -> int: ...

    def SetData(self, data: str, delimiter: str = ' ', smilesColumn: int = 0, nameColumn: int = 1, titleLine: bool = True, sanitize: bool = True) -> None:
        """Sets the text to be parsed."""

    def GetItemText(self, index: int) -> str:
        """Returns the text for an item."""

def SmilesMolSupplierFromText(text: str, delimiter: str = ' ', smilesColumn: int = 0, nameColumn: int = 1, titleLine: bool = True, sanitize: bool = True) -> SmilesMolSupplier: ...

class MaeMolSupplier:
    """
    A class which supplies molecules from a file-like object containing Maestro data.

      Usage examples:

        1) Lazy evaluation: the molecules are not constructed until we ask for them:

           >>> suppl = MaeMolSupplier(file('in.mae'))
           >>> for mol in suppl:
           ...    if mol is not None: mol.GetNumAtoms()

        2) we can also read from compressed files:

           >>> import gzip
           >>> suppl = MaeMolSupplier(gzip.open('in.maegz'))
           >>> for mol in suppl:
           ...   if mol is not None: print mol.GetNumAtoms()

      Properties in the Maestro file are used to set properties on each molecule.
      The properties are accessible using the mol.GetProp(propName) method.
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, filename: str | os.PathLike, sanitize: bool = True, removeHs: bool = True) -> None: ...

    @overload
    def __init__(self, fileobj: object, sanitize: bool = True, removeHs: bool = True) -> None: ...

    def __enter__(self) -> MaeMolSupplier: ...

    def __exit__(self, excType: object | None = None, excValue: object | None = None, traceback: object | None = None) -> bool: ...

    def __iter__(self) -> MaeMolSupplier: ...

    def __next__(self) -> rdkit.Chem.rdchem.Mol | None:
        """Returns the next molecule in the file. Raises _StopIteration_ on EOF."""

    def __getitem__(self, idx: int) -> rdkit.Chem.rdchem.Mol | None: ...

    def reset(self) -> None:
        """Resets our position in the file to the beginning."""

    def __len__(self) -> int: ...

    def SetData(self, data: str, sanitize: bool = True, removeHs: bool = True) -> None:
        """Sets the text to be parsed."""

    def atEnd(self) -> bool:
        """Returns whether or not we have hit EOF."""

class SmilesWriter:
    """A class for writing molecules to text files."""

    @overload
    def __init__(self, fileName: str | os.PathLike, delimiter: str = ' ', nameHeader: str = 'Name', includeHeader: bool = True, isomericSmiles: bool = True, kekuleSmiles: bool = False) -> None:
        """
        Constructor.

           ARGUMENTS:

             - fileName: name of the output file. ('-' to write to stdout)
             - delimiter: (optional) delimiter to be used to separate entries on each line.
             - nameHeader: (optional) text to use for the name column in the header line.
                           If this is blank, names will not be included in the output.
             - includeHeader: (optional) toggles inclusion of a header line in the output file.
             - isomericSmiles: (optional) toggles output of isomeric smiles
               (includes stereochem information).
             - kekuleSmiles: (optional) toggles output of kekule smiles (no aromatic
               bonds for molecules that have been kekulized).
        """

    @overload
    def __init__(self, fileObj: object, delimiter: str = ' ', nameHeader: str = 'Name', includeHeader: bool = True, isomericSmiles: bool = True, kekuleSmiles: bool = False) -> None: ...

    def __enter__(self) -> SmilesWriter: ...

    def __exit__(self, excType: object | None = None, excValue: object | None = None, traceback: object | None = None) -> bool: ...

    def SetProps(self, props: Sequence[str]) -> None:
        """
        Sets the properties to be written to the output file

          ARGUMENTS:

            - props: a list or tuple of property names
        """

    def write(self, mol: rdkit.Chem.rdchem.Mol, confId: int = -1) -> None:
        """
        Writes a molecule to the output file.

          ARGUMENTS:

            - mol: the Mol to be written
            - confId: (optional) ignored
        """

    def flush(self) -> None:
        """Flushes the output file (forces the disk file to be updated)."""

    def close(self) -> None:
        """
        Flushes the output file and closes it. The Writer cannot be used after this.
        """

    def NumMols(self) -> int:
        """Returns the number of molecules written so far."""

class SDWriter:
    """
    A class for writing molecules to SD files.

      Usage examples:

        1) writing to a named file:

           >>> writer = SDWriter('out.sdf')
           >>> for mol in list_of_mols:
           ...    writer.write(mol)

        2) writing to a file-like object:

           >>> import gzip
           >>> outf=gzip.open('out.sdf.gz','wt+')
           >>> writer = SDWriter(outf)
           >>> for mol in list_of_mols:
           ...   writer.write(mol)
           >>> writer.close()
           >>> outf.close()

      By default all non-private molecular properties are written to the SD file.
      This can be changed using the SetProps method:

           >>> writer = SDWriter('out.sdf')
           >>> writer.SetProps(['prop1','prop2'])
    """

    @overload
    def __init__(self, fileName: str | os.PathLike) -> None:
        """
        Constructor.

        If a string argument is provided, it will be treated as the name of the
        output file. If a file-like object is provided, output will be sent there.
        """

    @overload
    def __init__(self, fileObj: object) -> None: ...

    def __enter__(self) -> SDWriter: ...

    def __exit__(self, excType: object | None = None, excValue: object | None = None, traceback: object | None = None) -> bool: ...

    def SetProps(self, props: Sequence[str]) -> None:
        """
        Sets the properties to be written to the output file

          ARGUMENTS:

            - props: a list or tuple of property names
        """

    def write(self, mol: rdkit.Chem.rdchem.Mol, confId: int = -1) -> None:
        """
        Writes a molecule to the output file.

          ARGUMENTS:

            - mol: the Mol to be written
            - confId: (optional) ID of the conformation to write
        """

    def flush(self) -> None:
        """Flushes the output file (forces the disk file to be updated)."""

    def close(self) -> None:
        """
        Flushes the output file and closes it. The Writer cannot be used after this.
        """

    def NumMols(self) -> int:
        """Returns the number of molecules written so far."""

    def SetForceV3000(self, val: bool) -> None:
        """Sets whether or not V3000 mol file writing is being forced."""

    def GetForceV3000(self) -> bool:
        """Returns whether or not V3000 mol file writing is being forced."""

    def SetKekulize(self, val: bool) -> None:
        """Sets whether or not molecules are kekulized on writing."""

    def GetKekulize(self) -> bool:
        """Returns whether or not molecules are kekulized on writing."""

    @staticmethod
    def GetText(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, kekulize: bool = True, force_v3000: bool = False, molid: int = -1) -> str:
        """returns the SD text for a molecule"""

class TDTWriter:
    """A class for writing molecules to TDT files."""

    @overload
    def __init__(self, fileName: str | os.PathLike) -> None:
        """
        Constructor.

           If a string argument is provided, it will be treated as the name of the
           output file. If a file-like object is provided, output will be sent there.
        """

    @overload
    def __init__(self, fileObj: object) -> None: ...

    def __enter__(self) -> TDTWriter: ...

    def __exit__(self, excType: object | None = None, excValue: object | None = None, traceback: object | None = None) -> bool: ...

    def SetProps(self, props: Sequence[str]) -> None:
        """
        Sets the properties to be written to the output file

          ARGUMENTS:

            - props: a list or tuple of property names
        """

    def write(self, mol: rdkit.Chem.rdchem.Mol, confId: int = -1) -> None:
        """
        Writes a molecule to the output file.

          ARGUMENTS:

            - mol: the Mol to be written
            - confId: (optional) ID of the conformation to write
        """

    def flush(self) -> None:
        """Flushes the output file (forces the disk file to be updated)."""

    def close(self) -> None:
        """
        Flushes the output file and closes it. The Writer cannot be used after this.
        """

    def NumMols(self) -> int:
        """Returns the number of molecules written so far."""

    def SetWrite2D(self, state: bool = True) -> None:
        """Causes 2D conformations to be written (default is 3D conformations)."""

    def GetWrite2D(self) -> bool: ...

    def SetWriteNames(self, state: bool = True) -> None:
        """Causes names to be written to the output file as NAME records."""

    def GetWriteNames(self) -> bool: ...

    def SetNumDigits(self, numDigits: int) -> None:
        """Sets the number of digits to be written for coordinates."""

    def GetNumDigits(self) -> int: ...

class PDBWriter:
    """A class for writing molecules to PDB files."""

    @overload
    def __init__(self, fileName: str | os.PathLike, flavor: int = 0) -> None:
        """
        Constructor.

          ARGUMENTS:

            - fileName: name of the output file. ('-' to write to stdout)
            - flavor: (optional)
        """

    @overload
    def __init__(self, fileObj: object, flavor: int = 0) -> None: ...

    def __enter__(self) -> PDBWriter: ...

    def __exit__(self, excType: object | None = None, excValue: object | None = None, traceback: object | None = None) -> bool: ...

    def write(self, mol: rdkit.Chem.rdchem.Mol, confId: int = -1) -> None:
        """
        Writes a molecule to the output file.

          ARGUMENTS:

            - mol: the Mol to be written
            - confId: (optional) ignored
        """

    def flush(self) -> None:
        """Flushes the output file (forces the disk file to be updated)."""

    def close(self) -> None:
        """
        Flushes the output file and closes it. The Writer cannot be used after this.
        """

    def NumMols(self) -> int:
        """Returns the number of molecules written so far."""

class MaeWriter:
    """
    An experimental class for writing molecules to Maestro files.

      Usage examples:

        1) writing to a named file:

           >>> writer = MaeWriter('out.mae')
           >>> for mol in list_of_mols:
           ...    writer.write(mol)

        2) writing to a file-like object:

           >>> import gzip
           >>> outf=gzip.open('out.mae.gz','wt+')
           >>> writer = MaeWriter(outf)
           >>> for mol in list_of_mols:
           ...   writer.write(mol)
           >>> writer.close()
           >>> outf.close()

      By default all non-private molecule, atom and bond properties are written
      to the Maestro file. This can be changed using the SetProps method:

           >>> writer = MaeWriter('out.mae')
           >>> writer.SetProps(['prop1','prop2'])

      Properties that are specified, but are not present will be ignored.

      Kekulization is mandatory, as the Maestro format does not have
      the concept of an aromatic bond.

      As this is an experimental writer, many features are not supported yet,
      e.g. chirality and bond stereo labels, stereo groups, substance groups,
      isotopes, or even dummy atoms. Note that these are not supported by
      MaeMolSupplier either.
    """

    @overload
    def __init__(self, filename: str | os.PathLike) -> None: ...

    @overload
    def __init__(self, fileobj: object) -> None: ...

    def __enter__(self) -> MaeWriter: ...

    def __exit__(self, excType: object | None = None, excValue: object | None = None, traceback: object | None = None) -> bool: ...

    def SetProps(self, props_list: Sequence[str]) -> None:
        """
        Sets the atom and molecule properties to be written to the output file.

          ARGUMENTS:

            - props_list: a list of atom and molecule property names
        """

    def write(self, mol: rdkit.Chem.rdchem.Mol, confId: int = -1) -> None:
        """
        Writes a molecule to the output file.

          ARGUMENTS:

            - mol: the molecule to be written
            - confId: (optional) ID of the conformation to write
        """

    def flush(self) -> None:
        """Flushes the output file (forces the disk file to be updated)."""

    def close(self) -> None:
        """
        Flushes the output file and closes it. The writer cannot be used after this.
        """

    def NumMols(self) -> int:
        """Returns the number of molecules written so far."""

    @staticmethod
    def GetText(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, props_list: Sequence[str] = []) -> str:
        """Returns the Maestro CT block text for a molecule."""

class MultithreadedSmilesMolSupplier:
    """
    A class which concurrently supplies molecules from a text file.
      Please note that this class is still a bit experimental and the API may
      change in future releases.

      Usage examples:

        1) Lazy evaluation: the molecules might not be constructed until we ask for them:

           >>> suppl = MultithreadedSmilesMolSupplier('in.smi')
           >>> for mol in suppl:
           ...    if(mol):
           ...      mol.GetNumAtoms()

        2) Lazy evaluation 2:

           >>> suppl = MultithreadedSmilesMolSupplier('in.smi')
           >>> while (!suppl.atEnd()):
           ...    mol = next(mol)
           ...    if(mol):
           ...      mol.GetNumAtoms()
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, fileName: str | os.PathLike, delimiter: str = ' \t', smilesColumn: int = 0, nameColumn: int = 1, titleLine: bool = True, sanitize: bool = True, numWriterThreads: int = 1, sizeInputQueue: int = 5, sizeOutputQueue: int = 5) -> None:
        """
        Constructor

          ARGUMENTS:

            - fileName: name of the file to be read

            - delimiter: (optional) text delimiter (a string). Defaults to ' \\t'.

            - smilesColumn: (optional) index of the column containing the SMILES
              data. Defaults to 0.

            - nameColumn: (optional) index of the column containing molecule names.
              Defaults to 1.

            - titleLine: (optional) set this toggle if the file contains a title line.
              Defaults to true.

            - sanitize: (optional) toggles sanitization of molecules as they are read.
              Defaults to true.

            - numWriterThreads: (optional) number of writer threads. Defaults to 1.

            - sizeInputQueue: (optional) size of input/reader queue. Defaults to 5.

            - sizeOutputQueue: (optional) size of output/writer queue. Defaults to 5.
        """

    def __iter__(self) -> MultithreadedSmilesMolSupplier: ...

    def __enter__(self) -> MultithreadedSmilesMolSupplier: ...

    def __exit__(self, excType: object | None = None, excValue: object | None = None, traceback: object | None = None) -> bool: ...

    def __next__(self) -> rdkit.Chem.rdchem.Mol | None:
        """Returns the next molecule in the file. Raises _StopIteration_ on EOF."""

    def atEnd(self) -> bool:
        """Returns true if we have read all records else false."""

    def GetLastRecordId(self) -> int:
        """Returns the record id for the last extracted item."""

    def GetLastItemText(self) -> str:
        """Returns the text for the last extracted item."""

class MultithreadedSDMolSupplier:
    """
    A class which concurrently supplies molecules from an SD file.
      Please note that this class is still a bit experimental and the API may
      change in future releases.

      Usage examples:

        1) Lazy evaluation: the molecules might not be constructed until we ask for them:

           >>> suppl = MultithreadedSDMolSupplier('in.sdf')
           >>> for mol in suppl:
           ...    if(mol):
           ...      mol.GetNumAtoms()

        2) Lazy evaluation 2:

           >>> suppl = MultithreadedSDMolSupplier('in.sdf')
           >>> while (!suppl.atEnd()):
           ...    mol = next(mol)
           ...    if(mol):
           ...      mol.GetNumAtoms()
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, fileName: str | os.PathLike, sanitize: bool = True, removeHs: bool = True, strictParsing: bool = True, numWriterThreads: int = 1, sizeInputQueue: int = 5, sizeOutputQueue: int = 5) -> None:
        """
        Constructor

          ARGUMENTS:

            - fileName: name of the file to be read

            - sanitize: (optional) toggles sanitization of molecules as they are read.
              Defaults to true.

            - removeHs: (optional) removes Hs. Defaults to true.

            - strictParsing: (optional) allows strict or lax parsing. Defaults to true.

            - numWriterThreads: (optional) number of writer threads. Defaults to 1.

            - sizeInputQueue: (optional) size of input/reader queue. Defaults to 5.

            - sizeOutputQueue: (optional) size of output/writer queue. Defaults to 5.
        """

    def __iter__(self) -> MultithreadedSDMolSupplier: ...

    def __enter__(self) -> MultithreadedSDMolSupplier: ...

    def __exit__(self, excType: object | None = None, excValue: object | None = None, traceback: object | None = None) -> bool: ...

    def __next__(self) -> rdkit.Chem.rdchem.Mol | None:
        """Returns the next molecule in the file. Raises _StopIteration_ on EOF."""

    def atEnd(self) -> bool:
        """Returns true if we have read all records else false."""

    def GetLastRecordId(self) -> int:
        """Returns the record id for the last extracted item."""

    def GetLastItemText(self) -> str:
        """Returns the text for the last extracted item."""

    def GetProcessPropertyLists(self) -> bool:
        """
        Returns whether or not any property lists that are present will be processed when reading molecules.
        """

    def SetProcessPropertyLists(self, val: bool) -> None:
        """
        Sets whether or not any property lists that are present will be processed when reading molecules.
        """
