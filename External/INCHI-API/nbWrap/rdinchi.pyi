import rdkit.Chem.rdchem


def InchiToMol(inchi: str, sanitize: bool = True, removeHs: bool = True) -> tuple[rdkit.Chem.rdchem.Mol | None, int, str, str]:
    """
    return a ROMol for a InChI string
      Returns:
        a tuple with:
          - the molecule
          - the return code from the InChI conversion
          - a string with any messages from the InChI conversion
          - a string with any log messages from the InChI conversion
    """

def MolToInchi(mol: rdkit.Chem.rdchem.Mol, options: str = '') -> tuple[str, int, str, str, str]:
    """
    return the InChI for a ROMol molecule.

      Arguments:
        - mol: the molecule to use.
        - options: the InChI generation options.
          Options should be prefixed with either a - or a /
          Available options are explained in the InChI technical FAQ:
          https://www.inchi-trust.org/technical-faq/#15.14
          and the User Guide:
          https://github.com/IUPAC-InChI/InChI/blob/main/INCHI-1-DOC/UserGuide/InChI_UserGuide.pdf
      Returns:
        a tuple with:
          - the InChI
          - the return code from the InChI conversion
          - a string with any messages from the InChI conversion
          - a string with any log messages from the InChI conversion
          - a string with the InChI AuxInfo
    """

def MolBlockToInchi(molblock: str, options: str = '') -> tuple[str, int, str, str, str]:
    """
    return the InChI for a ROMol molecule.

      Arguments:
        - molblock: the mol block to use.
        - options: the InChI generation options.
          Options should be prefixed with either a - or a /
          Available options are explained in the InChI technical FAQ:
          https://www.inchi-trust.org/technical-faq/#15.14
          and the User Guide:
          https://github.com/IUPAC-InChI/InChI/blob/main/INCHI-1-DOC/UserGuide/InChI_UserGuide.pdf
      Returns:
        a tuple with:
          - the InChI
          - the return code from the InChI conversion
          - a string with any messages from the InChI conversion
          - a string with any log messages from the InChI conversion
          - a string with the InChI AuxInfo
    """

def InchiToInchiKey(inchi: str) -> str:
    """return the InChI key for an InChI string"""

def MolToInchiKey(mol: rdkit.Chem.rdchem.Mol, options: str = '') -> str:
    """
    return the InChI key for a ROMol molecule.

      Arguments:
        - mol: the molecule to use.
        - options: the InChI generation options.
          Options should be prefixed with either a - or a /
          Available options are explained in the InChI technical FAQ:
          https://www.inchi-trust.org/technical-faq/#15.14
          and the User Guide available from:
          https://github.com/IUPAC-InChI/InChI/blob/main/INCHI-1-DOC/UserGuide/InChI_UserGuide.pdf
      Returns: the InChI key
    """

def GetInchiVersion() -> str:
    """returns the version of the InChI software being used"""
