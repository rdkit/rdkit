"""
Wrapper function to parse SDFs containing multiple conformers
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def MultiConfSupplier(filename, propertyName=None, includeStereo=False):
    """
    Reads SDF file containing multiple conformers per molecule

    Adapted from https://greglandrum.github.io/rdkit-blog/posts/2022-10-28-dealing-with-multiconformer-sd-files.html

    ARGUMENTS:

        - filename: string with name of file to read
        - propertyName: string with the property to group by.
          Options:
            - countLine: checks Count line at top of block matches
            - SMILES
            - property: any attribute retrievable using mol.GetProp() (name, ID etc)
        - includeStereo: consider stereoisomers separately when grouping conformers (default False)

    RETURNS:
        generator object of molecules with cids assigned

    LoadMultiConfSDF useage:

    >>> suppl = MultiConfSupplier('test_confs.sdf')

    >>> mols = []
    >>> for mol in suppl:
    >>> mols.append(mol)

    """

    def CalcProp(m, propertyName):
        """
        Helper function to compute property to compare
        """

        if propertyName == 'countLine':

            prop = Chem.MolToMolBlock(m).split("\n")[3]

        elif propertyName == 'SMILES':

            for atm in m.GetAtoms():
                if atm.GetAtomMapNum() != 0:  # goes weird if already 0
                    atm.SetAtomMapNum(0)

                # canonicalise
                smiles = Chem.MolToSmiles(m, isomericSmiles=False)
                m = Chem.MolFromSmiles(smiles)

                for i, atm in enumerate(m.GetAtoms()):
                    atm.SetAtomMapNum(i + 1)

                prop = Chem.MolToSmiles(m)

        else:
            prop = m.GetProp(propertyName)

        return prop

    def matchConf(a, b, propertyName, includeStereo):
        """
        Helper function to check for conformers
        """

        # fail fast checks

        # number of atoms and bonds
        if rdMolDescriptors.CalcMolFormula(a) != rdMolDescriptors.CalcMolFormula(
                b) or a.GetNumAtoms() != b.GetNumAtoms() or a.GetNumBonds() != b.GetNumBonds():
            return False

        # atom type
        if [atom.GetSymbol() for atom in a.GetAtoms()] != [atom.GetSymbol() for atom in b.GetAtoms()]:
            return False

        # atom stereo
        if includeStereo:
            if [atom.GetChiralTag() for atom in a.GetAtoms()] != [atom.GetChiralTag() for atom in b.GetAtoms()]:
                return False

        # bond order
        if [bond.GetBondType() for bond in a.GetBonds()] != [bond.GetBondType() for bond in b.GetBonds()]:
            return False

        # bond stereo
        if [bond.GetStereo() for bond in a.GetBonds()] != [bond.GetStereo() for bond in b.GetBonds()]:
            return False

        # bonds by index - check first and last are same
        if a.GetBonds()[0].GetIdx() != b.GetBonds()[0].GetIdx() or a.GetBonds()[-1].GetIdx() != b.GetBonds()[
            -1].GetIdx(): \
                return False

        # if all pass, check supplied property
        if propertyName is not None:
            return CalcProp(a, propertyName) == CalcProp(b, propertyName)

        return True

    supplier = Chem.ForwardSDMolSupplier(filename)

    mol = None
    for itm in supplier:
        if itm is None:
            continue
        if mol is None:
            mol = itm
            continue

        if matchConf(mol, itm, propertyName, includeStereo):
            mol.AddConformer(itm.GetConformer(), assignId=True)
        else:
            # we're done with the last molecule, so let's restart the next one
            res = mol
            mol = itm
            yield res

    yield mol