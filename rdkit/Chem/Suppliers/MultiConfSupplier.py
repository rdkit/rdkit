"""
Wrapper function to parse SDFs containing multiple conformers
"""

from rdkit import Chem


def CanonSmilesMap(mol):
    """
    Helper function to produce canonicalised SMILES and atom maps independent of initial atom order
    """

    for atm in mol.GetAtoms():
        if atm.GetAtomMapNum() != 0:
            atm.SetAtomMapNum(0)

    # canonicalise
    smiles = Chem.MolToSmiles(mol, isomericSmiles=False)
    mol = Chem.MolFromSmiles(smiles)

    for i, atm in enumerate(mol.GetAtoms()):
        atm.SetAtomMapNum(i + 1)

    return Chem.MolToSmiles(mol)


def MultiConfSupplier(filename, propertyName='SMILES'):
    """
        reads SDF file containing multiple conformers per molecule

    Adapted from https://greglandrum.github.io/rdkit-blog/posts/2022-10-28-dealing-with-multiconformer-sd-files.html

    ARGUMENTS:

        - filename: string with name of file to read (including extension)
        - propertyName: string with the property to group by. Options: "SMILES", Count line block ("countLine")
      or any attribute in file retrievable using mol.GetProp() e.g. name, ID

    RETURNS:
        generator object of molecules with cids assigned

    LoadMultiConfSDF useage:

    >>> suppl = MultiConfSupplier('test_confs.sdf')

    >>> mols = []
    >>> for mol in suppl:
    >>> mols.append(mol)

    """

    supplier = Chem.ForwardSDMolSupplier(filename)

    mol = None
    for itm in supplier:
        if itm is None:
            continue
        if not itm.GetConformer().Is3D():
            raise ValueError('Molecule has no 3D coordinates')
        if mol is None:
            mol = itm
            if propertyName == 'SMILES':
                refVal = CanonSmilesMap(mol)
            elif propertyName == 'countLine':
                refVal = Chem.MolToMolBlock(mol).split("\n")[3]
            else:
                refVal = mol.GetProp(propertyName)
            continue
        if propertyName == 'SMILES':
            pVal = CanonSmilesMap(itm)
        elif propertyName == 'countLine':
            pVal = Chem.MolToMolBlock(itm).split("\n")[3]
        else:
            pVal = itm.GetProp(propertyName)
        if pVal == refVal:
            mol.AddConformer(itm.GetConformer(), assignId=True)
        else:
            # we're done with the last molecule, so let's restart the next one
            res = mol
            mol = itm
            refVal = pVal
            yield res

    yield mol