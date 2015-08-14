from rdkit import Chem
from rdkit.Chem import AllChem

suppl = Chem.SmilesMolSupplier('../Data/chembl_20_chiral.smi')
for i,mol in enumerate(suppl):
    csmi = Chem.MolToSmiles(mol,True)
    for j in range(100):
        mh = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mh,randomSeed=j+1)
        Chem.AssignAtomChiralTagsFromStructure(mh)
        nm = Chem.RemoveHs(mh)
        smi = Chem.MolToSmiles(nm,True)
        if smi!=csmi:
            print('%d %d %s:\n%s\n%s'%(i,j,mol.GetProp('_Name'),csmi,smi))
            
