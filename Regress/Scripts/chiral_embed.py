from rdkit import Chem
from rdkit.Chem import AllChem
import gzip

for i,line in enumerate(gzip.open('../Data/chembl_20_chiral.smi.gz')):
    line = line.strip().decode().split(' ')
    mol = Chem.MolFromSmiles(line[0])
    if not mol:
        continue
    nm = line[1]
    csmi = Chem.MolToSmiles(mol,True)
    for j in range(100):
        mh = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mh,randomSeed=j+1)
        Chem.AssignAtomChiralTagsFromStructure(mh)
        newm = Chem.RemoveHs(mh)
        smi = Chem.MolToSmiles(newm,True)
        if smi!=csmi:
            print('%d %d %s:\n%s\n%s'%(i,j,nm,csmi,smi))
            
