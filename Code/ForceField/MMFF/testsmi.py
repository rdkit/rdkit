from rdkit import Chem
from rdkit.Chem import AllChem

w = Chem.SDWriter('POP02A_smi.sdf')
mol = Chem.MolFromSmiles('O=[PH3]')
molH = Chem.AddHs(mol)
AllChem.EmbedMolecule(molH)
w.write(molH)
