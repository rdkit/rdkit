from rdkit import Chem
from rdkit.Chem import AllChem

w = Chem.SDWriter('MMFF94_hypervalent_smi.sdf')
suppl = Chem.SmilesMolSupplier('MMFF94_hypervalent.smi')
for mol in suppl:
  molH = Chem.AddHs(mol)
  AllChem.EmbedMolecule(molH)
  w.write(molH)

