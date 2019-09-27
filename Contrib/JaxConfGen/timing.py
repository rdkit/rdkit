from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import Geometry

import forcefield 
import numpy as np
import jax
import functools
import dg
import time
import optimizer

# mol = "O=C(C)Oc1ccccc1C(=O)O"
mol = "CCCCN(CCCC)C(=O)c1nn(c(C)c1Cl)-c1ccc(cc1C(=O)N1CCc2ccccc2C1)C(=O)NS(=O)(=O)c1ccc2ccc(I)cc2c1"
# # mol = "c1ccccc1"
mol = Chem.MolFromSmiles(mol)
mol = Chem.AddHs(mol)

n_confs = 1000

start_time = time.time()
minimized_conf = dg.generate_conformer(mol, n_confs, dims=3)
end_time = time.time()
print("custom timing", end_time - start_time)
# min_pos = np.argmin(minimized)
# min_conf = 
print(type(minimized_conf))
print(minimized_conf.shape)
print(np.asarray(minimized_conf[0]))
conf = Chem.Conformer(mol.GetNumAtoms())
for idx, pos in enumerate(minimized_conf[0]):
    # point = AllChem.Point3D()
    # print(idx, (pos[0], pos[1], pos[2]))
    # print(idx, (type(pos[0])))
    # conf.SetAtomPosition(idx, Geometry.Point3D(pos[0], pos[1], pos[2]))
    conf.SetAtomPosition(idx, (float(pos[0]), float(pos[1]), float(pos[2])))

print(type(conf))
mol.AddConformer(conf)


print(Chem.MolToMolBlock(mol))
# minimized_conf[0]

# custom timing 6.110635757446289

mol.RemoveAllConformers()
s_time = time.time()
AllChem.EmbedMultipleConfs(mol, numConfs=n_confs)
e_time = time.time()

# rdkit timing 44.17475509643555
print("rdkit timing", e_time - s_time)
print(Chem.MolToMolBlock(mol))