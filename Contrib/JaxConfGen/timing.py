from rdkit import Chem
from rdkit.Chem import AllChem

import forcefield 
import numpy as np
import jax
import functools
import dg
import time
import optimizer

# aspirin = "O=C(C)Oc1ccccc1C(=O)O"
aspirin = "CCCCN(CCCC)C(=O)c1nn(c(C)c1Cl)-c1ccc(cc1C(=O)N1CCc2ccccc2C1)C(=O)NS(=O)(=O)c1ccc2ccc(I)cc2c1"
mol = Chem.MolFromSmiles(aspirin)

n_confs = 1000

start_time = time.time()
dg.generate_conformer(mol, n_confs, dims=3)
end_time = time.time()

# custom timing 6.110635757446289
print("custom timing", end_time - start_time)

s_time = time.time()
AllChem.EmbedMultipleConfs(mol, numConfs=n_confs)
e_time = time.time()

# rdkit timing 44.17475509643555
print("rdkit timing", e_time - s_time)