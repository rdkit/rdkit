#!/usr/bin/env python

import rdkit
import rdkit.Chem.rdShape
#print rdkit.Chem.rdShape.Align.__doc__
from rdkit import Chem
mol1 = Chem.SDMolSupplier("reference.sdf").next()
supp = Chem.SDMolSupplier("database.sdf")
mol2 = supp.next()
score = rdkit.Chem.rdShape.Align(mol1, mol2, maxIter=10)
print score

supp = Chem.SDMolSupplier("database.sdf")
tanimoto = [rdkit.Chem.rdShape.Align(mol1, m, maxIter=0) for m in supp]


import pandas as pd
df1 = pd.read_csv("output_rdkit_1.tab", sep="\t")
df2 = pd.read_csv("output_rdkit_2.tab", sep="\t")

df = pd.DataFrame(tanimoto, columns=["Shape-it::Tanimoto",])

import numpy as np
np.testing.assert_array_almost_equal(df["Shape-it::Tanimoto"], df1["Shape-it::Tanimoto"], decimal=3)
np.testing.assert_array_almost_equal(df["Shape-it::Tanimoto"], df2["Shape-it::Tanimoto"], decimal=3)

