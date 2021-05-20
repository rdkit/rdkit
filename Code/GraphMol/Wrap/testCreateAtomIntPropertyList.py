#%%
import rdkit
print(rdkit.__file__)
print(rdkit.__version__)
#%%

from rdkit import Chem

m = Chem.MolFromSmiles("CC")
for a in m.GetAtoms():
    a.SetIntProp("foo", 1)

Chem.CreateAtomIntPropertyList(m, "foo")

print(m.GetProp("atom.iprop.foo"))
# output: 1 1
# %%
