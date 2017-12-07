import argparse
from rdkit import Chem
sdf = Chem.SDMolSupplier("cdk2.sdf")
f = open("cdk2.smi","w")
for mol in sdf:
    name = mol.GetProp("_Name")
    smi = Chem.MolToSmiles( mol )
    f.write( "{}\t{}\n".format(name,smi))
f.close()

