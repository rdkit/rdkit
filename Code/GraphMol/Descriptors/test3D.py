from rdkit import Chem
from rdkit import rdBase

from rdkit.Chem import rdMolDescriptors as rdMD
from rdkit.Chem import AllChem

print rdBase.rdkitVersion

m = Chem.MolFromSmiles('Cc1ccccc1')

m = Chem.AddHs(m)

AllChem.EmbedMolecule(m)
AllChem.MMFFOptimizeMolecule(m)

print dir(rdMD)


print "1"
print rdMD.CalcRDF(m)
print "2"
print rdMD.CalcMORSE(m)
print "3"
print rdMD.CalcWHIM(m)
print "4"
print rdMD.CalcGETAWAY(m)
print "5"
print rdMD.CalcAUTOCORR3D(m)
