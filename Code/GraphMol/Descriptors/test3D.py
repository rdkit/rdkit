from rdkit import Chem
from rdkit import rdBase

from rdkit.Chem import rdMolDescriptors as rdMD
from rdkit.Chem import AllChem

print rdBase.rdkitVersion

m = Chem.MolFromSmiles('Cc1ccccc1')

m = Chem.AddHs(m)

AllChem.EmbedMolecule(m)
AllChem.MMFFOptimizeMolecule(m)

print "1"
print rdMD.CalcRadiusOfGyration(m)
print "2"
print rdMD.calcRDFs(m)
print "3"
print rdMD.calcMORSEs(m)
print rdMD.calcWHIMs(m)
print rdMD.calcGETAWAYs(m)
print rdMD.calcAUTOCORR3Ds(m)
