import time

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors as rdMD
from rdkit.Chem.EState import AtomTypes, EStateIndices

print(rdBase.rdkitVersion)
print(rdBase.boostVersion)


def getEState(mol):

  return EStateIndices(mol)


def localopt(mol, steps=500):
  if mol.GetNumConformers() == 0:
    mol = make3D(mol)
  AllChem.MMFFOptimizeMolecule(mol, maxIters=steps)
  return mol


def make3D(mol, steps=50):
  mol = Chem.AddHs(mol)
  success = AllChem.EmbedMolecule(mol)
  if success == -1:  # Failed
    success = AllChem.EmbedMolecule(mol, useRandomCoords=True)
    if success == -1:
      raise (Error, "Embedding failed!")
  mol = localopt(mol, steps)
  return mol


def get3D(m, is3d):
  if not is3d:
    m = Chem.AddHs(m)
    AllChem.EmbedMolecule(m)
    AllChem.MMFFOptimizeMolecule(m)
  r = rdMD.CalcAUTOCORR3D(m) + rdMD.CalcRDF(m) + rdMD.CalcMORSE(m) + rdMD.CalcWHIM(
    m) + rdMD.CalcGETAWAY(m)
  return r


def generateALL():
  m = Chem.MolFromSmiles('Cc1ccccc1')
  thefile = open('testAC.txt', 'w')
  filename = "/Users/mbp/Github/rdkit_mine/Code/GraphMol/Descriptors/test_data/PBF_egfr.sdf"
  suppl = Chem.SDMolSupplier(filename, removeHs=False)
  mols = [x for x in suppl]
  start = time.time()
  for m in mols:
    r = get3D(m, True)
    for item in r:
      thefile.write("%.3f," % item)
    thefile.write("\n")

  end = time.time()
  print(end - start)


thefile = open('testSMWHIM.txt', 'w')
writer = Chem.SDWriter('3Dsmallmol.sdf')
A = [
  '[H][H]', 'B', 'O=O', 'C', 'CC', 'CCC', 'CCCC', 'CCCCC', 'CCCCCC', 'CO', 'CCO', 'CCCO', 'CCCCO',
  'CCCCCO', 'CCCCCCO', 'CCl', 'CCCl', 'CCCCl', 'CCCCCl', 'CCCCCCl', 'CCCCCCCl', 'CBr', 'CCBr',
  'CCCBr', 'CCCCBr', 'CCCCCBr', 'CCCCCCBr', 'CI', 'CCI', 'CCCI', 'CCCCI', 'CCCCCI', 'CCCCCCI', 'CF',
  'CCF', 'CCCF', 'CCCCF', 'CCCCCF', 'CCCCCCF', 'CS', 'CCS', 'CCCS', 'CCCCS', 'CCCCCS', 'CCCCCCS',
  'CN', 'CCN', 'CCCN', 'CCCCN', 'CCCCCN', 'CCCCCCN'
]
for smi in A:
  m = Chem.MolFromSmiles(smi)
  m = localopt(m, 100)
  #r=get3D(m,True)
  print(smi)
  print("---------")
  r = rdMD.CalcWHIM(m)
  print("Ei:" + str(r[0]) + "," + str(r[1]) + "," + str(r[2]) + "\n")
  print("Gi:" + str(r[5]) + "," + str(r[6]) + "," + str(r[7]) + "\n")
  print("SI:" + str(rdMD.CalcSpherocityIndex(m)))
  print("AS:" + str(rdMD.CalcAsphericity(m)))
  print("EX:" + str(rdMD.CalcEccentricity(m)))
  for item in r:
    thefile.write("%.3f," % item)
  thefile.write("\n")
  #m.SetProp("smi", smi)
  #writer.write(m)

thefile = open('testBPA.txt', 'w')
writer = Chem.SDWriter('3DBPAmol.sdf')
B = [
  'CN(C)CC(Br)c1ccccc1', 'CN(C)CC(Br)c1ccc(F)cc1', 'CN(C)CC(Br)c1ccc(Cl)cc1',
  'CN(C)CC(Br)c1ccc(Cl)cc1', 'CN(C)CC(Br)c1ccc(I)cc1', 'CN(C)CC(Br)c1ccc(C)cc1',
  'CN(C)CC(Br)c1cccc(F)c1', 'CN(C)CC(Br)c1cccc(Cl)c1', 'CN(C)CC(Br)c1cccc(Br)c1',
  'CN(C)CC(Br)c1cccc(I)c1', 'CN(C)CC(Br)c1cccc(C)c1', 'CN(C)CC(Br)c1ccc(F)c(Cl)c1',
  'CN(C)CC(Br)c1ccc(F)c(Br)c1', 'CN(C)CC(Br)c1ccc(F)c(C)c1', 'CN(C)CC(Br)c1ccc(Cl)c(Cl)c1',
  'CN(C)CC(Br)c1ccc(Cl)c(Br)c1', 'CN(C)CC(Br)c1ccc(Cl)c(C)c1', 'CN(C)CC(Br)c1ccc(Br)c(Cl)c1',
  'CN(C)CC(Br)c1ccc(Br)c(Br)c1', 'CN(C)CC(Br)c1ccc(Br)c(C)c1', 'CN(C)CC(Br)c1ccc(C)c(C)c1',
  'CN(C)CC(Br)c1ccc(C)c(Br)c1'
]
for smi in B:
  m = Chem.MolFromSmiles(smi)
  m = localopt(m, 100)
  #r=get3D(m,True)
  r = rdMD.CalcWHIM(m)
  for item in r:
    thefile.write("%.3f," % item)
  thefile.write("\n")
  #m.SetProp("smi", smi)
  #writer.write(m)

A = "G1w,G2w,G3w,Gw"
print(dir(rdMD))
