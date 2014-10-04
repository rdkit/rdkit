from __future__ import print_function
from Chem import rdmol
from Chem.rdmol import Atom,Bond,Mol


def ParseAtomBlock(lines,mol,nAtoms):
  for i in range(nAtoms):
    line = lines[i]
    pX = float(line[0:10])
    pY = float(line[10:20])
    pZ = float(line[20:30])
    symb = line[31:34].strip()
    newAt = rdmol.Atom(symb)
    newAt.setPos(pX,pY,pZ)

    chg = int(line[36:39])
    if chg in [1,2,3,5,6,7]:
      newAt.setFormalCharge(4-chg);

    # parse valence

    # parse rxn component

    mol.addAtom(newAt)

bondMap = {1:Bond.SINGLE,2:Bond.DOUBLE,3:Bond.TRIPLE,4:Bond.AROMATIC}

def ParseBondBlock(lines,mol,nBonds):
  for i in range(nBonds):
    line = lines[i]
    id1 = int(line[0:3])-1
    id2 = int(line[3:6])-1
    order = int(line[6:9])
    order = bondMap.get(order,Bond.OTHER);
    b = Bond(order)
    b.setOwningMol(mol)
    b.setBeginAtomIdx(id1)
    b.setEndAtomIdx(id2)
    mol.addBond(b)

def ParseMolBlock(lines,mol):
  header = lines[0:3]
  counts = lines[3]
  nAtoms = int(counts[0:3])
  nBonds = int(counts[3:6])
  nLists = int(counts[6:9])
  chiralFlag = int(counts[12:15])
  nsText = int(counts[15:18])
  nRxnComponents = int(counts[18:21])
  nReactants = int(counts[21:24])
  nProducts = int(counts[24:27])
  nIntermediates = int(counts[27:30])

  ParseAtomBlock(lines[4:],mol,nAtoms)
  ParseBondBlock(lines[4+nAtoms:],mol,nBonds)
  

if __name__ == '__main__':
  import sys
  fName = sys.argv[1]
  inF = open(fName,'r')
  lines = inF.readlines()
  m = rdmol.Mol()
  ParseMolBlock(lines,m)

  print(m.getNumAtoms())
  m.debugMol()

  print(rdmol.MolToCDXML(m))
