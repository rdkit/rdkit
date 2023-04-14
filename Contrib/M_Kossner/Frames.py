#!/usr/bin/python
# encoding: utf-8

#	Jan 2011	(markus kossner)	Cleaned up the code, added some documentation
#	somewhere around Aug 2008	(markus kossner)	created
#
#    This script extracts the molecular framework for a database of molecules.
#    You can use two modes (hard coded):
#    - Scaff:	The molecular frame is extracted
#    - RedScaff:	All linking chains between rings are deleted. The rings are directly connected.
#
#    You can comment in/out the code snippets indicated by the comments
#    to force each atom of the frame to be a Carbon.
#
#    Usage: Frames.py <database.sdf>
#    Output:
#    - sd files containing all molecules belonging to one frame (1.sdf, 2.sdf etc)
#    - frames.smi containing the (canonical) smiles and count of occurrence
#

import os
import sys

from Chem import AllChem as Chem


def flatten(x):
  """flatten(sequence) -> list
    Returns a single, flat list which contains all elements retrieved
    from the sequence and all nested sub-sequences (iterables).
    Examples:
    >>> [1, 2, [3,4], (5,6)]
    [1, 2, [3, 4], (5, 6)]
    >>> flatten([[[1,2,3], (42,None)], [4,5], [6], 7, MyVector(8,9,10)])
    [1, 2, 3, 42, None, 4, 5, 6, 7, 8, 9, 10]"""
  result = []
  for el in x:
    if hasattr(el, "__iter__") and not isinstance(el, str):
      result.extend(flatten(el))
    else:
      result.append(el)
  return result


def GetFrame(mol, mode='Scaff'):
  '''return a ganeric molecule defining the reduced scaffold of the input mol.
    mode can be 'Scaff' or 'RedScaff':

    Scaff	->	chop off the side chains and return the scaffold

    RedScaff	->	remove all linking chains and connect the rings
    directly at the atoms where the linker was
    '''

  ring = mol.GetRingInfo()
  RingAtoms = flatten(ring.AtomRings())
  NonRingAtoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIdx() not in RingAtoms]
  RingNeighbors = []
  Paths = []
  for NonRingAtom in NonRingAtoms:
    for neighbor in mol.GetAtomWithIdx(NonRingAtom).GetNeighbors():
      if neighbor.GetIdx() in RingAtoms:
        RingNeighbors.append(NonRingAtom)
        #The ring Atoms having a non ring Neighbor will be the start of a walk
        Paths.append([neighbor.GetIdx(), NonRingAtom])
        break
  PosConnectors = [x for x in NonRingAtoms if x not in RingNeighbors
                   ]  #Only these Atoms are potential starting points of a Linker chain
  #print 'PosConnectors:'
  #print PosConnectors
  Framework = [x for x in RingAtoms]
  #Start a list of pathways which we will have to walk
  #print 'Path atoms:'
  #print Paths
  Linkers = []
  while len(Paths) > 0:
    NewPaths = []
    for P in Paths:
      if P is None:
        print('ooh')
      else:
        for neighbor in mol.GetAtomWithIdx(P[-1]).GetNeighbors():
          if neighbor.GetIdx() not in P:
            if neighbor.GetIdx() in NonRingAtoms:
              n = P[:]
              n.append(neighbor.GetIdx())
              NewPaths.append(n[:])
            elif neighbor.GetIdx() in RingAtoms:
              #print 'adding the following path to Framework:'
              #print P
              n = P[:]
              n.append(neighbor.GetIdx())
              Linkers.append(n)
              Framework = Framework + P[:]

    Paths = NewPaths[:]
  #print 'Linkers:',Linkers
  #print 'RingAtoms:',RingAtoms
  #em.AddBond(3,4,Chem.BondType.SINGLE)
  if mode == 'RedScaff':
    Framework = list(set(Framework))
    todel = []
    NonRingAtoms.sort(reverse=True)
    em = Chem.EditableMol(mol)
    BondsToAdd = [sorted([i[0], i[-1]]) for i in Linkers]
    mem = []
    for i in BondsToAdd:
      if i not in mem:
        em.AddBond(i[0], i[1], Chem.BondType.SINGLE)
        mem.append(i)
    for i in NonRingAtoms:
      todel.append(i)
    for i in todel:
      em.RemoveAtom(i)
    m = em.GetMol()
    #===================================#
    #  Now do the flattening of atoms and bonds!
    #  Any heavy atom will become a carbon and any bond will become a single bond!	#
    #===================================#
    #		for atom in m.GetAtoms():                                                 #
    #			atom.SetAtomicNum(6)                                                    #
    #			atom.SetFormalCharge(0)                                                #
    #		for bond in m.GetBonds():                                                   #
    #			bond.SetBondType(Chem.BondType.SINGLE)                 #
    #		Chem.SanitizeMol(m)                                                          #
    #===================================#
    return m

  if mode == 'Scaff':
    Framework = list(set(Framework))
    todel = []
    NonRingAtoms.sort(reverse=True)
    for i in NonRingAtoms:
      if i is not None:
        if i not in Framework:
          todel.append(i)
    em = Chem.EditableMol(mol)
    for i in todel:
      em.RemoveAtom(i)
    m = em.GetMol()
    #===================================#
    #  Now do the flattening of atoms and bonds!
    #  Any heavy atom will become a carbon and any bond will become a single bond!!		#
    #===================================#
    #		for atom in m.GetAtoms():                                                 #
    #			atom.SetAtomicNum(6)                                                    #
    #			atom.SetFormalCharge(0)                                                #
    #		for bond in m.GetBonds():                                                   #
    #			bond.SetBondType(Chem.BondType.SINGLE)                 #
    #		Chem.SanitizeMol(m)                                                          #
    #===================================#
    return m


if __name__ == '__main__':
  if len(sys.argv) < 2:
    print("No input file provided: Frames.py filetosprocess.ext")
    sys.exit(1)

  suppl = Chem.SDMolSupplier(sys.argv[1])
  FrameDict = {}

  for mol in suppl:
    m = GetFrame(mol)
    cansmiles = Chem.MolToSmiles(m, isomericSmiles=True)
    if cansmiles in FrameDict:
      FrameDict[cansmiles].append(mol)
    else:
      FrameDict[cansmiles] = [
        mol,
      ]

  counter = 0
  w = open('frames.smi', 'w')
  for key, item in FrameDict.items():
    counter += 1
    d = Chem.SDWriter(str(counter) + '.sdf')
    for i in item:
      i.SetProp('Scaffold', key)
      i.SetProp('Cluster', str(counter))
      d.write(i)
    print(key, len(item))
    w.write(key + '\t' + str(len(item)) + '\n')
  w.close
  print('number of Clusters: %d' % (counter))
