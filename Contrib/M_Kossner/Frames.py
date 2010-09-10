#!/usr/bin/env python2.6
# encoding: utf-8
# $Id$
# original author: Markus Kossner

import os,sys
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
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result

def GetFrame(mol,mode='RedScaff'):
	'''return a ganeric molecule defining the reduced scaffold of the molecule.
	mode can be 'RedScaff' and 'Scaff'. The second mode will turn every atom to a carbon and every bond to a single bond!'''
	ring=mol.GetRingInfo()
	RingAtoms= flatten(ring.AtomRings())
	NonRingAtoms=[ atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIdx() not in RingAtoms ]
	RingNeighbors=[]
	Paths=[]
	for NonRingAtom in NonRingAtoms:
		for neighbor in mol.GetAtomWithIdx(NonRingAtom).GetNeighbors():
			if neighbor.GetIdx() in RingAtoms:
				RingNeighbors.append(NonRingAtom)
				Paths.append([neighbor.GetIdx(),NonRingAtom]) #The ring Atoms having a non ring Nieghbor will be the start of a walk
				break
	PosLinkers=[x for x in NonRingAtoms if x not in RingNeighbors] #Only these Atoms are Possible Linkers of two rings
	Framework=[ x for x in RingAtoms ]
	Linkers=[]
	while len(Paths)>0:
		NewPaths=[]
		for P in Paths:
			if P==None:
				print 'ooh, there is still a bug somewere '
			else:
				for neighbor in mol.GetAtomWithIdx(P[-1]).GetNeighbors():
					if neighbor.GetIdx() not in P:
						if neighbor.GetIdx() in NonRingAtoms:
							n=P[:]
							n.append(neighbor.GetIdx())
							NewPaths.append(n[:])
						elif neighbor.GetIdx() in RingAtoms:
							n=P[:]
							n.append(neighbor.GetIdx())
							Linkers.append(n)
							Framework=Framework+P[:]

		Paths=NewPaths[:]
	#print 'Linkers:',Linkers
	#print 'RingAtoms:',RingAtoms
	if mode=='Scaff':
		Framework=list(set(Framework))
		todel=[]
		NonRingAtoms.sort(reverse=True)
		em=Chem.EditableMol(mol)
		BondsToAdd=[ sorted([i[0],i[-1]]) for i in Linkers ]
		mem=[]
		for i in BondsToAdd:
			if i not in mem:
				em.AddBond(i[0],i[1],Chem.BondType.SINGLE)
				mem.append(i)
		for i in NonRingAtoms:
			todel.append(i)
		for i in todel:
			em.RemoveAtom(i)
		m=em.GetMol()
		return m

	if mode=='RedScaff':
		Framework=list(set(Framework))
		todel=[]
		NonRingAtoms.sort(reverse=True)
		for i in NonRingAtoms:
			if i != None:
				if i not in Framework:
					todel.append(i)
		em=Chem.EditableMol(mol)
		for i in todel:
			em.RemoveAtom(i)
		m=em.GetMol()
		#===================================#
		#Now do the flattening of atoms and bonds! 
		#Any heavy atom will become a carbon and any bond will become a single bond!!		#
		#===================================#
		for atom in m.GetAtoms():                                                 #
			atom.SetAtomicNum(6)                                                    #
			atom.SetFormalCharge(0)                                                #
		for bond in m.GetBonds():                                                   #
			bond.SetBondType(Chem.BondType.SINGLE)                 #
		Chem.SanitizeMol(m)                                                          #
		#===================================#
		return m

if len(sys.argv) < 2:
	print "No input file provided: Frames.py <file-to-process.sdf>"
	sys.exit(1)

suppl=Chem.SDMolSupplier(sys.argv[1])

FrameDict={}

for mol in suppl:
	if mol == 'None': continue
	try:
		m=GetFrame(mol,mode='RedScaff')
		if FrameDict.has_key(Chem.MolToSmiles(m)):
			FrameDict[Chem.MolToSmiles(m)].append(mol)
		else:
			FrameDict[Chem.MolToSmiles(m)]=[mol,]
	except:
		print '--------------------------------'
		print 'could not process the molecule with the following name:'
		print mol.GetProp('_Name')

counter=0
w=open('frames.smi','w')
d=Chem.SDWriter('frames.sdf')
for key,item in FrameDict.items():
	counter+=1
	#d=Chem.SDWriter(str(counter)+'.sdf')
	for i in item:
		i.SetProp('Scaffold',key)
		i.SetProp('Cluster',str(counter))
		d.write(i)
	print key,len(item)
	w.write(key+'\n')
w.close
print 'number of Frames found: %d' %(counter)
