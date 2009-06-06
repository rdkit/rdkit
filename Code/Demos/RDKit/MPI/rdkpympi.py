# $Id$
#
# Copyright (C) 2009 Greg Landrum
#  All rights reserved
#
# Demo for using boost.mpi with the RDKit
#
# run this with : mpirun -n 4 python rdkpympi.py
#
from boost import mpi
from rdkit import Chem
import sys

if mpi.world.rank==0:
    data = [Chem.MolFromSmiles('C'*x) for x in range(1,100)]
else:
    data=None
data=mpi.broadcast(mpi.world,data,0)


res=[]
allRes=[]
nProcs = mpi.world.size
chunkSize=len(data)//nProcs
extraBits =len(data)%nProcs


# handle extra bits on the root node:
if mpi.world.rank == 0:
    for i in range(extraBits):
      elem=data[i]
      res.append(elem.GetNumAtoms(False))
  
pos=extraBits+mpi.world.rank*chunkSize;
for i in range(chunkSize):
    elem=data[pos]
    pos += 1
    res.append(elem.GetNumAtoms(False))


if mpi.world.rank==0:
    allRes=mpi.gather(mpi.world,res,0)
else:
    mpi.gather(mpi.world,res,0)
    
# report:
if mpi.world.rank==0:
    for i in range(mpi.world.size):
        print "results from process:",i,": ",allRes[i]
