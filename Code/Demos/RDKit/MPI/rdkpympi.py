# $Id$
#
# Copyright (C) 2009 Greg Landrum
#  All rights reserved
#
# Demo for using boost.mpi with the RDKit
#
# run this with : mpirun -n 4 python rdkpympi.py
#
from __future__ import print_function
from boost import mpi
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.RDLogger import logger
logger = logger()

def dividetask(data,task,silent=True):
    data=mpi.broadcast(mpi.world,data,0)

    nProcs = mpi.world.size
    chunkSize=len(data)//nProcs
    extraBits =len(data)%nProcs

    res=[]
    allRes=[]
    # the root node handles the extra pieces:
    if mpi.world.rank == 0:
        for i in range(extraBits):
          elem=data[i]
          res.append(task(elem))
          if not silent:
              logger.info('task(%d) done %d'%(mpi.world.rank,i+1))
    pos=extraBits+mpi.world.rank*chunkSize;
    for i in range(chunkSize):
        elem=data[pos]
        pos += 1
        res.append(task(elem))
        if not silent:
            logger.info('task(%d) done %d'%(mpi.world.rank,i+1))
    if mpi.world.rank==0:
        tmp=mpi.gather(mpi.world,res,0)
        for res in tmp: allRes.extend(res)
    else:
        mpi.gather(mpi.world,res,0)
    return allRes
    
if __name__=='__main__':
    from rdkit import RDConfig
    import os
    fName = os.path.join(RDConfig.RDBaseDir,'Projects','DbCLI','testData','bzr.sdf')
    if mpi.world.rank==0:
        data = [x for x in Chem.SDMolSupplier(fName)][:50]
    else:
        data=None

    def generateconformations(m):
        m = Chem.AddHs(m)
        ids=AllChem.EmbedMultipleConfs(m,numConfs=10)
        for id in ids:
            AllChem.UFFOptimizeMolecule(m,confId=id)
        return m

    cms=dividetask(data,generateconformations,silent=False)
    # report:
    if mpi.world.rank==0:
        for i,mol in enumerate(cms):
            print(i,mol.GetNumConformers())
