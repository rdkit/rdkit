# $Id$
#
#  Copyright (C) 2005-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import random
from rdkit import Chem


def RandomizeMolBlock(molB):
    splitB = molB.split('\n')
    res = []
    res.extend(splitB[0:3])
    idx = 3
    inL = splitB[idx]
    res.append(inL)
    nAts = int(inL[0:3])
    nBonds = int(inL[3:6])

    idx += 1
    atLines = splitB[idx:idx + nAts]

    order = list(range(nAts))
    random.shuffle(order, random=random.random)

    for i in order:
        res.append(atLines[i])

    # print 'ORDER:',order
    idx += nAts
    for i in range(nBonds):
        inL = splitB[idx]
        idx1 = int(inL[0:3]) - 1
        idx2 = int(inL[3:6]) - 1
        idx1 = order.index(idx1)
        idx2 = order.index(idx2)
        inL = '% 3d% 3d' % (idx1 + 1, idx2 + 1) + inL[6:]
        res.append(inL)
        idx += 1
    res.append('M  END')
    return '\n'.join(res)


def RandomizeMol(mol):
    mb = Chem.MolToMolBlock(mol)
    # print '-----------------'
    # print mb
    mb = RandomizeMolBlock(mb)
    # print mb
    return Chem.MolFromMolBlock(mb)


def CheckCanonicalization(mol, nReps=10):
    refSmi = Chem.MolToSmiles(mol, False)
    for i in range(nReps):
        m2 = RandomizeMol(mol)
        smi = Chem.MolToSmiles(m2, False)
        if smi != refSmi:
            raise ValueError('\nRef: %s\n   : %s' % (refSmi, smi))


if __name__ == '__main__':
    from rdkit.Chem import Randomize
    CheckCanonicalization(Chem.MolFromSmiles('CON'))
    CheckCanonicalization(Chem.MolFromSmiles('c1ccccn1'))
    CheckCanonicalization(Chem.MolFromSmiles('C/C=C/F'))
