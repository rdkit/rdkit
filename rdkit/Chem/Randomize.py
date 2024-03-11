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


def _shuffle(indices):
    # Shuffle without in-place mutation.
    # See https://docs.python.org/3/library/random.html#random.shuffle.
    return random.sample(indices, len(indices))


def RandomizeMolBlock(molblock):
    mol = Chem.MolFromMolBlock(molblock, sanitize=False, removeHs=False)
    atom_indices = [atom.GetIdx() for atom in mol.GetAtoms()]
    atom_indices_randomized = _shuffle(atom_indices)
    if len(atom_indices) > 1:
        # Enforce randomization.
        while atom_indices_randomized == atom_indices:
            atom_indices_randomized = _shuffle(atom_indices)

    return Chem.MolToMolBlock(Chem.RenumberAtoms(mol, atom_indices_randomized))


def RandomizeMol(mol):
    molblock = Chem.MolToMolBlock(mol)
    molblock_randomized = RandomizeMolBlock(molblock)

    return Chem.MolFromMolBlock(molblock_randomized)
