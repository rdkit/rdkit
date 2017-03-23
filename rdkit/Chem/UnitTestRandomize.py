#
#  Copyright (C) 2017 Esben Jannik Bjerrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
from __future__ import print_function

import unittest

from rdkit import Chem
from rdkit.Chem import Randomize
import random


def MolBlockToSmiles(mb):
  return Chem.MolToSmiles(Chem.MolFromMolBlock(mb))

class TestCase(unittest.TestCase):

  def test1(self):
    random.seed(42)
    smiles = ["C[N+](C)(C)C", "CC(=O)[O-]", "[NH3+]CC(=O)[O-]",'[Cl-].[Cl-].[Mg+2]','[O-][O-]','[O-]C[O-]','[Na+].[Na+].[O-2]']
    msg = "Could not randomize charged mol %s"
    for i in range(10):
      for smi in smiles:
        mb = Chem.MolToMolBlock(Chem.MolFromSmiles(smi))
        rmb = Randomize.RandomizeMolBlock(mb)
        #Test if fails conversion
        self.assertNotEqual(None,Chem.MolFromMolBlock(rmb), msg=msg % (smi))
        #Test if canonical Smile is same
        self.assertEqual(MolBlockToSmiles(mb),MolBlockToSmiles(mb),msg=msg % (smi))

  def test2(self):
    random.seed(42)
    smiles = ['CON', 'c1ccccn1', 'C/C=C/F']
    msg = "\nRef: %s\n   : %s"
    nReps = 10
    for smi in smiles:
      mol = Chem.MolFromSmiles(smi)
      refSmi = Chem.MolToSmiles(mol, False)
      for i in range(nReps):
        m2 = Randomize.RandomizeMol(mol)
        randomizedsmi = Chem.MolToSmiles(m2, False)
        self.assertEqual(refSmi, randomizedsmi, msg=msg % (refSmi, smi))

if __name__ == '__main__':  # pragma: nocover
  unittest.main()
