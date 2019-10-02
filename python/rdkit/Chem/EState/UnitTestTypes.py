# $Id$
#
#  Copyright (C) 2002-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""unit testing code for the EState atom typing

validation values are from the paper (JCICS _35_ 1039-1045 (1995))

"""

import unittest
from rdkit import Chem
from rdkit.Chem.EState import AtomTypes


class TestCase(unittest.TestCase):

  def _validate(self, vals, tol=1e-2, show=False):
    for smi, ans in vals:
      mol = Chem.MolFromSmiles(smi)
      types = AtomTypes.TypeAtoms(mol)
      if show:  # pragma: nocover
        print(types)
      self.assertEqual(len(ans), len(types), 'bad type len for smiles: %s' % (smi))
      lens = [len(x) for x in types]
      self.assertEqual(max(lens), 1, 'atom matched multiple types for smiles: %s' % (smi))
      for a, b in zip(ans, [x[0] for x in types]):
        self.assertEqual(a, b, 'bad type for SMILES: %s' % (smi))

  def test1_simpleMolecules(self):
    data = [
      ('CC', ['sCH3', 'sCH3']),
      ('CCC', ['sCH3', 'ssCH2', 'sCH3']),
      ('CCOC', ['sCH3', 'ssCH2', 'ssO', 'sCH3']),
      ('c1ccccc1[NH3+]', ['aaCH', 'aaCH', 'aaCH', 'aaCH', 'aaCH', 'aasC', 'sNH3']),
      ('c1ccccc1N', ['aaCH', 'aaCH', 'aaCH', 'aaCH', 'aaCH', 'aasC', 'sNH2']),
      ('C#C', ['tCH', 'tCH']),
      ('C=C=C', ['dCH2', 'ddC', 'dCH2']),
    ]
    self._validate(data, show=False)

  def test2_complexMolecules(self):
    data = [
      ('c1[nH]cnc1CC(N)C(O)=O', ['aaCH', 'aaNH', 'aaCH', 'aaN', 'aasC', 'ssCH2', 'sssCH', 'sNH2',
                                 'dssC', 'sOH', 'dO']),
      ('c1nc[n-]c1CC(N)C(O)=O', ['aaCH', 'aaN', 'aaCH', 'aaN', 'aasC', 'ssCH2', 'sssCH', 'sNH2',
                                 'dssC', 'sOH', 'dO']),
      ('c1cccc2c1cccc2',
       ['aaCH', 'aaCH', 'aaCH', 'aaCH', 'aaaC', 'aaaC', 'aaCH', 'aaCH', 'aaCH', 'aaCH']),
    ]
    self._validate(data, show=False)


if __name__ == '__main__':
  unittest.main()
