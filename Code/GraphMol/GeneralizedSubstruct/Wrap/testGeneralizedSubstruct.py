#
#  Copyright (C) 2023 Greg Landrum and other RDKit contributors
#   @@ All Rights Reserved @@
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.

import unittest

#
from rdkit import Chem
from rdkit.Chem import rdGeneralizedSubstruct


class TestCase(unittest.TestCase):

  def test1CreationAndSubstruct(self):
    m = Chem.MolFromSmiles('COCc1n[nH]c(F)c1 |LN:1:1.3|')
    xqm = rdGeneralizedSubstruct.CreateExtendedQueryMol(m)
    matches = rdGeneralizedSubstruct.MolGetXQMSubstructMatches(
      Chem.MolFromSmiles('COOCc1n[nH]c(F)c1'), xqm)
    self.assertEqual(len(matches), 1)
    self.assertEqual(tuple(matches[0]), (0, 1, 2, 3, 4, 5, 9, 6, 7, 8))

    match = rdGeneralizedSubstruct.MolGetXQMSubstructMatch(Chem.MolFromSmiles('COOCc1n[nH]c(F)c1'),
                                                           xqm)
    self.assertEqual(tuple(match), (0, 1, 2, 3, 4, 5, 9, 6, 7, 8))

    #     CHECK(SubstructMatch(*"FCN(C)CC"_smiles, xqm).size() == 1);
    # CHECK(SubstructMatch(*"FCN(O)N(C)CC"_smiles, xqm).size() == 1);


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
