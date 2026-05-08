#
#  Copyright (C) 2023 Greg Landrum and other RDKit contributors
#   @@ All Rights Reserved @@
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.

import unittest

from rdkit import Chem
from rdkit.Chem import rdGeneralizedSubstruct
from rdkit import rdBase


@unittest.skipIf(not rdBase._serializationEnabled, "not built with serialization support")
class TestCase(unittest.TestCase):

  def test1CreationAndSubstruct(self):
    m = Chem.MolFromSmiles('COCc1n[nH]c(F)c1 |LN:1:1.3|')
    xqm = rdGeneralizedSubstruct.CreateExtendedQueryMol(m)

    mol = Chem.MolFromSmiles('COOCc1n[nH]c(F)c1')
    matches = rdGeneralizedSubstruct.MolGetSubstructMatches(mol, xqm)
    self.assertEqual(len(matches), 1)
    self.assertEqual(tuple(matches[0]), (0, 1, 2, 3, 4, 5, 9, 6, 7, 8))

    match = rdGeneralizedSubstruct.MolGetSubstructMatch(mol, xqm)
    self.assertEqual(tuple(match), (0, 1, 2, 3, 4, 5, 9, 6, 7, 8))

    self.assertTrue(rdGeneralizedSubstruct.MolHasSubstructMatch(mol, xqm))

  def test2Params(self):
    m = Chem.MolFromSmiles('COCc1n[nH]c(F)c1[C@H](F)Cl |LN:1:1.3|')
    xqm = rdGeneralizedSubstruct.CreateExtendedQueryMol(m)

    mol1 = Chem.MolFromSmiles('COOCc1n[nH]c(F)c1[C@H](F)Cl')
    mol2 = Chem.MolFromSmiles('COOCc1n[nH]c(F)c1[C@@H](F)Cl')
    self.assertEqual(len(rdGeneralizedSubstruct.MolGetSubstructMatches(mol1, xqm)), 1)
    self.assertEqual(len(rdGeneralizedSubstruct.MolGetSubstructMatches(mol2, xqm)), 1)

    params = Chem.SubstructMatchParameters()
    params.useChirality = True
    self.assertEqual(len(rdGeneralizedSubstruct.MolGetSubstructMatches(mol1, xqm, params)), 1)
    self.assertEqual(len(rdGeneralizedSubstruct.MolGetSubstructMatches(mol2, xqm, params)), 0)

  def test3MultipleResults(self):
    m = Chem.MolFromSmiles('COC |LN:1:1.3|')
    xqm = rdGeneralizedSubstruct.CreateExtendedQueryMol(m)

    mol1 = Chem.MolFromSmiles('COOCC')
    mol2 = Chem.MolFromSmiles('COCC')

    params = Chem.SubstructMatchParameters()
    params.uniquify = False
    self.assertEqual(len(rdGeneralizedSubstruct.MolGetSubstructMatches(mol1, xqm, params)), 2)
    self.assertEqual(len(rdGeneralizedSubstruct.MolGetSubstructMatches(mol2, xqm, params)), 2)

  def test4AdjustQueryProperties(self):
    m = Chem.MolFromSmiles('COCc1n[nH]cc1 |LN:1:1.3|')
    xqm1 = rdGeneralizedSubstruct.CreateExtendedQueryMol(m)
    xqm2 = rdGeneralizedSubstruct.CreateExtendedQueryMol(m, adjustQueryProperties=True)
    nops = Chem.AdjustQueryParameters.NoAdjustments()
    xqm3 = rdGeneralizedSubstruct.CreateExtendedQueryMol(m, adjustQueryProperties=True,
                                                         adjustQueryParameters=nops)

    mol1 = Chem.MolFromSmiles('COCc1n[nH]cc1')
    mol2 = Chem.MolFromSmiles('COCc1n[nH]c(F)c1')
    self.assertEqual(len(rdGeneralizedSubstruct.MolGetSubstructMatches(mol1, xqm1)), 1)
    self.assertEqual(len(rdGeneralizedSubstruct.MolGetSubstructMatches(mol1, xqm2)), 1)
    self.assertEqual(len(rdGeneralizedSubstruct.MolGetSubstructMatches(mol1, xqm3)), 1)
    self.assertEqual(len(rdGeneralizedSubstruct.MolGetSubstructMatches(mol2, xqm1)), 1)
    self.assertEqual(len(rdGeneralizedSubstruct.MolGetSubstructMatches(mol2, xqm2)), 0)
    self.assertEqual(len(rdGeneralizedSubstruct.MolGetSubstructMatches(mol2, xqm3)), 1)

  def test5ControlSteps(self):
    m = Chem.MolFromSmiles('COCC1OC(N)=N1 |LN:1:1.3|')
    xqm1 = rdGeneralizedSubstruct.CreateExtendedQueryMol(m)
    xqm2 = rdGeneralizedSubstruct.CreateExtendedQueryMol(m, doEnumeration=False)
    xqm3 = rdGeneralizedSubstruct.CreateExtendedQueryMol(m, doTautomers=False)
    xqm4 = rdGeneralizedSubstruct.CreateExtendedQueryMol(m, doEnumeration=False, doTautomers=False)

    mol1 = Chem.MolFromSmiles('COCC1OC(N)=N1')
    self.assertEqual(len(rdGeneralizedSubstruct.MolGetSubstructMatches(mol1, xqm1)), 1)
    self.assertEqual(len(rdGeneralizedSubstruct.MolGetSubstructMatches(mol1, xqm2)), 1)
    self.assertEqual(len(rdGeneralizedSubstruct.MolGetSubstructMatches(mol1, xqm3)), 1)
    self.assertEqual(len(rdGeneralizedSubstruct.MolGetSubstructMatches(mol1, xqm4)), 1)

    mol2 = Chem.MolFromSmiles('COCC1OC(=N)N1')
    self.assertEqual(len(mol2.GetSubstructMatches(m)), 0)
    self.assertEqual(len(rdGeneralizedSubstruct.MolGetSubstructMatches(mol2, xqm1)), 1)
    self.assertEqual(len(rdGeneralizedSubstruct.MolGetSubstructMatches(mol2, xqm2)), 1)
    self.assertEqual(len(rdGeneralizedSubstruct.MolGetSubstructMatches(mol2, xqm3)), 0)
    self.assertEqual(len(rdGeneralizedSubstruct.MolGetSubstructMatches(mol2, xqm4)), 0)

    mol3 = Chem.MolFromSmiles('COOCC1OC(N)=N1')
    self.assertEqual(len(rdGeneralizedSubstruct.MolGetSubstructMatches(mol3, xqm1)), 1)
    self.assertEqual(len(rdGeneralizedSubstruct.MolGetSubstructMatches(mol3, xqm2)), 0)
    self.assertEqual(len(rdGeneralizedSubstruct.MolGetSubstructMatches(mol3, xqm3)), 1)
    self.assertEqual(len(rdGeneralizedSubstruct.MolGetSubstructMatches(mol3, xqm4)), 0)

  def test6Serialization(self):
    m = Chem.MolFromSmiles('COCc1n[nH]c(F)c1 |LN:1:1.3|')
    xqm = rdGeneralizedSubstruct.CreateExtendedQueryMol(m)

    mol = Chem.MolFromSmiles('COOCc1n[nH]c(F)c1')
    matches = rdGeneralizedSubstruct.MolGetSubstructMatches(mol, xqm)
    self.assertEqual(len(matches), 1)
    self.assertEqual(tuple(matches[0]), (0, 1, 2, 3, 4, 5, 9, 6, 7, 8))

    xqm2 = rdGeneralizedSubstruct.ExtendedQueryMol(xqm.ToBinary())
    matches = rdGeneralizedSubstruct.MolGetSubstructMatches(mol, xqm2)
    self.assertEqual(len(matches), 1)

  def test7GenericMatchers(self):
    m = Chem.MolFromSmiles('COC1=NNC(*)=C1 |$;;;;;;AEL_p;$,LN:1:1.3|')
    qps = Chem.AdjustQueryParameters.NoAdjustments()
    qps.makeDummiesQueries = True
    m = Chem.AdjustQueryProperties(m, qps)
    Chem.SetGenericQueriesFromProperties(m)
    ps = Chem.SubstructMatchParameters()
    ps.useGenericMatchers = True

    mol1 = Chem.MolFromSmiles('COC1=NNC(C=C)=C1')
    mol2 = Chem.MolFromSmiles('COC1=NNC(CC)=C1')
    mol3 = Chem.MolFromSmiles('COOC1=NNC(C=C)=C1')
    mol4 = Chem.MolFromSmiles('COOC1=NNC(CC)=C1')

    self.assertTrue(mol1.HasSubstructMatch(m, ps))
    self.assertFalse(mol2.HasSubstructMatch(m, ps))
    self.assertFalse(mol3.HasSubstructMatch(m, ps))
    self.assertFalse(mol4.HasSubstructMatch(m, ps))

    xqm = rdGeneralizedSubstruct.CreateExtendedQueryMol(m)
    self.assertTrue(rdGeneralizedSubstruct.MolHasSubstructMatch(mol1, xqm, ps))
    self.assertFalse(rdGeneralizedSubstruct.MolHasSubstructMatch(mol2, xqm, ps))
    self.assertTrue(rdGeneralizedSubstruct.MolHasSubstructMatch(mol3, xqm, ps))
    self.assertFalse(rdGeneralizedSubstruct.MolHasSubstructMatch(mol4, xqm, ps))

    xqm2 = rdGeneralizedSubstruct.ExtendedQueryMol(xqm.ToBinary())
    self.assertTrue(rdGeneralizedSubstruct.MolHasSubstructMatch(mol1, xqm2, ps))
    self.assertFalse(rdGeneralizedSubstruct.MolHasSubstructMatch(mol2, xqm2, ps))
    self.assertTrue(rdGeneralizedSubstruct.MolHasSubstructMatch(mol3, xqm2, ps))
    self.assertFalse(rdGeneralizedSubstruct.MolHasSubstructMatch(mol4, xqm2, ps))


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
