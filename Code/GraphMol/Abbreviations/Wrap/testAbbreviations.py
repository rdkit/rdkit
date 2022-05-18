#
#  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
#   @@ All Rights Reserved @@
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.

#
from rdkit import Chem
from rdkit.Chem import rdAbbreviations

import unittest


class TestCase(unittest.TestCase):

  def setUp(self):
    self.defaultAbbrevs = rdAbbreviations.GetDefaultAbbreviations()
    self.defaultLinkers = rdAbbreviations.GetDefaultLinkers()
    self.customLinkers = rdAbbreviations.ParseLinkers('''PEG3  *OCCOCCOCC* PEG3
Pent  *CCCCC*
Cy   *C1CCC(*)CC1  Cy''')
    self.parseMapping = lambda mapping: tuple(map(int, mapping.strip('[],').split(',')))

  def testParsingAbbrevs(self):
    defn = '''CO2Et    C(=O)OCC
COOEt    C(=O)OCC
OiBu     OCC(C)C
tBu      C(C)(C)C'''
    abbrevs = rdAbbreviations.ParseAbbreviations(defn)
    m = Chem.MolFromSmiles('CCC(=O)OCC')
    nm = rdAbbreviations.CondenseMolAbbreviations(m, abbrevs, maxCoverage=1.0)
    self.assertEqual(Chem.MolToCXSmiles(nm), '*CC |$CO2Et;;$|')
    self.assertTrue(nm.HasProp('_origAtomMapping'))
    self.assertEqual(self.parseMapping(nm.GetProp('_origAtomMapping')), (0, 1, 2))
    self.assertTrue(nm.HasProp('_origBondMapping'))
    self.assertEqual(self.parseMapping(nm.GetProp('_origBondMapping')), (0, 1))

  def testCondense(self):
    m = Chem.MolFromSmiles('FC(F)(F)CC(=O)O')
    nm = rdAbbreviations.CondenseMolAbbreviations(m, self.defaultAbbrevs, maxCoverage=1.0)
    self.assertEqual(Chem.MolToCXSmiles(nm), '*C* |$CF3;;CO2H$|')
    self.assertTrue(nm.HasProp('_origAtomMapping'))
    self.assertEqual(self.parseMapping(nm.GetProp('_origAtomMapping')), (1, 4, 5))
    self.assertTrue(nm.HasProp('_origBondMapping'))
    self.assertEqual(self.parseMapping(nm.GetProp('_origBondMapping')), (3, 4))

    m = Chem.MolFromSmiles('CCC(F)(F)F')
    nm = rdAbbreviations.CondenseMolAbbreviations(m, self.defaultAbbrevs)
    self.assertEqual(Chem.MolToCXSmiles(nm), '*C(F)(F)F |$Et;;;;$|')
    self.assertTrue(nm.HasProp('_origAtomMapping'))
    self.assertEqual(self.parseMapping(nm.GetProp('_origAtomMapping')), (1, 2, 3, 4, 5))
    self.assertTrue(nm.HasProp('_origBondMapping'))
    self.assertEqual(self.parseMapping(nm.GetProp('_origBondMapping')), (1, 2, 3, 4))

    # make sure we don't mess up chirality
    m = Chem.MolFromSmiles('FC(F)(F)[C@](Cl)(F)I')
    nm = rdAbbreviations.CondenseMolAbbreviations(m, self.defaultAbbrevs, maxCoverage=1.0)
    self.assertEqual(Chem.MolToCXSmiles(nm), '*[C@@](F)(Cl)I |$CF3;;;;$|')
    self.assertTrue(nm.HasProp('_origAtomMapping'))
    self.assertEqual(self.parseMapping(nm.GetProp('_origAtomMapping')), (1, 4, 5, 6, 7))
    self.assertTrue(nm.HasProp('_origBondMapping'))
    self.assertEqual(self.parseMapping(nm.GetProp('_origBondMapping')), (3, 4, 5, 6))

  def testLabel(self):
    m = Chem.MolFromSmiles('CC(C)CC(F)(F)F')
    nm = rdAbbreviations.LabelMolAbbreviations(m, self.defaultAbbrevs, maxCoverage=1.0)
    sgs = Chem.GetMolSubstanceGroups(nm)
    self.assertEqual(len(sgs), 2)
    self.assertEqual(sgs[0].GetProp('TYPE'), "SUP")
    self.assertEqual(sgs[0].GetProp('LABEL'), "iPr")
    self.assertEqual(list(sgs[0].GetAtoms()), [1, 0, 2])
    self.assertEqual(list(sgs[0].GetBonds()), [2])
    aps = sgs[0].GetAttachPoints()
    self.assertEqual(len(aps), 1)
    self.assertEqual(aps[0].aIdx, 1)
    self.assertEqual(aps[0].lvIdx, 3)

    self.assertEqual(sgs[1].GetProp('TYPE'), "SUP")
    self.assertEqual(sgs[1].GetProp('LABEL'), "CF3")
    self.assertEqual(list(sgs[1].GetAtoms()), [4, 5, 6, 7])
    self.assertEqual(list(sgs[1].GetBonds()), [3])
    aps = sgs[1].GetAttachPoints()
    self.assertEqual(len(aps), 1)
    self.assertEqual(aps[0].aIdx, 4)
    self.assertEqual(aps[0].lvIdx, 3)

  def testCondenseLinkers(self):
    m = Chem.MolFromSmiles('FCOCCOCCOCCCCCCCCCCl')
    nm = rdAbbreviations.CondenseMolAbbreviations(m, self.defaultLinkers, maxCoverage=1.0)
    self.assertEqual(Chem.MolToCXSmiles(nm), 'FC**Cl |$;;PEG3;Hept;$|')
    self.assertTrue(nm.HasProp('_origAtomMapping'))
    self.assertEqual(self.parseMapping(nm.GetProp('_origAtomMapping')), (0, 1, 2, 11, 18))
    self.assertTrue(nm.HasProp('_origBondMapping'))
    self.assertEqual(self.parseMapping(nm.GetProp('_origBondMapping')), (0, 1, 10, 17))

    m = Chem.MolFromSmiles('COC1CCC(C)CC1')
    nm = rdAbbreviations.CondenseMolAbbreviations(m, self.customLinkers, maxCoverage=1.0)
    self.assertEqual(Chem.MolToCXSmiles(nm), 'C*OC |$;Cy;;$|')
    self.assertTrue(nm.HasProp('_origAtomMapping'))
    self.assertEqual(self.parseMapping(nm.GetProp('_origAtomMapping')), (0, 1, 2, 6))
    self.assertTrue(nm.HasProp('_origBondMapping'))
    self.assertEqual(self.parseMapping(nm.GetProp('_origBondMapping')), (0, 1, 5))

  def testAbbreviationsAndLinkers(self):
    m = Chem.MolFromSmiles('COC1CCC(C)CC1')
    # wouldn't normally do this in this order:
    nm = rdAbbreviations.CondenseMolAbbreviations(m, self.defaultAbbrevs, maxCoverage=1.0)
    self.assertEqual(Chem.MolToCXSmiles(nm), '*C1CCC(C)CC1 |$OMe;;;;;;;$|')
    self.assertTrue(nm.HasProp('_origAtomMapping'))
    self.assertEqual(self.parseMapping(nm.GetProp('_origAtomMapping')), (1, 2, 3, 4, 5, 6, 7, 8))
    self.assertTrue(nm.HasProp('_origBondMapping'))
    self.assertEqual(self.parseMapping(nm.GetProp('_origBondMapping')), (1, 2, 3, 4, 5, 6, 7, 8))
    nm = rdAbbreviations.CondenseMolAbbreviations(nm, self.customLinkers, maxCoverage=1.0)
    self.assertEqual(Chem.MolToCXSmiles(nm), '**C |$OMe;Cy;$|')
    self.assertTrue(nm.HasProp('_origAtomMapping'))
    self.assertEqual(self.parseMapping(nm.GetProp('_origAtomMapping')), (1, 2, 6))
    self.assertTrue(nm.HasProp('_origBondMapping'))
    self.assertEqual(self.parseMapping(nm.GetProp('_origBondMapping')), (1, 5))

    # This is a more logical order
    nm = rdAbbreviations.CondenseMolAbbreviations(m, self.customLinkers, maxCoverage=1.0)
    self.assertEqual(Chem.MolToCXSmiles(nm), 'C*OC |$;Cy;;$|')
    self.assertTrue(nm.HasProp('_origAtomMapping'))
    self.assertEqual(self.parseMapping(nm.GetProp('_origAtomMapping')), (0, 1, 2, 6))
    self.assertTrue(nm.HasProp('_origBondMapping'))
    self.assertEqual(self.parseMapping(nm.GetProp('_origBondMapping')), (0, 1, 5))
    nm = rdAbbreviations.CondenseMolAbbreviations(nm, self.defaultAbbrevs, maxCoverage=1.0)
    self.assertEqual(Chem.MolToCXSmiles(nm), 'C*OC |$;Cy;;$|')
    self.assertTrue(nm.HasProp('_origAtomMapping'))
    self.assertEqual(self.parseMapping(nm.GetProp('_origAtomMapping')), (0, 1, 2, 6))
    self.assertTrue(nm.HasProp('_origBondMapping'))
    self.assertEqual(self.parseMapping(nm.GetProp('_origBondMapping')), (0, 1, 5))

  def testAbbreviationsSubstanceGroups(self):
    m = Chem.MolFromMolBlock('''
  Mrv2014 09152006492D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 7 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 5.25 -5.9858 0 0
M  V30 2 C 4.48 -7.3196 0 0
M  V30 3 C 6.02 -7.3196 0 0
M  V30 4 F 8.6873 -8.8596 0 0
M  V30 5 C 7.3537 -8.0896 0 0
M  V30 6 F 6.02 -8.8596 0 0
M  V30 7 F 7.3537 -6.5496 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 3 1
M  V30 3 1 2 3
M  V30 4 1 3 5
M  V30 5 1 4 5
M  V30 6 1 5 6
M  V30 7 1 5 7
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SUP 0 ATOMS=(4 4 5 6 7) SAP=(3 5 3 1) XBONDS=(1 4) LABEL=CF3
M  V30 END SGROUP
M  V30 END CTAB
M  END''')
    nm = rdAbbreviations.CondenseAbbreviationSubstanceGroups(m)
    nm.RemoveAllConformers()  # avoid coords in CXSMILES
    self.assertEqual(Chem.MolToCXSmiles(nm), '*C1CC1 |$CF3;;;$|')
    self.assertTrue(nm.HasProp('_origAtomMapping'))
    self.assertEqual(self.parseMapping(nm.GetProp('_origAtomMapping')), (0, 1, 2, 4))
    self.assertTrue(nm.HasProp('_origBondMapping'))
    self.assertEqual(self.parseMapping(nm.GetProp('_origBondMapping')), (0, 1, 2, 3))

  def testGithub3692(self):
    defaults = rdAbbreviations.GetDefaultAbbreviations()
    self.assertIsNotNone(defaults[0].mol)
    lbls = [x.label for x in defaults]
    self.assertIn('CO2Et', lbls)
    idx = lbls.index('CO2Et')
    self.assertEqual(Chem.MolToSmiles(defaults[idx].mol), '*C(=O)OCC')


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
