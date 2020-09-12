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

  def testParsingAbbrevs(self):
    defn = '''CO2Et    C(=O)OCC
COOEt    C(=O)OCC
OiBu     OCC(C)C
tBu      C(C)(C)C'''
    abbrevs = rdAbbreviations.ParseAbbreviations(defn)
    m = Chem.MolFromSmiles('CCC(=O)OCC')
    nm = rdAbbreviations.CondenseMolAbbreviations(m, abbrevs, maxCoverage=1.0)
    self.assertEqual(Chem.MolToCXSmiles(nm), '*CC |$CO2Et;;$|')

  def testCondense(self):
    m = Chem.MolFromSmiles('FC(F)(F)CC(=O)O')
    nm = rdAbbreviations.CondenseMolAbbreviations(m, self.defaultAbbrevs, maxCoverage=1.0)
    self.assertEqual(Chem.MolToCXSmiles(nm), '*C* |$CF3;;CO2H$|')
    m = Chem.MolFromSmiles('CCC(F)(F)F')
    nm = rdAbbreviations.CondenseMolAbbreviations(m, self.defaultAbbrevs)
    self.assertEqual(Chem.MolToCXSmiles(nm), '*C(F)(F)F |$Et;;;;$|')

    # make sure we don't mess up chirality
    m = Chem.MolFromSmiles('FC(F)(F)[C@](Cl)(F)I')
    nm = rdAbbreviations.CondenseMolAbbreviations(m, self.defaultAbbrevs, maxCoverage=1.0)
    self.assertEqual(Chem.MolToCXSmiles(nm), '*[C@@](F)(Cl)I |$CF3;;;;$|')

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

    m = Chem.MolFromSmiles('COC1CCC(C)CC1')
    nm = rdAbbreviations.CondenseMolAbbreviations(m, self.customLinkers, maxCoverage=1.0)
    self.assertEqual(Chem.MolToCXSmiles(nm), 'C*OC |$;Cy;;$|')

  def testAbbreviationsAndLinkers(self):
    m = Chem.MolFromSmiles('COC1CCC(C)CC1')
    # wouldn't normally do this in this order:
    nm = rdAbbreviations.CondenseMolAbbreviations(m, self.defaultAbbrevs, maxCoverage=1.0)
    self.assertEqual(Chem.MolToCXSmiles(nm), '*C1CCC(C)CC1 |$OMe;;;;;;;$|')
    nm = rdAbbreviations.CondenseMolAbbreviations(nm, self.customLinkers, maxCoverage=1.0)
    self.assertEqual(Chem.MolToCXSmiles(nm), '**C |$OMe;Cy;$|')

    # This is a more logical order
    nm = rdAbbreviations.CondenseMolAbbreviations(m, self.customLinkers, maxCoverage=1.0)
    self.assertEqual(Chem.MolToCXSmiles(nm), 'C*OC |$;Cy;;$|')
    nm = rdAbbreviations.CondenseMolAbbreviations(nm, self.defaultAbbrevs, maxCoverage=1.0)
    self.assertEqual(Chem.MolToCXSmiles(nm), 'C*OC |$;Cy;;$|')


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
