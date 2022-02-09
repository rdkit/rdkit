#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import doctest
import unittest

from rdkit.Chem import inchi
from rdkit.TestRunner import redirect_stderr
import io
if inchi.INCHI_AVAILABLE:
  from rdkit.Chem.MolKey.InchiInfo import InchiInfo

try:
  from rdkit.Chem.MolKey import MolKey
  _testMolKey = True
except ImportError:
  _testMolKey = False


def load_tests(loader, tests, ignore):
  """ Add the Doctests from the module """
  if _testMolKey:
    tests.addTests(doctest.DocTestSuite(MolKey, optionflags=doctest.ELLIPSIS))
  return tests


@unittest.skipUnless(_testMolKey, 'Avalon tools and Inchi required')
class TestMolKey(unittest.TestCase):

  def test_GetKeyForCTAB(self):
    f = io.StringIO()
    with redirect_stderr(f):
      res = MolKey.GetKeyForCTAB('IncorrectCTAB')
      self.assertNotEqual(res.error, 0)
    s = f.getvalue()
    self.assertIn('WARNING:', s)

  def test_CheckCTAB(self):
    self.assertRaises(MolKey.BadMoleculeException, MolKey.CheckCTAB, None)
    self.assertRaises(MolKey.BadMoleculeException, MolKey.CheckCTAB, '')

    ok, _ = MolKey.CheckCTAB('CCincorrect', isSmiles=True)
    self.assertEqual(ok, 1)

    ok, _ = MolKey.CheckCTAB('NO_STRUCTURE', isSmiles=True)
    self.assertEqual(ok, MolKey.ERROR_DICT['NULL_MOL'])

    ok, ctab = MolKey.CheckCTAB('CC', isSmiles=True)
    self.assertEqual(ok, 0)
    ok, ctab2 = MolKey.CheckCTAB(ctab, isSmiles=False)
    self.assertEqual(ok, 0)
    self.assertEqual(ctab, ctab2)

  def test_GetInchiForCTAB(self):
    self.assertRaises(MolKey.BadMoleculeException, MolKey.GetInchiForCTAB, 'IncorrectCTAB')

  def test_ErrorBitsToText(self):
    errors = MolKey.ErrorBitsToText(3)
    self.assertIn('BAD_MOLECULE', errors)
    self.assertIn('ALIAS_CONVERSION_FAILED', errors)
    for k, v in MolKey.ERROR_DICT.items():
      errors = MolKey.ErrorBitsToText(v)
      self.assertEqual(len(errors), 1)
      self.assertIn(k, errors)

  def test_get_chiral_identification_string(self):
    cases = [((0, 0), 'S_ACHIR'),  # No stereo centers
             ((0, 1), 'R_ONE'),  # One undefined stereo centers
             ((0, 2), 'S_UNKN'),  # More than one undefined stereo centers
             ((0, 3), 'S_UNKN'),  # More than one undefined stereo centers
             ((1, 0), 'S_ABS'),  # Fully defined stereo center
             ((2, 0), 'S_ABS'),  # Fully defined stereo centers
             ((1, 1), 'S_PART'),  # Partially defined stereo centers
             ((2, 1), 'S_PART'),  # Partially defined stereo centers
             ]
    for (nDefined, nUndefined), expected in cases:
      self.assertEqual(MolKey._get_chiral_identification_string(nDefined, nUndefined), expected)


GUANINE = 'InChI=1S/C5H5N5O/c6-5-9-3-2(4(11)10-5)7-1-8-3/h1H0,(H4,6,7,8,9,10,11)'
# 'N=C(-O)N', '/FixedH /SUU'
UREA1 = 'InChI=1/CH4N2O/c2-1(3)4/h(H4,2,3,4)/f/h2,4H,3H2/b2-1?'
# 'NC(=O)N', '/FixedH /SUU'
UREA2 = 'InChI=1/CH4N2O/c2-1(3)4/h(H4,2,3,4)/f/h2-3H2'
TRITIATED_UREA = 'InChI=1S/CH4N2O/c2-1(3)4/h(H4,2,3,4)/i/hT3'
DEUTERATED_UREA = 'InChI=1S/CH4N2O/c2-1(3)4/h(H4,2,3,4)/i/hD2'
ACETIC_ACID = 'InChI=1S/C3H6O2/c1-2-3(4)5/h2H2,1H3,(H,4,5)'
ACETATE = 'InChI=1S/C3H6O2/c1-2-3(4)5/h2H2,1H3,(H,4,5)/p-1'
mobile1 = 'InChI=1S/C5H5N3O2/c6-4(9)3-1-7-2-8-5(3)10/h1-2H,(H2,6,9)(H,7,8,10)'  # invented
mobile2 = 'InChI=1S/C7H10N4O/c1-4-2-5(3-6(8)12)11-7(9)10-4/h2H,3H2,1H3,(H2,8,12)(H2,9,10,11)'

# sp3 stereo
sugar1 = 'InChI=1S/C14H20O9/c1-6-11(20-7(2)15)12(21-8(3)16)13(22-9(4)17)14(19-6)23-10(5)18/h6,11-14H,1-5H3/t6-,11-,12+,13+,14?/m0/s1'  # L-rhamnopyranose (source: chemspider)
sugar2 = 'InChI=1S/C12H20O6/c1-11(2)14-5-6(16-11)8-7(13)9-10(15-8)18-12(3,4)17-9/h6-10,13H,5H2,1-4H3/t6-,7-,8-,9-,10-/m1/s1'  # MFCD00135634 (Diacetone-D-Glucose, souce: chemspider)
sp3_unk = 'InChI=1S/C12H21NO4/c1-8(2)10(12(15)16-3)13-11(14)9-5-4-6-17-7-9/h8-10H,4-7H2,1-3H3,(H,13,14)/t9?,10-/m0/s1'  # derived from ChemSpider 34044335


@unittest.skipUnless(inchi.INCHI_AVAILABLE, 'Inchi required')
class TestInchiInfo(unittest.TestCase):

  def doTest(self, inchi, numSp3=0, numUndefSp3=0, numMobileHGroups=0, layer='non-isotopic'):
    ii = InchiInfo(inchi)
    nSp3, nUndefSp3, _, _ = ii.get_sp3_stereo()['main'][layer]
    self.assertEqual(nSp3, numSp3)
    self.assertEqual(nUndefSp3, numUndefSp3)

    nMobileHGroups, _ = ii.get_mobile_h()['main'][layer]
    self.assertEqual(nMobileHGroups, numMobileHGroups)

  def testGuanine(self):
    self.doTest(GUANINE, 0, 0, 1)

  def testTritiatedUrea(self):
    self.doTest(TRITIATED_UREA, 0, 0, 1)

  def testDeuteratedUrea(self):
    self.doTest(DEUTERATED_UREA, 0, 0, 1)

  def testAceticAcid(self):
    self.doTest(ACETIC_ACID, 0, 0, 1)

  def testAcetate(self):
    self.doTest(ACETATE, 0, 0, 1)

  def testMobile1(self):
    self.doTest(mobile1, 0, 0, 2)

  def testMobile2(self):
    self.doTest(mobile2, 0, 0, 2)

  # sp3 stereo
  def testSugar1(self):
    self.doTest(sugar1, 5, 1, 0)

  def testSugar2(self):
    self.doTest(sugar2, 5, 0, 0)

  def testSP3_unk(self):
    self.doTest(sp3_unk, 2, 1, 1)


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
