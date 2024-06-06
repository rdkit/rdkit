#
#  Copyright (c) 2007-2023, Novartis Institutes for BioMedical Research Inc. and
#  other RDKit contributors
#
#  All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

import doctest
import importlib.util
import os
import pickle
import sys
import unittest
import tempfile

from rdkit import Chem, Geometry, RDConfig, rdBase
from rdkit.Chem import AllChem, rdChemReactions
from rdkit.Chem.SimpleEnum import Enumerator


def feq(v1, v2, tol2=1e-4):
  return abs(v1 - v2) <= tol2


def ptEq(pt1, pt2, tol=1e-4):
  return feq(pt1.x, pt2.x, tol) and feq(pt1.y, pt2.y, tol) and feq(pt1.z, pt2.z, tol)


# Boost functions are NOT found by doctest, this "fixes" them
#  by adding the doctests to a fake module
spec = importlib.util.spec_from_loader("TestPreprocess", loader=None)
TestPreprocess = importlib.util.module_from_spec(spec)
code = """
from rdkit.Chem import rdChemReactions
def PreprocessReaction(*a, **kw):
    '''%s
    '''
    return rdChemReactions.PreprocessReaction(*a, **kw)
""" % "\n".join([x.lstrip() for x in rdChemReactions.PreprocessReaction.__doc__.split("\n")])
exec(code, TestPreprocess.__dict__)


def load_tests(loader, tests, ignore):
  tests.addTests(doctest.DocTestSuite(Enumerator))
  tests.addTests(doctest.DocTestSuite(TestPreprocess))
  return tests


class TestCase(unittest.TestCase):

  def setUp(self):
    self.dataDir = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'ChemReactions', 'testData')

  def test1Basics(self):
    rxna = rdChemReactions.ChemicalReaction()
    # also tests empty copy constructor
    for rxn in [rxna, rdChemReactions.ChemicalReaction(rxna)]:
      self.assertTrue(rxn.GetNumReactantTemplates() == 0)
      self.assertTrue(rxn.GetNumProductTemplates() == 0)

      r1 = Chem.MolFromSmarts('[C:1](=[O:2])O')
      rxn.AddReactantTemplate(r1)
      self.assertTrue(rxn.GetNumReactantTemplates() == 1)

      r1 = Chem.MolFromSmarts('[N:3]')
      rxn.AddReactantTemplate(r1)
      self.assertTrue(rxn.GetNumReactantTemplates() == 2)

      r1 = Chem.MolFromSmarts('[C:1](=[O:2])[N:3]')
      rxn.AddProductTemplate(r1)
      self.assertTrue(rxn.GetNumProductTemplates() == 1)

      reacts = (Chem.MolFromSmiles('C(=O)O'), Chem.MolFromSmiles('N'))
      ps = rxn.RunReactants(reacts)
      self.assertTrue(len(ps) == 1)
      self.assertTrue(len(ps[0]) == 1)
      self.assertTrue(ps[0][0].GetNumAtoms() == 3)

      ps = rxn.RunReactants(list(reacts))
      self.assertTrue(len(ps) == 1)
      self.assertTrue(len(ps[0]) == 1)
      self.assertTrue(ps[0][0].GetNumAtoms() == 3)

  def test2DaylightParser(self):
    rxna = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]')
    for rxn in [rxna, rdChemReactions.ChemicalReaction(rxna)]:
      self.assertTrue(rxn)
      self.assertTrue(rxn.GetNumReactantTemplates() == 2)
      self.assertTrue(rxn.GetNumProductTemplates() == 1)
      self.assertTrue(rxn._getImplicitPropertiesFlag())

      reacts = (Chem.MolFromSmiles('C(=O)O'), Chem.MolFromSmiles('N'))
      ps = rxn.RunReactants(reacts)
      self.assertTrue(len(ps) == 1)
      self.assertTrue(len(ps[0]) == 1)
      self.assertTrue(ps[0][0].GetNumAtoms() == 3)

      reacts = (Chem.MolFromSmiles('CC(=O)OC'), Chem.MolFromSmiles('CN'))
      ps = rxn.RunReactants(reacts)
      self.assertTrue(len(ps) == 1)
      self.assertTrue(len(ps[0]) == 1)
      self.assertTrue(ps[0][0].GetNumAtoms() == 5)

  def test3MDLParsers(self):
    fileN = os.path.join(self.dataDir, 'AmideBond.rxn')
    rxna = rdChemReactions.ReactionFromRxnFile(fileN)
    print("*" * 44)
    print(fileN)
    print(rxna)
    for rxn in [rxna, rdChemReactions.ChemicalReaction(rxna)]:
      self.assertTrue(rxn)
      self.assertFalse(rxn._getImplicitPropertiesFlag())

      self.assertTrue(rxn.GetNumReactantTemplates() == 2)
      self.assertTrue(rxn.GetNumProductTemplates() == 1)

      reacts = (Chem.MolFromSmiles('C(=O)O'), Chem.MolFromSmiles('N'))
      ps = rxn.RunReactants(reacts)
      self.assertTrue(len(ps) == 1)
      self.assertTrue(len(ps[0]) == 1)
      self.assertTrue(ps[0][0].GetNumAtoms() == 3)

      with open(fileN, 'r') as rxnF:
        rxnBlock = rxnF.read()
      rxn = rdChemReactions.ReactionFromRxnBlock(rxnBlock)
      self.assertTrue(rxn)

      self.assertTrue(rxn.GetNumReactantTemplates() == 2)
      self.assertTrue(rxn.GetNumProductTemplates() == 1)

      reacts = (Chem.MolFromSmiles('C(=O)O'), Chem.MolFromSmiles('N'))
      ps = rxn.RunReactants(reacts)
      self.assertTrue(len(ps) == 1)
      self.assertTrue(len(ps[0]) == 1)
      self.assertTrue(ps[0][0].GetNumAtoms() == 3)

  def test4ErrorHandling(self):
    self.assertRaises(
      ValueError,
      lambda x='[C:1](=[O:2])Q.[N:3]>>[C:1](=[O:2])[N:3]': rdChemReactions.ReactionFromSmarts(x))
    self.assertRaises(
      ValueError,
      lambda x='[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]Q': rdChemReactions.ReactionFromSmarts(x))
    self.assertRaises(
      ValueError, lambda x='[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]>>CC': rdChemReactions.
      ReactionFromSmarts(x))

    block = """$RXN

      ISIS     082120061354

  3  1
$MOL

  -ISIS-  08210613542D

  3  2  0  0  0  0  0  0  0  0999 V2000
   -1.4340   -0.6042    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0
   -0.8639   -0.9333    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4340    0.0542    0.0000 O   0  0  0  0  0  0  0  0  0  1  0  0
  1  2  1  0  0  0  0
  1  3  2  0  0  0  0
M  END
$MOL

  -ISIS-  08210613542D

  1  0  0  0  0  0  0  0  0  0999 V2000
    2.2125   -0.7833    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0
M  END
$MOL

  -ISIS-  08210613542D

  3  2  0  0  0  0  0  0  0  0999 V2000
    9.5282   -0.8083    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0
    8.9579   -0.4792    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0
    8.9579    0.1792    0.0000 O   0  0  0  0  0  0  0  0  0  1  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
M  END
    """
    self.assertRaises(ValueError, lambda x=block: rdChemReactions.ReactionFromRxnBlock(x))

    block = """$RXN

      ISIS     082120061354

  2  1
$MOL

  -ISIS-  08210613542D

  4  2  0  0  0  0  0  0  0  0999 V2000
   -1.4340   -0.6042    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0
   -0.8639   -0.9333    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4340    0.0542    0.0000 O   0  0  0  0  0  0  0  0  0  1  0  0
  1  2  1  0  0  0  0
  1  3  2  0  0  0  0
M  END
$MOL

  -ISIS-  08210613542D

  1  0  0  0  0  0  0  0  0  0999 V2000
    2.2125   -0.7833    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0
M  END
$MOL

  -ISIS-  08210613542D

  3  2  0  0  0  0  0  0  0  0999 V2000
    9.5282   -0.8083    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0
    8.9579   -0.4792    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0
    8.9579    0.1792    0.0000 O   0  0  0  0  0  0  0  0  0  1  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
M  END
    """
    #self.assertRaises(ValueError,lambda x=block:rdChemReactions.ReactionFromRxnBlock(x))

    block = """$RXN

      ISIS     082120061354

  2  1
$MOL

  -ISIS-  08210613542D

  3  2  0  0  0  0  0  0  0  0999 V2000
   -1.4340   -0.6042    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0
   -0.8639   -0.9333    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4340    0.0542    0.0000 O   0  0  0  0  0  0  0  0  0  1  0  0
  1  2  1  0  0  0  0
  1  3  2  0  0  0  0
M  END
$MOL

  -ISIS-  08210613542D

  1  0  0  0  0  0  0  0  0  0999 V2000
    2.2125   -0.7833    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0
M  END
$MOL

  -ISIS-  08210613542D

  3  1  0  0  0  0  0  0  0  0999 V2000
    9.5282   -0.8083    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0
    8.9579   -0.4792    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0
    8.9579    0.1792    0.0000 O   0  0  0  0  0  0  0  0  0  1  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
M  END
    """
    #self.assertRaises(ValueError,lambda x=block:rdChemReactions.ReactionFromRxnBlock(x))

  def test5Validation(self):
    rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]')
    self.assertTrue(rxn)
    self.assertTrue(rxn.Validate() == (0, 0))

    rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:1])O.[N:3]>>[C:1](=[O:2])[N:3]')
    self.assertTrue(rxn)
    self.assertTrue(rxn.Validate() == (1, 1))

    rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])[O:4].[N:3]>>[C:1](=[O:2])[N:3]')
    self.assertTrue(rxn)
    self.assertTrue(rxn.Validate() == (1, 0))

    rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3][C:5]')
    self.assertTrue(rxn)
    self.assertTrue(rxn.Validate() == (1, 0))

  def test6Exceptions(self):
    rxn = rdChemReactions.ReactionFromSmarts('[C:1]Cl>>[C:1]')
    self.assertTrue(rxn)
    self.assertRaises(ValueError, lambda x=rxn: x.RunReactants(()))
    self.assertRaises(
      ValueError, lambda x=rxn: x.RunReactants((Chem.MolFromSmiles('CC'), Chem.MolFromSmiles('C'))))
    ps = rxn.RunReactants((Chem.MolFromSmiles('CCCl'), ))
    self.assertTrue(len(ps) == 1)
    self.assertTrue(len(ps[0]) == 1)

  def _test7Leak(self):
    rxn = rdChemReactions.ReactionFromSmarts('[C:1]Cl>>[C:1]')
    self.assertTrue(rxn)
    print('running: ')
    for i in range(1e5):
      ps = rxn.RunReactants((Chem.MolFromSmiles('CCCl'), ))
      self.assertTrue(len(ps) == 1)
      self.assertTrue(len(ps[0]) == 1)
      if not i % 1000:
        print(i)

  def test8Properties(self):
    rxn = rdChemReactions.ReactionFromSmarts('[O:1]>>[O:1][3#0]')
    self.assertTrue(rxn)
    ps = rxn.RunReactants((Chem.MolFromSmiles('CO'), ))
    self.assertTrue(len(ps) == 1)
    self.assertTrue(len(ps[0]) == 1)
    Chem.SanitizeMol(ps[0][0])
    self.assertEqual(ps[0][0].GetAtomWithIdx(1).GetIsotope(), 3)

  def test9AromaticityTransfer(self):
    # this was issue 2664121
    mol = Chem.MolFromSmiles('c1ccc(C2C3(Cc4c(cccc4)C2)CCCC3)cc1')
    rxn = rdChemReactions.ReactionFromSmarts(
      '[A:1]1~[*:2]~[*:3]~[*:4]~[*:5]~[A:6]-;@1>>[*:1]~[*:2]~[*:3]~[*:4]~[*:5]~[*:6]')
    products = rxn.RunReactants([mol])
    self.assertEqual(len(products), 6)
    for p in products:
      self.assertEqual(len(p), 1)
      Chem.SanitizeMol(p[0])

  def test10DotSeparation(self):
    # 08/05/14
    # This test is changed due to a new behavior of the smarts
    # reaction parser which now allows using parenthesis in products
    # as well. original smiles: '[C:1]1[O:2][N:3]1>>[C:1]1[O:2].[N:3]1'
    rxn = rdChemReactions.ReactionFromSmarts('[C:1]1[O:2][N:3]1>>([C:1]1[O:2].[N:3]1)')
    mol = Chem.MolFromSmiles('C1ON1')
    products = rxn.RunReactants([mol])
    self.assertEqual(len(products), 1)
    for p in products:
      self.assertEqual(len(p), 1)
      self.assertEqual(p[0].GetNumAtoms(), 3)
      self.assertEqual(p[0].GetNumBonds(), 2)

  def test11ImplicitProperties(self):
    rxn = rdChemReactions.ReactionFromSmarts('[C:1]O>>[C:1]')
    mol = Chem.MolFromSmiles('CCO')
    products = rxn.RunReactants([mol])
    self.assertEqual(len(products), 1)
    for p in products:
      self.assertEqual(len(p), 1)
      self.assertEqual(Chem.MolToSmiles(p[0]), 'CC')
    mol2 = Chem.MolFromSmiles('C[CH-]O')
    products = rxn.RunReactants([mol2])
    self.assertEqual(len(products), 1)
    for p in products:
      self.assertEqual(len(p), 1)
      self.assertEqual(Chem.MolToSmiles(p[0]), '[CH2-]C')

    rxn._setImplicitPropertiesFlag(False)
    products = rxn.RunReactants([mol])
    self.assertEqual(len(products), 1)
    for p in products:
      self.assertEqual(len(p), 1)
      self.assertEqual(Chem.MolToSmiles(p[0]), 'CC')
    products = rxn.RunReactants([mol2])
    self.assertEqual(len(products), 1)
    for p in products:
      self.assertEqual(len(p), 1)
      self.assertEqual(Chem.MolToSmiles(p[0]), 'CC')

  def test12Pickles(self):
    # 08/05/14
    # This test is changed due to a new behavior of the smarts
    # reaction parser which now allows using parenthesis in products
    # as well. original smiles: '[C:1]1[O:2][N:3]1>>[C:1]1[O:2].[N:3]1'
    rxn = rdChemReactions.ReactionFromSmarts('[C:1]1[O:2][N:3]1>>([C:1]1[O:2].[N:3]1)')
    pkl = pickle.dumps(rxn)
    rxn = pickle.loads(pkl)
    mol = Chem.MolFromSmiles('C1ON1')
    products = rxn.RunReactants([mol])
    self.assertEqual(len(products), 1)
    for p in products:
      self.assertEqual(len(p), 1)
      self.assertEqual(p[0].GetNumAtoms(), 3)
      self.assertEqual(p[0].GetNumBonds(), 2)

    rxn = rdChemReactions.ChemicalReaction(rxn.ToBinary())
    products = rxn.RunReactants([mol])
    self.assertEqual(len(products), 1)
    for p in products:
      self.assertEqual(len(p), 1)
      self.assertEqual(p[0].GetNumAtoms(), 3)
      self.assertEqual(p[0].GetNumBonds(), 2)

  def test13GetTemplates(self):
    rxn = rdChemReactions.ReactionFromSmarts('[C:1]1[O:2][N:3]1>>[C:1][O:2].[N:3]')
    r1 = rxn.GetReactantTemplate(0)
    sma = Chem.MolToSmarts(r1)
    self.assertEqual(sma, '[C:1]1[O:2][N:3]1')
    p1 = rxn.GetProductTemplate(0)
    sma = Chem.MolToSmarts(p1)
    self.assertEqual(sma, '[C:1][O:2]')

    p2 = rxn.GetProductTemplate(1)
    sma = Chem.MolToSmarts(p2)
    self.assertEqual(sma, '[N:3]')

    self.assertRaises(ValueError, lambda: rxn.GetProductTemplate(2))
    self.assertRaises(ValueError, lambda: rxn.GetReactantTemplate(1))

  def test14Matchers(self):
    rxn = rdChemReactions.ReactionFromSmarts(
      '[C;!$(C(-O)-O):1](=[O:2])[O;H,-1].[N;!H0:3]>>[C:1](=[O:2])[N:3]')
    self.assertTrue(rxn)
    rxn.Initialize()
    self.assertTrue(rxn.IsMoleculeReactant(Chem.MolFromSmiles('OC(=O)C')))
    self.assertFalse(rxn.IsMoleculeReactant(Chem.MolFromSmiles('OC(=O)O')))
    self.assertTrue(rxn.IsMoleculeReactant(Chem.MolFromSmiles('CNC')))
    self.assertFalse(rxn.IsMoleculeReactant(Chem.MolFromSmiles('CN(C)C')))
    self.assertTrue(rxn.IsMoleculeProduct(Chem.MolFromSmiles('NC(=O)C')))
    self.assertTrue(rxn.IsMoleculeProduct(Chem.MolFromSmiles('CNC(=O)C')))
    self.assertFalse(rxn.IsMoleculeProduct(Chem.MolFromSmiles('COC(=O)C')))

  def test15Replacements(self):
    rxn = rdChemReactions.ReactionFromSmarts(
      '[{amine}:1]>>[*:1]-C',
      replacements={'{amine}': '$([N;!H0;$(N-[#6]);!$(N-[!#6;!#1]);!$(N-C=[O,N,S])])'})
    self.assertTrue(rxn)
    rxn.Initialize()
    reactants = (Chem.MolFromSmiles('CCN'), )
    ps = rxn.RunReactants(reactants)
    self.assertEqual(len(ps), 1)
    self.assertEqual(len(ps[0]), 1)
    self.assertEqual(ps[0][0].GetNumAtoms(), 4)

  def test16GetReactingAtoms(self):
    rxn = rdChemReactions.ReactionFromSmarts("[O:1][C:2].[N:3]>>[N:1][C:2].[N:3]")
    self.assertTrue(rxn)
    rxn.Initialize()
    rAs = rxn.GetReactingAtoms()
    self.assertEqual(len(rAs), 2)
    self.assertEqual(len(rAs[0]), 1)
    self.assertEqual(len(rAs[1]), 0)

    rxn = rdChemReactions.ReactionFromSmarts("[O:1]C>>[O:1]C")
    self.assertTrue(rxn)
    rxn.Initialize()
    rAs = rxn.GetReactingAtoms()
    self.assertEqual(len(rAs), 1)
    self.assertEqual(len(rAs[0]), 2)
    rAs = rxn.GetReactingAtoms(True)
    self.assertEqual(len(rAs), 1)
    self.assertEqual(len(rAs[0]), 1)

  def test17AddRecursiveQueriesToReaction(self):
    rxn = rdChemReactions.ReactionFromSmarts("[C:1][O:2].[N:3]>>[C:1][N:2]")
    self.assertTrue(rxn)
    rxn.Initialize()
    qs = {'aliphatic': Chem.MolFromSmiles('CC')}
    rxn.GetReactantTemplate(0).GetAtomWithIdx(0).SetProp('query', 'aliphatic')
    rxn.AddRecursiveQueriesToReaction(qs, 'query')
    q = rxn.GetReactantTemplate(0)
    m = Chem.MolFromSmiles('CCOC')
    self.assertTrue(m.HasSubstructMatch(q))
    m = Chem.MolFromSmiles('CO')
    self.assertFalse(m.HasSubstructMatch(q))

    rxn = rdChemReactions.ReactionFromSmarts("[C:1][O:2].[N:3]>>[C:1][N:2]")
    rxn.Initialize()
    rxn.GetReactantTemplate(0).GetAtomWithIdx(0).SetProp('query', 'aliphatic')
    labels = rxn.AddRecursiveQueriesToReaction(qs, 'query', getLabels=True)
    self.assertTrue(len(labels), 1)

  def test17bAddRecursiveQueriesToReaction(self):
    from rdkit.Chem import FilterCatalog
    rxn = rdChemReactions.ReactionFromSmarts("[C:1][O:2].[N:3]>>[C:1][N:2]")
    self.assertTrue(rxn)
    rxn.Initialize()
    rxn.GetReactantTemplate(0).GetAtomWithIdx(0).SetProp('query', 'carboxylicacid')
    querydefs = {
      k.lower(): v
      for k, v in FilterCatalog.GetFlattenedFunctionalGroupHierarchy().items()
    }

    self.assertTrue('CarboxylicAcid' in FilterCatalog.GetFlattenedFunctionalGroupHierarchy())
    rxn.AddRecursiveQueriesToReaction(querydefs, 'query')
    q = rxn.GetReactantTemplate(0)
    m = Chem.MolFromSmiles('C(=O)[O-].N')
    self.assertTrue(m.HasSubstructMatch(q))
    m = Chem.MolFromSmiles('C.N')
    self.assertFalse(m.HasSubstructMatch(q))

  def test18GithubIssue16(self):
    rxn = rdChemReactions.ReactionFromSmarts("[F:1]>>[Cl:1]")
    self.assertTrue(rxn)
    rxn.Initialize()
    self.assertRaises(ValueError, lambda: rxn.RunReactants((None, )))

  def test19RemoveUnmappedMoleculesToAgents(self):
    rxn = rdChemReactions.ReactionFromSmarts(
      "[C:1]=[O:2].[N:3].C(=O)O>[OH2].[Na].[Cl]>[N:3]~[C:1]=[O:2]")
    self.assertTrue(rxn)
    rxn.Initialize()
    self.assertTrue(rxn.GetNumReactantTemplates() == 3)
    self.assertTrue(rxn.GetNumProductTemplates() == 1)
    self.assertTrue(rxn.GetNumAgentTemplates() == 3)

    rxn.RemoveUnmappedReactantTemplates()
    rxn.RemoveUnmappedProductTemplates()

    self.assertTrue(rxn.GetNumReactantTemplates() == 2)
    self.assertTrue(rxn.GetNumProductTemplates() == 1)
    self.assertTrue(rxn.GetNumAgentTemplates() == 4)

    rxn = rdChemReactions.ReactionFromSmarts("[C:1]=[O:2].[N:3].C(=O)O>>[N:3]~[C:1]=[O:2].[OH2]")
    self.assertTrue(rxn)
    rxn.Initialize()
    self.assertTrue(rxn.GetNumReactantTemplates() == 3)
    self.assertTrue(rxn.GetNumProductTemplates() == 2)
    self.assertTrue(rxn.GetNumAgentTemplates() == 0)

    agentList = []
    rxn.RemoveUnmappedReactantTemplates(moveToAgentTemplates=False, targetList=agentList)
    rxn.RemoveUnmappedProductTemplates(targetList=agentList)

    self.assertTrue(rxn.GetNumReactantTemplates() == 2)
    self.assertTrue(rxn.GetNumProductTemplates() == 1)
    self.assertTrue(rxn.GetNumAgentTemplates() == 1)
    self.assertTrue(len(agentList) == 2)

  def test20CheckCopyConstructedReactionAtomProps(self):
    RLABEL = "_MolFileRLabel"
    amine_rxn = '$RXN\n\n      ISIS     090220091541\n\n  2  1\n$MOL\n\n  -ISIS-  09020915412D\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n   -2.9083   -0.4708    0.0000 R#  0  0  0  0  0  0  0  0  0  1  0  0\n   -2.3995   -0.1771    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0\n   -2.4042    0.4125    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  2  3  2  0  0  0  0\nV    2 aldehyde\nM  RGP  1   1   1\nM  END\n$MOL\n\n  -ISIS-  09020915412D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    2.8375   -0.2500    0.0000 R#  0  0  0  0  0  0  0  0  0  3  0  0\n    3.3463    0.0438    0.0000 N   0  0  0  0  0  0  0  0  0  4  0  0\n  1  2  1  0  0  0  0\nV    2 amine\nM  RGP  1   1   2\nM  END\n$MOL\n\n  -ISIS-  09020915412D\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n   13.3088    0.9436    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0\n   13.8206    1.2321    0.0000 R#  0  0  0  0  0  0  0  0  0  1  0  0\n   13.3028    0.3561    0.0000 N   0  0  0  0  0  0  0  0  0  4  0  0\n   12.7911    0.0676    0.0000 R#  0  0  0  0  0  0  0  0  0  3  0  0\n  1  3  1  0  0  0  0\n  1  2  1  0  0  0  0\n  3  4  1  0  0  0  0\nM  RGP  2   2   1   4   2\nM  END\n'
    rxn = rdChemReactions.ReactionFromRxnBlock(amine_rxn)
    res = []
    for atom in rxn.GetReactantTemplate(0).GetAtoms():
      if atom.HasProp(RLABEL):
        res.append((atom.GetIdx(), atom.GetProp(RLABEL)))
    rxn2 = rdChemReactions.ChemicalReaction(rxn)
    res2 = []

    for atom in rxn2.GetReactantTemplate(0).GetAtoms():
      if atom.HasProp(RLABEL):
        res2.append((atom.GetIdx(), atom.GetProp(RLABEL)))
    self.assertEqual(res, res2)

    # currently ToBinary does not save atom props
    # rxn2 = rdChemReactions.ChemicalReaction(rxn.ToBinary())

  def test21CheckRawIters(self):
    RLABEL = "_MolFileRLabel"
    amine_rxn = '$RXN\n\n      ISIS     090220091541\n\n  2  1\n$MOL\n\n  -ISIS-  09020915412D\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n   -2.9083   -0.4708    0.0000 R#  0  0  0  0  0  0  0  0  0  1  0  0\n   -2.3995   -0.1771    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0\n   -2.4042    0.4125    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  2  3  2  0  0  0  0\nV    2 aldehyde\nM  RGP  1   1   1\nM  END\n$MOL\n\n  -ISIS-  09020915412D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    2.8375   -0.2500    0.0000 R#  0  0  0  0  0  0  0  0  0  3  0  0\n    3.3463    0.0438    0.0000 N   0  0  0  0  0  0  0  0  0  4  0  0\n  1  2  1  0  0  0  0\nV    2 amine\nM  RGP  1   1   2\nM  END\n$MOL\n\n  -ISIS-  09020915412D\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n   13.3088    0.9436    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0\n   13.8206    1.2321    0.0000 R#  0  0  0  0  0  0  0  0  0  1  0  0\n   13.3028    0.3561    0.0000 N   0  0  0  0  0  0  0  0  0  4  0  0\n   12.7911    0.0676    0.0000 R#  0  0  0  0  0  0  0  0  0  3  0  0\n  1  3  1  0  0  0  0\n  1  2  1  0  0  0  0\n  3  4  1  0  0  0  0\nM  RGP  2   2   1   4   2\nM  END\n'
    rxn = rdChemReactions.ReactionFromRxnBlock(amine_rxn)
    reactants = rxn.GetReactants()
    self.assertEqual(len(reactants), rxn.GetNumReactantTemplates())
    products = rxn.GetProducts()
    self.assertEqual(len(products), rxn.GetNumProductTemplates())
    agents = rxn.GetAgents()
    self.assertEqual(len(agents), rxn.GetNumAgentTemplates())

    for i in range(rxn.GetNumReactantTemplates()):
      p = rxn.GetReactantTemplate(i)
      mb1 = Chem.MolToMolBlock(p)
      mb2 = Chem.MolToMolBlock(reactants[i])
      self.assertEqual(mb1, mb2)

  def test22RunSingleReactant(self):
    # from
    # A Collection of Robust Organic Synthesis Reactions for In Silico Molecule Design
    # Markus Hartenfeller,*, Martin Eberle, Peter Meier, Cristina Nieto-Oberhuber,
    # Karl-Heinz Altmann, Gisbert Schneider, Edgar Jacoby, and Steffen Renner
    # Novartis Institutes for BioMedical Research, Novartis Pharma AG, Forum 1,
    # Novartis Campus, CH-4056 Basel, Switzerland Swiss Federal Institute of Technology (ETH)
    #  Zurich, Switzerland
    smirks_thiourea = "[N;$(N-[#6]):3]=[C;$(C=S):1].[N;$(N[#6]);!$(N=*);!$([N-]);!$(N#*);!$([ND3]);!$([ND4]);!$(N[O,N]);!$(N[C,S]=[S,O,N]):2]>>[N:3]-[C:1]-[N+0:2]"
    rxn = rdChemReactions.ReactionFromSmarts(smirks_thiourea)
    reagents = [Chem.MolFromSmiles(x) for x in ['C=CCN=C=S', 'NCc1ncc(Cl)cc1Br']]
    res = rxn.RunReactants(reagents)
    self.assertTrue(res)
    expected_result = [Chem.MolToSmiles(Chem.MolFromSmiles("C=CCNC(N)=S"))]
    expected_result.sort()
    sidechains_expected_result = [
      Chem.MolToSmiles(Chem.MolFromSmiles("[*:1]=S.[*:3]CC=C"), isomericSmiles=True)
    ]
    sidechains_nodummy_expected_result = [[0, [
      3,
    ], [
      1,
    ]], [3, [
      1,
    ], [
      2,
    ]]]
    sidechains_nodummy = []

    sidechains_expected_result.sort()

    for addDummy in [True, False]:
      res = rxn.RunReactant(reagents[0], 0)
      assert res
      result = []
      sidechains = []
      for match in res:
        for mol in match:
          result.append(Chem.MolToSmiles(mol, isomericSmiles=True))
          sidechain = rdChemReactions.ReduceProductToSideChains(mol, addDummy)
          sidechains.append(Chem.MolToSmiles(sidechain, isomericSmiles=True))
          if not addDummy:
            for atom in sidechain.GetAtoms():
              if atom.HasProp("_rgroupAtomMaps"):
                sidechains_nodummy.append([
                  atom.GetIdx(),
                  eval(atom.GetProp("_rgroupAtomMaps")),
                  eval(atom.GetProp("_rgroupBonds")),
                ])
          result.sort()
          sidechains.sort()

      if addDummy:
        self.assertEqual(result, expected_result)
        self.assertEqual(sidechains, sidechains_expected_result)
      else:
        self.assertEqual(sidechains_nodummy, sidechains_nodummy_expected_result)

    expected_result = [Chem.MolToSmiles(Chem.MolFromSmiles("NCNCc1ncc(Cl)cc1Br"))]
    expected_result.sort()
    sidechains_expected_result = [
      Chem.MolToSmiles(Chem.MolFromSmiles("[*:2]Cc1ncc(Cl)cc1Br"), isomericSmiles=True)
    ]
    sidechains_expected_result.sort()

    res = rxn.RunReactant(reagents[1], 1)
    result = []
    sidechains = []
    for match in res:
      for mol in match:
        result.append(Chem.MolToSmiles(mol, isomericSmiles=True))
        sidechains.append(
          Chem.MolToSmiles(rdChemReactions.ReduceProductToSideChains(mol), isomericSmiles=True))

    result.sort()
    self.assertEqual(result, expected_result)
    self.assertEqual(sidechains, sidechains_expected_result)

    self.assertFalse(rxn.RunReactant(reagents[0], 1))
    self.assertFalse(rxn.RunReactant(reagents[1], 0))

    # try a broken ring based side-chain
    sidechains_expected_result = ['c1ccc2c(c1)nc1n2CC[*:2]1']
    reactant = Chem.MolFromSmiles('c1ccc2c(c1)nc1n2CCN1')
    res = rxn.RunReactant(reactant, 1)
    result = []
    sidechains = []
    for match in res:
      for mol in match:
        result.append(Chem.MolToSmiles(mol, isomericSmiles=True))
        sidechains.append(
          Chem.MolToSmiles(rdChemReactions.ReduceProductToSideChains(mol), isomericSmiles=True))
        sidechain = rdChemReactions.ReduceProductToSideChains(mol, addDummyAtoms=False)

    self.assertEqual(sidechains, sidechains_expected_result)

  def test23CheckNonProduct(self):
    smirks_thiourea = "[N;$(N-[#6]):3]=[C;$(C=S):1].[N;$(N[#6]);!$(N=*);!$([N-]);!$(N#*);!$([ND3]);!$([ND4]);!$(N[O,N]);!$(N[C,S]=[S,O,N]):2]>>[N:3]-[C:1]-[N+0:2]"
    rxn = rdChemReactions.ReactionFromSmarts(smirks_thiourea)
    mol = Chem.MolFromSmiles("CCCCCCCC")
    m = rdChemReactions.ReduceProductToSideChains(mol)
    self.assertTrue(m.GetNumAtoms() == 0)
    mol = Chem.AddHs(mol)
    m = rdChemReactions.ReduceProductToSideChains(mol)
    self.assertTrue(m.GetNumAtoms() == 0)

  def testPreprocess(self):
    testFile = os.path.join(RDConfig.RDCodeDir, 'Chem', 'SimpleEnum', 'test_data', 'boronic1.rxn')
    rxn = rdChemReactions.ReactionFromRxnFile(testFile)
    rxn.Initialize()
    res = rdChemReactions.PreprocessReaction(rxn)
    self.assertEqual(res,
                     (0, 0, 2, 1, (((0, 'halogen.bromine.aromatic'), ), ((1, 'boronicacid'), ))))

  def testProperties(self):
    smirks_thiourea = "[N;$(N-[#6]):3]=[C;$(C=S):1].[N;$(N[#6]);!$(N=*);!$([N-]);!$(N#*);!$([ND3]);!$([ND4]);!$(N[O,N]);!$(N[C,S]=[S,O,N]):2]>>[N:3]-[C:1]-[N+0:2]"
    rxn = rdChemReactions.ReactionFromSmarts(smirks_thiourea)
    self.assertFalse(rxn.HasProp("fooprop"))
    rxn.SetProp("fooprop", "bar", computed=True)
    rxn.SetIntProp("intprop", 3)
    self.assertTrue(rxn.HasProp("fooprop"))
    self.assertTrue(rxn.HasProp("intprop"))
    self.assertEqual(rxn.GetIntProp("intprop"), 3)
    nrxn = rdChemReactions.ChemicalReaction(rxn.ToBinary())
    self.assertFalse(nrxn.HasProp("fooprop"))
    nrxn = rdChemReactions.ChemicalReaction(rxn.ToBinary(Chem.PropertyPickleOptions.AllProps))
    self.assertTrue(nrxn.HasProp("fooprop"))
    nrxn.ClearComputedProps()
    self.assertFalse(nrxn.HasProp("fooprop"))
    self.assertTrue(nrxn.HasProp("intprop"))
    self.assertEqual(nrxn.GetIntProp("intprop"), 3)

  def testRoundTripException(self):
    smarts = '[C:1]([C@:3]1([OH:24])[CH2:8][CH2:7][C@H:6]2[C@H:9]3[C@H:19]([C@@H:20]([F:22])[CH2:21][C@:4]12[CH3:5])[C@:17]1([CH3:18])[C:12](=[CH:13][C:14](=[O:23])[CH2:15][CH2:16]1)[CH:11]=[CH:10]3)#[CH:2].C(Cl)CCl.ClC1C=CC=C(C(OO)=[O:37])C=1.C(O)(C)(C)C>C(OCC)(=O)C>[C:1]([C@:3]1([OH:24])[CH2:8][CH2:7][C@H:6]2[C@H:9]3[C@H:19]([C@@H:20]([F:22])[CH2:21][C@:4]12[CH3:5])[C@:17]1([CH3:18])[C:12](=[CH:13][C:14](=[O:23])[CH2:15][CH2:16]1)[C@H:11]1[O:37][C@@H:10]31)#[CH:2]'
    rxn = rdChemReactions.ReactionFromSmarts(smarts)
    # this shouldn't throw an exception
    smarts = rdChemReactions.ReactionToSmarts(rxn)

  def testMaxProducts(self):
    smarts = "[c:1]1[c:2][c:3][c:4][c:5][c:6]1>>[c:1]1[c:2][c:3][c:4][c:5][c:6]1"
    rxn = rdChemReactions.ReactionFromSmarts(smarts)
    m = Chem.MolFromSmiles("c1ccccc1")
    prods = rxn.RunReactants([m])
    self.assertEqual(len(prods), 12)

    prods = rxn.RunReactants([m], 1)
    self.assertEqual(len(prods), 1)

  def testGitHub1868(self):
    fileN = os.path.join(self.dataDir, 'v3k.AmideBond.rxn')
    for i in range(100):
      _rxn = rdChemReactions.ReactionFromRxnFile(fileN)
      _rxn.Initialize()
      _reacts = [Chem.MolToSmarts(r) for r in _rxn.GetReactants()]
      _prods = [Chem.MolToSmarts(p) for p in _rxn.GetProducts()]

  @unittest.skipUnless(hasattr(rdChemReactions, 'ReactionFromPNGFile'),
                       "RDKit not built with iostreams support")
  def test_PNGMetadata(self):
    fname = os.path.join(self.dataDir, 'reaction1.smarts.png')
    rxn = rdChemReactions.ReactionFromPNGFile(fname)
    self.assertFalse(rxn is None)
    self.assertEqual(rxn.GetNumReactantTemplates(), 2)
    self.assertEqual(rxn.GetNumProductTemplates(), 1)

    fname = os.path.join(self.dataDir, 'reaction1.no_metadata.png')
    npng1 = rdChemReactions.ReactionMetadataToPNGFile(rxn, fname)
    png = open(fname, 'rb').read()
    npng2 = rdChemReactions.ReactionMetadataToPNGString(rxn, png)
    self.assertEqual(npng1, npng2)
    nrxn = rdChemReactions.ReactionFromPNGString(npng2)
    self.assertFalse(nrxn is None)
    self.assertEqual(nrxn.GetNumReactantTemplates(), 2)
    self.assertEqual(nrxn.GetNumProductTemplates(), 1)
    opts = [
      {
        'includePkl': True,
        'includeSmiles': False,
        'includeSmarts': False,
        'includeRxn': False
      },
      {
        'includePkl': False,
        'includeSmiles': True,
        'includeSmarts': False,
        'includeRxn': False
      },
      {
        'includePkl': False,
        'includeSmiles': False,
        'includeSmarts': True,
        'includeRxn': False
      },
      {
        'includePkl': False,
        'includeSmiles': False,
        'includeSmarts': False,
        'includeRxn': True
      },
    ]
    for opt in opts:
      npng = rdChemReactions.ReactionMetadataToPNGString(rxn, png, **opt)
      nrxn = rdChemReactions.ReactionFromPNGString(npng)
      self.assertFalse(nrxn is None)
      self.assertEqual(nrxn.GetNumReactantTemplates(), 2)
      self.assertEqual(nrxn.GetNumProductTemplates(), 1)


def _getProductCXSMILES(product):
  """
  Clear properties.

  Mapping properties show up in CXSMILES make validation less readable.
  """
  for a in product.GetAtoms():
    for k in a.GetPropsAsDict():
      a.ClearProp(k)
  return Chem.MolToCXSmiles(product)


def _reactAndSummarize(rxn_smarts, *smiles):
  """
  Run a reaction and combine the products in a single string.  

  Makes errors readable ish
  """
  rxn = rdChemReactions.ReactionFromSmarts(rxn_smarts)
  mols = [Chem.MolFromSmiles(s) for s in smiles]
  products = []
  for prods in rxn.RunReactants(mols):
    products.append(' + '.join(map(_getProductCXSMILES, prods)))
  products = ' OR '.join(products)
  return products


class StereoGroupTests(unittest.TestCase):

  def test_reaction_preserves_stereo(self):
    """
    StereoGroup atoms are in the reaction, but the reaction doesn't affect
    the chirality at the stereo centers
    -> preserve stereo group
    """
    reaction = '[C@:1]>>[C@:1]'
    reactants = ['F[C@H](Cl)Br |o1:1|', 'F[C@H](Cl)Br |&1:1|', 'FC(Cl)Br']
    for reactant in reactants:
      products = _reactAndSummarize(reaction, reactant)
      self.assertEqual(products, reactant)

  def test_reaction_ignores_stereo(self):
    """
    StereoGroup atoms are in the reaction, but the reaction doesn't specify the
    chirality at the stereo centers
    -> preserve stereo group
    """
    reaction = '[C:1]>>[C:1]'
    reactants = ['F[C@H](Cl)Br |o1:1|', 'F[C@H](Cl)Br |&1:1|', 'FC(Cl)Br']
    for reactant in reactants:
      products = _reactAndSummarize(reaction, reactant)
      self.assertEqual(products, reactant)

  def test_reaction_inverts_stereo(self):
    """
    StereoGroup atoms are in the reaction, and the reaction inverts the specified
    chirality at the stereo centers.
    -> preserve stereo group
    """
    reaction = '[C@:1]>>[C@@:1]'

    products = _reactAndSummarize(reaction, 'F[C@H](Cl)Br |o1:1|')
    self.assertEqual(products, 'F[C@H](Cl)Br |o1:1|')
    products = _reactAndSummarize(reaction, 'F[C@@H](Cl)Br |&1:1|')
    self.assertEqual(products, 'F[C@H](Cl)Br |&1:1|')
    products = _reactAndSummarize(reaction, 'FC(Cl)Br')
    self.assertEqual(products, 'FC(Cl)Br')

  def test_reaction_destroys_stereo(self):
    """
    StereoGroup atoms are in the reaction, and the reaction destroys the specified
    chirality at the stereo centers
    -> invalidate stereo center, preserve the rest of the stereo group.
    """
    reaction = '[C@:1]>>[C:1]'
    products = _reactAndSummarize(reaction, 'F[C@H](Cl)Br |o1:1|')
    self.assertEqual(products, 'FC(Cl)Br')
    products = _reactAndSummarize(reaction, 'F[C@@H](Cl)Br |&1:1|')
    self.assertEqual(products, 'FC(Cl)Br')
    products = _reactAndSummarize(reaction, 'FC(Cl)Br')
    self.assertEqual(products, 'FC(Cl)Br')

    reaction = '[C@:1]F>>[C:1]F'
    # Reaction destroys stereo (but preserves unaffected group
    products = _reactAndSummarize(reaction, 'F[C@H](Cl)[C@@H](Cl)Br |o1:1,&2:3|')
    self.assertEqual(products, 'FC(Cl)[C@H](Cl)Br |&1:3|')
    # Reaction destroys stereo (but preserves the rest of the group
    products = _reactAndSummarize(reaction, 'F[C@H](Cl)[C@@H](Cl)Br |&1:1,3|')
    self.assertEqual(products, 'FC(Cl)[C@H](Cl)Br |&1:3|')

  def test_reaction_defines_stereo(self):
    """
    StereoGroup atoms are in the reaction, and the reaction creates the specified
    chirality at the stereo centers
    -> remove the stereo center from 
    -> invalidate stereo group
    """
    products = _reactAndSummarize('[C:1]>>[C@@:1]', 'F[C@H](Cl)Br |o1:1|')
    self.assertEqual(products, 'F[C@@H](Cl)Br')
    products = _reactAndSummarize('[C:1]>>[C@@:1]', 'F[C@@H](Cl)Br |&1:1|')
    self.assertEqual(products, 'F[C@@H](Cl)Br')
    products = _reactAndSummarize('[C:1]>>[C@@:1]', 'FC(Cl)Br')
    self.assertEqual(products, 'F[C@@H](Cl)Br')

    # Remove group with defined stereo
    products = _reactAndSummarize('[C:1]F>>[C@@:1]F', 'F[C@H](Cl)[C@@H](Cl)Br |o1:1,&2:3|')
    self.assertEqual(products, 'F[C@@H](Cl)[C@H](Cl)Br |&1:3|')

    # Remove atoms with defined stereo from group
    products = _reactAndSummarize('[C:1]F>>[C@@:1]F', 'F[C@H](Cl)[C@@H](Cl)Br |o1:1,3|')
    self.assertEqual(products, 'F[C@@H](Cl)[C@H](Cl)Br |o1:3|')

  def test_stereogroup_is_spectator_to_reaction(self):
    """
    StereoGroup atoms are not in the reaction
    -> stereo group is unaffected
    """
    # 5a. Reaction preserves unrelated stereo
    products = _reactAndSummarize('[C@:1]F>>[C@:1]F', 'F[C@H](Cl)[C@@H](Cl)Br |o1:3|')
    self.assertEqual(products, 'F[C@H](Cl)[C@H](Cl)Br |o1:3|')
    # 5b. Reaction ignores unrelated stereo'
    products = _reactAndSummarize('[C:1]F>>[C:1]F', 'F[C@H](Cl)[C@@H](Cl)Br |o1:3|')
    self.assertEqual(products, 'F[C@H](Cl)[C@H](Cl)Br |o1:3|')
    # 5c. Reaction inverts unrelated stereo'
    products = _reactAndSummarize('[C@:1]F>>[C@@:1]F', 'F[C@H](Cl)[C@@H](Cl)Br |o1:3|')
    self.assertEqual(products, 'F[C@@H](Cl)[C@H](Cl)Br |o1:3|')
    # 5d. Reaction destroys unrelated stereo' 1:3|
    products = _reactAndSummarize('[C@:1]F>>[C:1]F', 'F[C@H](Cl)[C@@H](Cl)Br |o1:3|')
    self.assertEqual(products, 'FC(Cl)[C@H](Cl)Br |o1:3|')
    # 5e. Reaction assigns unrelated stereo'
    products = _reactAndSummarize('[C:1]F>>[C@@:1]F', 'F[C@H](Cl)[C@@H](Cl)Br |o1:3|')
    self.assertEqual(products, 'F[C@@H](Cl)[C@H](Cl)Br |o1:3|')

  def test_reaction_splits_stereogroup(self):
    """
    StereoGroup atoms are split into two products by the reaction
    -> Should the group be invalidated or trimmed?
    """
    products = _reactAndSummarize('[C:1]OO[C:2]>>[C:2]O.O[C:1]',
                                  'F[C@H](Cl)OO[C@@H](Cl)Br |o1:1,5|')
    # Two product sets, each with two mols:
    self.assertEqual(products.count('|o1:1|'), 4)

  def test_reaction_copies_stereogroup(self):
    """
    If multiple copies of an atom in StereoGroup show up in the product, they
    should all be part of the same product StereoGroup.
    """
    # Stereogroup atoms are in the reaction with multiple copies in the product
    products = _reactAndSummarize('[O:1].[C:2]=O>>[O:1][C:2][O:1]',
                                  'Cl[C@@H](Br)C[C@H](Br)CCO |&1:1,4|', 'CC(=O)C')
    # stereogroup manually checked, product SMILES assumed correct.
    self.assertEqual(products,
                     'CC(C)(OCC[C@H](Br)C[C@H](Cl)Br)OCC[C@H](Br)C[C@H](Cl)Br |&1:6,9,15,18|')

    # Stereogroup atoms are not in the reaction, but have multiple copies in the
    # product.
    products = _reactAndSummarize('[O:1].[C:2]=O>>[O:1][C:2][O:1]',
                                  'Cl[C@@H](Br)C[C@H](Br)CCO |&1:1,4|', 'CC(=O)C')
    # stereogroup manually checked, product SMILES assumed correct.
    self.assertEqual(products,
                     'CC(C)(OCC[C@H](Br)C[C@H](Cl)Br)OCC[C@H](Br)C[C@H](Cl)Br |&1:6,9,15,18|')

  def test_github(self):
    rxn = rdChemReactions.ReactionFromSmarts('[C:2][C:3]>>[C:2][C][*:3]')
    mol = Chem.MolFromSmiles("C/C=C/C")
    #  this shouldn't raise
    rxn.RunReactants([mol])

  def test_rxnblock_removehs(self):
    rxnblock = """$RXN
Dummy 0
  Dummy        0123456789

  1  1
$MOL

  Dummy   01234567892D

 10 10  0  0  0  0            999 V2000
    7.0222  -11.1783    0.0000 C   0  0  0  0  0  0  0  0  0  1  0  0
    8.0615  -11.7783    0.0000 O   0  0  0  0  0  0  0  0  0  2  0  0
    7.0222   -9.6783    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0
    5.7231   -8.9283    0.0000 C   0  0  0  0  0  0  0  0  0  4  0  0
    5.7231   -7.7283    0.0000 A   0  0  0  0  0  0  0  0  0  5  0  0
    4.4242   -9.6783    0.0000 C   0  0  0  0  0  0  0  0  0  6  0  0
    4.4242  -11.1783    0.0000 C   0  0  0  0  0  0  0  0  0  7  0  0
    3.3849  -11.7783    0.0000 A   0  0  0  0  0  0  0  0  0  8  0  0
    5.7231  -11.9283    0.0000 N   0  0  0  0  0  0  0  0  0  9  0  0
    5.7231  -13.1094    0.0000 H   0  0
  1  2  2  0  0  0  8
  1  3  1  0  0  0  8
  3  4  2  0  0  0  8
  4  5  1  0  0  0  2
  4  6  1  0  0  0  8
  6  7  2  0  0  0  8
  7  8  1  0  0  0  2
  7  9  1  0  0  0  8
  9  1  1  0  0  0  8
  9 10  1  0
M  SUB  1   9   2
M  END
$MOL

  Dummy   01234567892D

  9  9  0  0  0  0            999 V2000
   17.0447  -11.1783    0.0000 C   0  0  0  0  0  0  0  0  0  1  0  0
   18.0840  -11.7783    0.0000 O   0  0  0  0  0  0  0  0  0  2  0  0
   17.0447   -9.6783    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0
   15.7457   -8.9283    0.0000 C   0  0  0  0  0  0  0  0  0  4  0  0
   15.7457   -7.7283    0.0000 A   0  0  0  0  0  0  0  0  0  5  0  0
   14.4467   -9.6783    0.0000 C   0  0  0  0  0  0  0  0  0  6  0  0
   14.4467  -11.1783    0.0000 C   0  0  0  0  0  0  0  0  0  7  0  0
   13.4074  -11.7783    0.0000 A   0  0  0  0  0  0  0  0  0  8  0  0
   15.7457  -11.9283    0.0000 N   0  0  0  0  0  0  0  0  0  9  0  0
  1  2  1  0  0  0  8
  1  3  1  0  0  0  8
  3  4  2  0  0  0  8
  4  5  1  0  0  0  2
  4  6  1  0  0  0  8
  6  7  2  0  0  0  8
  7  8  1  0  0  0  2
  7  9  1  0  0  0  8
  9  1  2  0  0  0  8
M  END
"""

    mol = Chem.MolFromSmiles("c1(=O)nc([Cl])cc([F])[nH]1")
    rxn = rdChemReactions.ReactionFromRxnBlock(rxnblock)
    self.assertIsNotNone(rxn)
    prods = rxn.RunReactants((mol, ))
    # if the explicit hydrogen is not removed and the reactant template
    # is not sanitized, the reactant template is not aromatic and our
    # aromatic reactant won't match
    self.assertEqual(len(prods), 0)

    rxn = rdChemReactions.ReactionFromRxnBlock(rxnblock, removeHs=True, sanitize=True)
    self.assertIsNotNone(rxn)
    prods = rxn.RunReactants((mol, ))
    self.assertEqual(len(prods), 2)

  def testReactionInPlace(self):
    rxn = rdChemReactions.ReactionFromSmarts('[C:1][N:2][C:3]=[O:4]>>[C:1][O:2][C:3]=[O:4]')
    self.assertIsNotNone(rxn)
    reactant = Chem.MolFromSmiles('O=C(C)NCC')
    self.assertTrue(rxn.RunReactantInPlace(reactant))
    Chem.SanitizeMol(reactant)
    self.assertEqual(Chem.MolToSmiles(reactant), 'CCOC(C)=O')
    self.assertFalse(rxn.RunReactantInPlace(reactant))
    self.assertEqual(Chem.MolToSmiles(reactant), 'CCOC(C)=O')

    rxn = rdChemReactions.ReactionFromSmarts('CC[N:1]>>[N:1]')
    self.assertIsNotNone(rxn)
    reactant = Chem.MolFromSmiles('CCCN.Cl')
    self.assertTrue(rxn.RunReactantInPlace(reactant))
    Chem.SanitizeMol(reactant)
    self.assertEqual(Chem.MolToSmiles(reactant), 'N')

    reactant = Chem.MolFromSmiles('CCCN.Cl')
    self.assertTrue(rxn.RunReactantInPlace(reactant, removeUnmatchedAtoms=False))
    Chem.SanitizeMol(reactant)
    self.assertEqual(Chem.MolToSmiles(reactant), 'C.Cl.N')

  def testGithub4651(self):
    mol_sulfonylchloride = Chem.MolFromSmiles("Nc1c(CCCSNCC)cc(cc1)S(=O)(=O)Cl")
    mol_amine = Chem.MolFromSmiles("Nc2cc(C)on2")
    mol_sulfonamide = Chem.MolFromSmiles("CCNSCCCc1cc(S(=O)(=O)Nc2cc(C)on2)ccc1N")

    smirks_fwd = (
      "[S;$(S(=O)(=O)[C,c,N]):1](=O)(=O)(-[Cl])"
      "."
      "[N;$([N&H2&D1,N&H1&D2])"
      # N aliphatic and not aromatic bond to carbon
      "&$(N(-&!@[#6]))"
      "&!$(N-C=[O,N,S]):2]"
      ">>"
      "[S:1](=O)(=O)-[N+0:2]")

    smirks_reverse = ">>".join(smirks_fwd.split(">>")[::-1])

    reaction_fwd = rdChemReactions.ReactionFromSmarts(smirks_fwd)
    reaction_reverse = rdChemReactions.ReactionFromSmarts(smirks_reverse)

    product = reaction_fwd.RunReactants((mol_sulfonylchloride, mol_amine))[0][0]

    reagent_sulfonylchloride, reagent_amine = reaction_reverse.RunReactants((mol_sulfonamide, ))[0]

    # trigger bug
    with self.assertRaises(RuntimeError):
      reaction_fwd.RunReactants((reagent_sulfonylchloride, reagent_amine))

    # The second call used to deadlock
    with self.assertRaises(RuntimeError):
      reaction_fwd.RunReactants((reagent_sulfonylchloride, reagent_amine))

  def testV3000Reactions(self):
    reaction = rdChemReactions.ReactionFromSmarts(
      '[cH:1]1[cH:2][cH:3][cH:4][cH:5][c:6]1-[Br].[#6:7]B(O)O>>[cH:1]1[cH:2][cH:3][cH:4][cH:5][c:6]1-[#6:7]'
    )
    self.assertIsNotNone(reaction)

    mb1 = rdChemReactions.ReactionToRxnBlock(reaction, forceV3000=True)
    mb2 = rdChemReactions.ReactionToV3KRxnBlock(reaction)
    self.assertEqual(mb1, mb2)
    reaction2 = rdChemReactions.ReactionFromRxnBlock(mb1)
    self.assertIsNotNone(reaction2)

    self.assertEqual(reaction.GetNumReactantTemplates(), reaction2.GetNumReactantTemplates())

  def testGithub6138(self):
    mol = Chem.MolFromSmiles("COc1ccccc1Oc1nc(Nc2cc(C)[nH]n2)cc2ccccc12")
    rxn = AllChem.ReactionFromSmarts(
      "([c:1]:[n&H1&+0&D2:3]:[n:2])>>([c:1]:[n&H0&+0&D3:3](:[n:2])-C1-C-C-C-C-O-1)")

    def run(r):
      return Chem.MolToSmiles(r.RunReactants((mol, ))[0][0])

    rxn_reloaded = pickle.loads(pickle.dumps(rxn))

    res1 = Chem.MolToSmiles(rxn.RunReactants((mol, ))[0][0])
    res2 = Chem.MolToSmiles(rxn_reloaded.RunReactants((mol, ))[0][0])
    rxn_reloaded_after_use = pickle.loads(pickle.dumps(rxn))
    res3 = Chem.MolToSmiles(rxn_reloaded_after_use.RunReactants((mol, ))[0][0])
    self.assertEqual(res1, res2)
    self.assertEqual(res1, res3)

  def testGithub6211(self):
    rxn = AllChem.ReactionFromSmarts("[C:1][C@:2]([N:3])[O:4]>>[C:1][C@@:2]([N:3])[O:4]")
    mol = Chem.MolFromSmiles("CC[C@@H](N)O")
    self.assertEqual(Chem.MolToSmiles(rxn.RunReactants((mol, ))[0][0]), "CC[C@H](N)O")
    rxn.GetSubstructParams().useChirality = True
    self.assertEqual(len(rxn.RunReactants((mol, ))), 0)

  def testMrvBlockContainsReaction(self):
    fn1 = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'MarvinParse', 'test_data',
                       'aspirin.mrv')
    with open(fn1, 'r') as inf:
      ind1 = inf.read()
    fn2 = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'MarvinParse', 'test_data',
                       'aspirineSynthesisWithAttributes.mrv')
    with open(fn2, 'r') as inf:
      ind2 = inf.read()

    self.assertFalse(rdChemReactions.MrvFileIsReaction(fn1))
    self.assertTrue(rdChemReactions.MrvFileIsReaction(fn2))

    self.assertFalse(rdChemReactions.MrvBlockIsReaction(ind1))
    self.assertTrue(rdChemReactions.MrvBlockIsReaction(ind2))

  def testMrvOutput(self):
    fn2 = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'MarvinParse', 'test_data',
                       'aspirineSynthesisWithAttributes.mrv')
    rxn = rdChemReactions.ReactionFromMrvFile(fn2)
    self.assertIsNotNone(rxn)
    rxnb = rdChemReactions.ReactionToMrvBlock(rxn)
    self.assertTrue('<reaction>' in rxnb)

    fName = tempfile.NamedTemporaryFile(suffix='.mrv').name
    self.assertFalse(os.path.exists(fName))
    rdChemReactions.ReactionToMrvFile(rxn, fName)
    self.assertTrue(os.path.exists(fName))
    os.unlink(fName)

  def testCDXML(self):
    fname = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'test_data', 'CDXML', 'rxn2.cdxml')
    rxns = AllChem.ReactionsFromCDXMLFile(fname)
    self.assertEqual(len(rxns), 1)
    self.assertEqual(
      AllChem.ReactionToSmarts(rxns[0]),
      "[#6&D2:2]1:[#6&D2:3]:[#6&D2:4]:[#6&D3:1](:[#6&D2:5]:[#6&D2:6]:1)-[#17&D1].[#6&D3](-[#5&D2]-[#6&D3:7]1:[#6&D2:8]:[#6&D2:9]:[#6&D2:10]:[#6&D2:11]:[#6&D2:12]:1)(-[#8&D1])-[#8&D1]>>[#6&D2:1]1:[#6&D2:5]:[#6&D3:6](:[#6&D2:2]:[#6&D2:3]:[#6&D2:4]:1)-[#6&D3:7]1:[#6&D2:8]:[#6&D2:9]:[#6&D2:10]:[#6&D2:11]:[#6&D2:12]:1"
    )

    rxns = AllChem.ReactionsFromCDXMLFile('does-not-exist.cdxml')
    self.assertEqual(len(rxns), 0)

    with open(fname, 'r') as inf:
      cdxml = inf.read()
    rxns = AllChem.ReactionsFromCDXMLBlock(cdxml)
    self.assertEqual(len(rxns), 1)
    self.assertEqual(
      AllChem.ReactionToSmarts(rxns[0]),
      "[#6&D2:2]1:[#6&D2:3]:[#6&D2:4]:[#6&D3:1](:[#6&D2:5]:[#6&D2:6]:1)-[#17&D1].[#6&D3](-[#5&D2]-[#6&D3:7]1:[#6&D2:8]:[#6&D2:9]:[#6&D2:10]:[#6&D2:11]:[#6&D2:12]:1)(-[#8&D1])-[#8&D1]>>[#6&D2:1]1:[#6&D2:5]:[#6&D3:6](:[#6&D2:2]:[#6&D2:3]:[#6&D2:4]:1)-[#6&D3:7]1:[#6&D2:8]:[#6&D2:9]:[#6&D2:10]:[#6&D2:11]:[#6&D2:12]:1"
    )

    rxns = AllChem.ReactionsFromCDXMLBlock('')
    self.assertEqual(len(rxns), 0)

  def testSanitizeRxnAsMols(self):
    rxn = AllChem.ReactionFromSmarts("C1=CC=CC=C1>CN(=O)=O>C1=CC=CC=N1", useSmiles=True)
    self.assertIsNotNone(rxn)
    self.assertFalse(rxn.GetReactantTemplate(0).GetBondWithIdx(0).GetIsAromatic())
    self.assertFalse(rxn.GetProductTemplate(0).GetBondWithIdx(0).GetIsAromatic())
    self.assertEqual(rxn.GetAgentTemplate(0).GetAtomWithIdx(1).GetFormalCharge(), 0)

    AllChem.SanitizeRxnAsMols(rxn)
    self.assertTrue(rxn.GetReactantTemplate(0).GetBondWithIdx(0).GetIsAromatic())
    self.assertTrue(rxn.GetProductTemplate(0).GetBondWithIdx(0).GetIsAromatic())
    self.assertEqual(rxn.GetAgentTemplate(0).GetAtomWithIdx(1).GetFormalCharge(), 1)

    rxn = AllChem.ReactionFromSmarts("C1=CC=CC=C1>CN(=O)=O>C1=CC=CC=N1", useSmiles=True)
    self.assertIsNotNone(rxn)
    self.assertFalse(rxn.GetReactantTemplate(0).GetBondWithIdx(0).GetIsAromatic())
    self.assertFalse(rxn.GetProductTemplate(0).GetBondWithIdx(0).GetIsAromatic())
    self.assertEqual(rxn.GetAgentTemplate(0).GetAtomWithIdx(1).GetFormalCharge(), 0)

    AllChem.SanitizeRxnAsMols(rxn, Chem.SanitizeFlags.SANITIZE_CLEANUP)
    self.assertFalse(rxn.GetReactantTemplate(0).GetBondWithIdx(0).GetIsAromatic())
    self.assertFalse(rxn.GetProductTemplate(0).GetBondWithIdx(0).GetIsAromatic())
    self.assertEqual(rxn.GetAgentTemplate(0).GetAtomWithIdx(1).GetFormalCharge(), 1)

  def testSmilesWriteParams(self):
    rxn = AllChem.ReactionFromSmarts(
      "[C:1]-[C:2].[NH3:3]->[Fe:4]-[NH2:5]>>[C:1]=[C:2].[NH3:3]->[Fe:4]-[NH2:5]")
    self.assertIsNotNone(rxn)
    params = AllChem.SmilesWriteParams()
    self.assertEqual(
      AllChem.ReactionToSmiles(rxn, params),
      "[CH3:1][CH3:2].[NH3:3]->[Fe:4][NH2:5]>>[CH2:1]=[CH2:2].[NH3:3]->[Fe:4][NH2:5]")
    self.assertEqual(
      AllChem.ReactionToSmarts(rxn, params),
      "[C:1]-[C:2].[N&H3:3]->[#26:4]-[N&H2:5]>>[C:1]=[C:2].[N&H3:3]->[#26:4]-[N&H2:5]")
    params.includeDativeBonds = False
    self.assertEqual(AllChem.ReactionToSmiles(rxn, params),
                     "[CH3:1][CH3:2].[NH3:3][Fe:4][NH2:5]>>[CH2:1]=[CH2:2].[NH3:3][Fe:4][NH2:5]")
    self.assertEqual(
      AllChem.ReactionToSmarts(rxn, params),
      "[C:1]-[C:2].[N&H3:3]-[#26:4]-[N&H2:5]>>[C:1]=[C:2].[N&H3:3]-[#26:4]-[N&H2:5]")


if __name__ == '__main__':
  unittest.main(verbosity=True)
