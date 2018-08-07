#  $Id$
#
#  Copyright (c) 2007-2014, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
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
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
from __future__ import print_function

import unittest, doctest
import os, sys
from rdkit.six import exec_
from rdkit.six.moves import cPickle

from rdkit import rdBase
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit import Geometry
from rdkit import RDConfig
from rdkit.Chem.SimpleEnum import Enumerator


def feq(v1, v2, tol2=1e-4):
  return abs(v1 - v2) <= tol2


def ptEq(pt1, pt2, tol=1e-4):
  return feq(pt1.x, pt2.x, tol) and feq(pt1.y, pt2.y, tol) and feq(pt1.z, pt2.z, tol)

# Boost functions are NOT found by doctest, this "fixes" them
#  by adding the doctests to a fake module
import imp
TestPreprocess = imp.new_module("TestPreprocess")
code = """
from rdkit.Chem import rdChemReactions
def PreprocessReaction(*a, **kw):
    '''%s
    '''
    return rdChemReactions.PreprocessReaction(*a, **kw)
""" % "\n".join([x.lstrip() for x in rdChemReactions.PreprocessReaction.__doc__.split("\n")])
exec_(code, TestPreprocess.__dict__)


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
      ValueError,
      lambda x='[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]>>CC': rdChemReactions.ReactionFromSmarts(x))

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
    pkl = cPickle.dumps(rxn)
    rxn = cPickle.loads(pkl)
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
    querydefs = {k.lower(): v
                 for k, v in FilterCatalog.GetFlattenedFunctionalGroupHierarchy().items()}

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
    self.failUnless(rxn)
    rxn.Initialize()
    self.failUnless(rxn.GetNumReactantTemplates() == 3)
    self.failUnless(rxn.GetNumProductTemplates() == 1)
    self.failUnless(rxn.GetNumAgentTemplates() == 3)

    rxn.RemoveUnmappedReactantTemplates()
    rxn.RemoveUnmappedProductTemplates()

    self.failUnless(rxn.GetNumReactantTemplates() == 2)
    self.failUnless(rxn.GetNumProductTemplates() == 1)
    self.failUnless(rxn.GetNumAgentTemplates() == 4)

    rxn = rdChemReactions.ReactionFromSmarts("[C:1]=[O:2].[N:3].C(=O)O>>[N:3]~[C:1]=[O:2].[OH2]")
    self.failUnless(rxn)
    rxn.Initialize()
    self.failUnless(rxn.GetNumReactantTemplates() == 3)
    self.failUnless(rxn.GetNumProductTemplates() == 2)
    self.failUnless(rxn.GetNumAgentTemplates() == 0)

    agentList = []
    rxn.RemoveUnmappedReactantTemplates(moveToAgentTemplates=False, targetList=agentList)
    rxn.RemoveUnmappedProductTemplates(targetList=agentList)

    self.failUnless(rxn.GetNumReactantTemplates() == 2)
    self.failUnless(rxn.GetNumProductTemplates() == 1)
    self.failUnless(rxn.GetNumAgentTemplates() == 1)
    self.failUnless(len(agentList) == 2)

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
    self.assertEquals(res, res2)

    # currently ToBinary does not save atom props
    # rxn2 = rdChemReactions.ChemicalReaction(rxn.ToBinary())

  def test21CheckRawIters(self):
    RLABEL = "_MolFileRLabel"
    amine_rxn = '$RXN\n\n      ISIS     090220091541\n\n  2  1\n$MOL\n\n  -ISIS-  09020915412D\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n   -2.9083   -0.4708    0.0000 R#  0  0  0  0  0  0  0  0  0  1  0  0\n   -2.3995   -0.1771    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0\n   -2.4042    0.4125    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  2  3  2  0  0  0  0\nV    2 aldehyde\nM  RGP  1   1   1\nM  END\n$MOL\n\n  -ISIS-  09020915412D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    2.8375   -0.2500    0.0000 R#  0  0  0  0  0  0  0  0  0  3  0  0\n    3.3463    0.0438    0.0000 N   0  0  0  0  0  0  0  0  0  4  0  0\n  1  2  1  0  0  0  0\nV    2 amine\nM  RGP  1   1   2\nM  END\n$MOL\n\n  -ISIS-  09020915412D\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n   13.3088    0.9436    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0\n   13.8206    1.2321    0.0000 R#  0  0  0  0  0  0  0  0  0  1  0  0\n   13.3028    0.3561    0.0000 N   0  0  0  0  0  0  0  0  0  4  0  0\n   12.7911    0.0676    0.0000 R#  0  0  0  0  0  0  0  0  0  3  0  0\n  1  3  1  0  0  0  0\n  1  2  1  0  0  0  0\n  3  4  1  0  0  0  0\nM  RGP  2   2   1   4   2\nM  END\n'
    rxn = rdChemReactions.ReactionFromRxnBlock(amine_rxn)
    reactants = rxn.GetReactants()
    self.assertEquals(len(reactants), rxn.GetNumReactantTemplates())
    products = rxn.GetProducts()
    self.assertEquals(len(products), rxn.GetNumProductTemplates())
    agents = rxn.GetAgents()
    self.assertEquals(len(agents), rxn.GetNumAgentTemplates())

    for i in range(rxn.GetNumReactantTemplates()):
      p = rxn.GetReactantTemplate(i)
      mb1 = Chem.MolToMolBlock(p)
      mb2 = Chem.MolToMolBlock(reactants[i])
      self.assertEquals(mb1, mb2)

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
    sidechains_expected_result = [Chem.MolToSmiles(
      Chem.MolFromSmiles("[*:1]=S.[*:3]CC=C"), isomericSmiles=True)]
    sidechains_nodummy_expected_result = [[0, [3, ], [1, ]], [3, [1, ], [2, ]]]
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
                sidechains_nodummy.append([atom.GetIdx(),
                                           eval(atom.GetProp("_rgroupAtomMaps")),
                                           eval(atom.GetProp("_rgroupBonds")), ])
          result.sort()
          sidechains.sort()

      if addDummy:
        self.assertEquals(result, expected_result)
        self.assertEquals(sidechains, sidechains_expected_result)
      else:
        self.assertEquals(sidechains_nodummy, sidechains_nodummy_expected_result)

    expected_result = [Chem.MolToSmiles(Chem.MolFromSmiles("NCNCc1ncc(Cl)cc1Br"))]
    expected_result.sort()
    sidechains_expected_result = [Chem.MolToSmiles(
      Chem.MolFromSmiles("[*:2]Cc1ncc(Cl)cc1Br"), isomericSmiles=True)]
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
    self.assertEquals(result, expected_result)
    self.assertEquals(sidechains, sidechains_expected_result)

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

    self.assertEquals(sidechains, sidechains_expected_result)

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
    self.assertEquals(res, (0, 0, 2, 1, (((0, 'halogen.bromine.aromatic'), ), (
      (1, 'boronicacid'), ))))

  def testProperties(self):
    smirks_thiourea = "[N;$(N-[#6]):3]=[C;$(C=S):1].[N;$(N[#6]);!$(N=*);!$([N-]);!$(N#*);!$([ND3]);!$([ND4]);!$(N[O,N]);!$(N[C,S]=[S,O,N]):2]>>[N:3]-[C:1]-[N+0:2]"
    rxn = rdChemReactions.ReactionFromSmarts(smirks_thiourea)
    self.assertFalse(rxn.HasProp("fooprop"))
    rxn.SetProp("fooprop","bar",computed=True)
    rxn.SetIntProp("intprop",3)
    self.assertTrue(rxn.HasProp("fooprop"))
    self.assertTrue(rxn.HasProp("intprop"))
    self.assertEquals(rxn.GetIntProp("intprop"),3)
    nrxn = rdChemReactions.ChemicalReaction(rxn.ToBinary())
    self.assertFalse(nrxn.HasProp("fooprop"))
    nrxn = rdChemReactions.ChemicalReaction(rxn.ToBinary(Chem.PropertyPickleOptions.AllProps))
    self.assertTrue(nrxn.HasProp("fooprop"))
    nrxn.ClearComputedProps()
    self.assertFalse(nrxn.HasProp("fooprop"))
    self.assertTrue(nrxn.HasProp("intprop"))
    self.assertEquals(nrxn.GetIntProp("intprop"),3)

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
    
    prods = rxn.RunReactants([m],1)
    self.assertEqual(len(prods), 1)
    
if __name__ == '__main__':
  unittest.main(verbosity=True)
