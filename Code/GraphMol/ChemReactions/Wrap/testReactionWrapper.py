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

import unittest
import os,sys

from rdkit.six.moves import cPickle

from rdkit import rdBase
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit import Geometry
from rdkit import RDConfig

def feq(v1,v2,tol2=1e-4):
  return abs(v1-v2)<=tol2

def ptEq(pt1, pt2, tol=1e-4):
  return feq(pt1.x,pt2.x,tol) and feq(pt1.y,pt2.y,tol) and feq(pt1.z,pt2.z,tol)

class TestCase(unittest.TestCase) :
  def setUp(self):
    self.dataDir = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','ChemReactions','testData')


  def test1Basics(self):
    rxn = rdChemReactions.ChemicalReaction()
    self.assertTrue(rxn.GetNumReactantTemplates()==0)
    self.assertTrue(rxn.GetNumProductTemplates()==0)

    r1= Chem.MolFromSmarts('[C:1](=[O:2])O')
    rxn.AddReactantTemplate(r1)
    self.assertTrue(rxn.GetNumReactantTemplates()==1)

    r1= Chem.MolFromSmarts('[N:3]')
    rxn.AddReactantTemplate(r1)
    self.assertTrue(rxn.GetNumReactantTemplates()==2)

    r1= Chem.MolFromSmarts('[C:1](=[O:2])[N:3]')
    rxn.AddProductTemplate(r1)
    self.assertTrue(rxn.GetNumProductTemplates()==1)

    reacts = (Chem.MolFromSmiles('C(=O)O'),Chem.MolFromSmiles('N'))
    ps = rxn.RunReactants(reacts)
    self.assertTrue(len(ps)==1)
    self.assertTrue(len(ps[0])==1)
    self.assertTrue(ps[0][0].GetNumAtoms()==3)
        
    ps = rxn.RunReactants(list(reacts))
    self.assertTrue(len(ps)==1)
    self.assertTrue(len(ps[0])==1)
    self.assertTrue(ps[0][0].GetNumAtoms()==3)    

  def test2DaylightParser(self):
    rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]')
    self.assertTrue(rxn)
    self.assertTrue(rxn.GetNumReactantTemplates()==2)
    self.assertTrue(rxn.GetNumProductTemplates()==1)
    self.assertTrue(rxn._getImplicitPropertiesFlag())
    
    reacts = (Chem.MolFromSmiles('C(=O)O'),Chem.MolFromSmiles('N'))
    ps = rxn.RunReactants(reacts)
    self.assertTrue(len(ps)==1)
    self.assertTrue(len(ps[0])==1)
    self.assertTrue(ps[0][0].GetNumAtoms()==3)  
    
    reacts = (Chem.MolFromSmiles('CC(=O)OC'),Chem.MolFromSmiles('CN'))
    ps = rxn.RunReactants(reacts)
    self.assertTrue(len(ps)==1)
    self.assertTrue(len(ps[0])==1)
    self.assertTrue(ps[0][0].GetNumAtoms()==5)  
    
  def test3MDLParsers(self):
    fileN = os.path.join(self.dataDir,'AmideBond.rxn')
    rxn = rdChemReactions.ReactionFromRxnFile(fileN) 
    self.assertTrue(rxn)
    self.assertFalse(rxn._getImplicitPropertiesFlag())

    self.assertTrue(rxn.GetNumReactantTemplates()==2)
    self.assertTrue(rxn.GetNumProductTemplates()==1)
    
    reacts = (Chem.MolFromSmiles('C(=O)O'),Chem.MolFromSmiles('N'))
    ps = rxn.RunReactants(reacts)
    self.assertTrue(len(ps)==1)
    self.assertTrue(len(ps[0])==1)
    self.assertTrue(ps[0][0].GetNumAtoms()==3)  
    
    with open(fileN, 'r') as rxnF:
      rxnBlock = rxnF.read()
    rxn = rdChemReactions.ReactionFromRxnBlock(rxnBlock) 
    self.assertTrue(rxn)

    self.assertTrue(rxn.GetNumReactantTemplates()==2)
    self.assertTrue(rxn.GetNumProductTemplates()==1)
    
    reacts = (Chem.MolFromSmiles('C(=O)O'),Chem.MolFromSmiles('N'))
    ps = rxn.RunReactants(reacts)
    self.assertTrue(len(ps)==1)
    self.assertTrue(len(ps[0])==1)
    self.assertTrue(ps[0][0].GetNumAtoms()==3)  

  def test4ErrorHandling(self):
    self.assertRaises(ValueError,lambda x='[C:1](=[O:2])Q.[N:3]>>[C:1](=[O:2])[N:3]':rdChemReactions.ReactionFromSmarts(x))
    self.assertRaises(ValueError,lambda x='[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]Q':rdChemReactions.ReactionFromSmarts(x))
    self.assertRaises(ValueError,lambda x='[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]>>CC':rdChemReactions.ReactionFromSmarts(x))

    block="""$RXN

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
    self.assertRaises(ValueError,lambda x=block:rdChemReactions.ReactionFromRxnBlock(x))
        
    block="""$RXN

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

    block="""$RXN

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
    self.assertTrue(rxn.Validate()==(0,0))

    rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:1])O.[N:3]>>[C:1](=[O:2])[N:3]')
    self.assertTrue(rxn)
    self.assertTrue(rxn.Validate()==(1,1))

    rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])[O:4].[N:3]>>[C:1](=[O:2])[N:3]')
    self.assertTrue(rxn)
    self.assertTrue(rxn.Validate()==(1,0))

    rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3][C:5]')
    self.assertTrue(rxn)
    self.assertTrue(rxn.Validate()==(1,0))

  def test6Exceptions(self):
    rxn = rdChemReactions.ReactionFromSmarts('[C:1]Cl>>[C:1]')
    self.assertTrue(rxn)
    self.assertRaises(ValueError,lambda x=rxn:x.RunReactants(()))
    self.assertRaises(ValueError,lambda x=rxn:x.RunReactants((Chem.MolFromSmiles('CC'),Chem.MolFromSmiles('C'))))
    ps=rxn.RunReactants((Chem.MolFromSmiles('CCCl'),))
    self.assertTrue(len(ps)==1)
    self.assertTrue(len(ps[0])==1)

  def _test7Leak(self):
    rxn = rdChemReactions.ReactionFromSmarts('[C:1]Cl>>[C:1]')
    self.assertTrue(rxn)
    print('running: ')
    for i in range(1e5):
      ps=rxn.RunReactants((Chem.MolFromSmiles('CCCl'),))
      self.assertTrue(len(ps)==1)
      self.assertTrue(len(ps[0])==1)
      if not i%1000: print(i)

  def test8Properties(self):
    rxn = rdChemReactions.ReactionFromSmarts('[O:1]>>[O:1][3#0]')
    self.assertTrue(rxn)
    ps=rxn.RunReactants((Chem.MolFromSmiles('CO'),))
    self.assertTrue(len(ps)==1)
    self.assertTrue(len(ps[0])==1)
    Chem.SanitizeMol(ps[0][0])
    self.assertEqual(ps[0][0].GetAtomWithIdx(1).GetIsotope(),3);

  def test9AromaticityTransfer(self):
    # this was issue 2664121
    mol = Chem.MolFromSmiles('c1ccc(C2C3(Cc4c(cccc4)C2)CCCC3)cc1')
    rxn = rdChemReactions.ReactionFromSmarts('[A:1]1~[*:2]~[*:3]~[*:4]~[*:5]~[A:6]-;@1>>[*:1]~[*:2]~[*:3]~[*:4]~[*:5]~[*:6]')
    products = rxn.RunReactants([mol])
    self.assertEqual(len(products),6)
    for p in products:
      self.assertEqual(len(p),1)
      Chem.SanitizeMol(p[0])

  def test10DotSeparation(self):
    # 08/05/14
    # This test is changed due to a new behavior of the smarts
    # reaction parser which now allows using parenthesis in products
    # as well. original smiles: '[C:1]1[O:2][N:3]1>>[C:1]1[O:2].[N:3]1'
    rxn = rdChemReactions.ReactionFromSmarts('[C:1]1[O:2][N:3]1>>([C:1]1[O:2].[N:3]1)')
    mol = Chem.MolFromSmiles('C1ON1')
    products = rxn.RunReactants([mol])
    self.assertEqual(len(products),1)
    for p in products:
      self.assertEqual(len(p),1)
      self.assertEqual(p[0].GetNumAtoms(),3)
      self.assertEqual(p[0].GetNumBonds(),2)

  def test11ImplicitProperties(self):
    rxn = rdChemReactions.ReactionFromSmarts('[C:1]O>>[C:1]')
    mol = Chem.MolFromSmiles('CCO')
    products = rxn.RunReactants([mol])
    self.assertEqual(len(products),1)
    for p in products:
      self.assertEqual(len(p),1)
      self.assertEqual(Chem.MolToSmiles(p[0]),'CC')
    mol2 = Chem.MolFromSmiles('C[CH-]O')
    products = rxn.RunReactants([mol2])
    self.assertEqual(len(products),1)
    for p in products:
      self.assertEqual(len(p),1)
      self.assertEqual(Chem.MolToSmiles(p[0]),'[CH2-]C')

    rxn._setImplicitPropertiesFlag(False)
    products = rxn.RunReactants([mol])
    self.assertEqual(len(products),1)
    for p in products:
      self.assertEqual(len(p),1)
      self.assertEqual(Chem.MolToSmiles(p[0]),'CC')
    products = rxn.RunReactants([mol2])
    self.assertEqual(len(products),1)
    for p in products:
      self.assertEqual(len(p),1)
      self.assertEqual(Chem.MolToSmiles(p[0]),'CC')


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
    self.assertEqual(len(products),1)
    for p in products:
      self.assertEqual(len(p),1)
      self.assertEqual(p[0].GetNumAtoms(),3)
      self.assertEqual(p[0].GetNumBonds(),2)

    rxn = rdChemReactions.ChemicalReaction(rxn.ToBinary())
    products = rxn.RunReactants([mol])
    self.assertEqual(len(products),1)
    for p in products:
      self.assertEqual(len(p),1)
      self.assertEqual(p[0].GetNumAtoms(),3)
      self.assertEqual(p[0].GetNumBonds(),2)

  def test13GetTemplates(self):
    rxn = rdChemReactions.ReactionFromSmarts('[C:1]1[O:2][N:3]1>>[C:1][O:2].[N:3]')
    r1 = rxn.GetReactantTemplate(0)
    sma=Chem.MolToSmarts(r1)
    self.assertEqual(sma,'[C:1]1-,:[O:2]-,:[N:3]-,:1')
    p1 = rxn.GetProductTemplate(0)
    sma=Chem.MolToSmarts(p1)
    self.assertEqual(sma,'[C:1]-,:[O:2]')
    
    p2 = rxn.GetProductTemplate(1)
    sma=Chem.MolToSmarts(p2)
    self.assertEqual(sma,'[N:3]')

    self.assertRaises(ValueError,lambda :rxn.GetProductTemplate(2))
    self.assertRaises(ValueError,lambda :rxn.GetReactantTemplate(1))

  def test14Matchers(self):
    rxn = rdChemReactions.ReactionFromSmarts('[C;!$(C(-O)-O):1](=[O:2])[O;H,-1].[N;!H0:3]>>[C:1](=[O:2])[N:3]')
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
    rxn = rdChemReactions.ReactionFromSmarts('[{amine}:1]>>[*:1]-C',
                                             replacements={'{amine}':'$([N;!H0;$(N-[#6]);!$(N-[!#6;!#1]);!$(N-C=[O,N,S])])'})
    self.assertTrue(rxn)
    rxn.Initialize()
    reactants = (Chem.MolFromSmiles('CCN'),)
    ps = rxn.RunReactants(reactants)
    self.assertEqual(len(ps),1)
    self.assertEqual(len(ps[0]),1)
    self.assertEqual(ps[0][0].GetNumAtoms(),4)
    
  def test16GetReactingAtoms(self):
    rxn = rdChemReactions.ReactionFromSmarts("[O:1][C:2].[N:3]>>[N:1][C:2].[N:3]")
    self.assertTrue(rxn)
    rxn.Initialize()
    rAs = rxn.GetReactingAtoms()
    self.assertEqual(len(rAs),2)
    self.assertEqual(len(rAs[0]),1)    
    self.assertEqual(len(rAs[1]),0)    

    rxn = rdChemReactions.ReactionFromSmarts("[O:1]C>>[O:1]C")
    self.assertTrue(rxn)
    rxn.Initialize()
    rAs = rxn.GetReactingAtoms()
    self.assertEqual(len(rAs),1)
    self.assertEqual(len(rAs[0]),2)    
    rAs = rxn.GetReactingAtoms(True)
    self.assertEqual(len(rAs),1)
    self.assertEqual(len(rAs[0]),1)   

  def test17AddRecursiveQueriesToReaction(self):
    rxn = rdChemReactions.ReactionFromSmarts("[C:1][O:2].[N:3]>>[C:1][N:2]")
    self.assertTrue(rxn)
    rxn.Initialize()
    qs = {'aliphatic':Chem.MolFromSmiles('CC')}
    rxn.GetReactantTemplate(0).GetAtomWithIdx(0).SetProp('query', 'aliphatic')
    rxn.AddRecursiveQueriesToReaction(qs,'query')
    q = rxn.GetReactantTemplate(0)
    m = Chem.MolFromSmiles('CCOC')
    self.assertTrue(m.HasSubstructMatch(q))
    m = Chem.MolFromSmiles('CO')
    self.assertFalse(m.HasSubstructMatch(q))

    rxn = rdChemReactions.ReactionFromSmarts("[C:1][O:2].[N:3]>>[C:1][N:2]")
    rxn.Initialize()
    rxn.GetReactantTemplate(0).GetAtomWithIdx(0).SetProp('query', 'aliphatic')
    labels = rxn.AddRecursiveQueriesToReaction(qs,'query', getLabels=True)
    self.assertTrue(len(labels), 1)

  def test18GithubIssue16(self):
    rxn = rdChemReactions.ReactionFromSmarts("[F:1]>>[Cl:1]")
    self.assertTrue(rxn)
    rxn.Initialize()
    self.assertRaises(ValueError,lambda : rxn.RunReactants((None,)))

  def test19RemoveUnmappedMoleculesToAgents(self):
    rxn = rdChemReactions.ReactionFromSmarts("[C:1]=[O:2].[N:3].C(=O)O>[OH2].[Na].[Cl]>[N:3]~[C:1]=[O:2]")
    self.failUnless(rxn)
    rxn.Initialize()
    self.failUnless(rxn.GetNumReactantTemplates()==3)
    self.failUnless(rxn.GetNumProductTemplates()==1)
    self.failUnless(rxn.GetNumAgentTemplates()==3)

    rxn.RemoveUnmappedReactantTemplates()
    rxn.RemoveUnmappedProductTemplates()

    self.failUnless(rxn.GetNumReactantTemplates()==2)
    self.failUnless(rxn.GetNumProductTemplates()==1)
    self.failUnless(rxn.GetNumAgentTemplates()==4)

    rxn = rdChemReactions.ReactionFromSmarts("[C:1]=[O:2].[N:3].C(=O)O>>[N:3]~[C:1]=[O:2].[OH2]")
    self.failUnless(rxn)
    rxn.Initialize()
    self.failUnless(rxn.GetNumReactantTemplates()==3)
    self.failUnless(rxn.GetNumProductTemplates()==2)
    self.failUnless(rxn.GetNumAgentTemplates()==0)
    
    agentList=[]
    rxn.RemoveUnmappedReactantTemplates(moveToAgentTemplates=False, targetList=agentList)
    rxn.RemoveUnmappedProductTemplates(targetList=agentList)

    self.failUnless(rxn.GetNumReactantTemplates()==2)
    self.failUnless(rxn.GetNumProductTemplates()==1)
    self.failUnless(rxn.GetNumAgentTemplates()==1)
    self.failUnless(len(agentList)==2)

if __name__ == '__main__':
  unittest.main()
