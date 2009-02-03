#  $Id$
#
#  Copyright (c) 2007, Novartis Institutes for BioMedical Research Inc.
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
from rdkit import rdBase
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit import Geometry
from rdkit import RDConfig
import unittest
import os,sys
import cPickle as pickle

def feq(v1,v2,tol2=1e-4):
  return abs(v1-v2)<=tol2

def ptEq(pt1, pt2, tol=1e-4):
  return feq(pt1.x,pt2.x,tol) and feq(pt1.y,pt2.y,tol) and feq(pt1.z,pt2.z,tol)

class TestCase(unittest.TestCase) :
  def setUp(self):
    self.dataDir = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','ChemReactions','testData')


  def test1Basics(self):
    rxn = rdChemReactions.ChemicalReaction()
    self.failUnless(rxn.GetNumReactantTemplates()==0)
    self.failUnless(rxn.GetNumProductTemplates()==0)

    r1= Chem.MolFromSmarts('[C:1](=[O:2])O')
    rxn.AddReactantTemplate(r1)
    self.failUnless(rxn.GetNumReactantTemplates()==1)

    r1= Chem.MolFromSmarts('[N:3]')
    rxn.AddReactantTemplate(r1)
    self.failUnless(rxn.GetNumReactantTemplates()==2)

    r1= Chem.MolFromSmarts('[C:1](=[O:2])[N:3]')
    rxn.AddProductTemplate(r1)
    self.failUnless(rxn.GetNumProductTemplates()==1)

    reacts = (Chem.MolFromSmiles('C(=O)O'),Chem.MolFromSmiles('N'))
    ps = rxn.RunReactants(reacts)
    self.failUnless(len(ps)==1)
    self.failUnless(len(ps[0])==1)
    self.failUnless(ps[0][0].GetNumAtoms()==3)
        
    ps = rxn.RunReactants(list(reacts))
    self.failUnless(len(ps)==1)
    self.failUnless(len(ps[0])==1)
    self.failUnless(ps[0][0].GetNumAtoms()==3)    

  def test2DaylightParser(self):
    rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]')
    self.failUnless(rxn)
    self.failUnless(rxn.GetNumReactantTemplates()==2)
    self.failUnless(rxn.GetNumProductTemplates()==1)
    
    reacts = (Chem.MolFromSmiles('C(=O)O'),Chem.MolFromSmiles('N'))
    ps = rxn.RunReactants(reacts)
    self.failUnless(len(ps)==1)
    self.failUnless(len(ps[0])==1)
    self.failUnless(ps[0][0].GetNumAtoms()==3)  
    
    reacts = (Chem.MolFromSmiles('CC(=O)OC'),Chem.MolFromSmiles('CN'))
    ps = rxn.RunReactants(reacts)
    self.failUnless(len(ps)==1)
    self.failUnless(len(ps[0])==1)
    self.failUnless(ps[0][0].GetNumAtoms()==5)  
    
  def test3MDLParsers(self):
    fileN = os.path.join(self.dataDir,'AmideBond.rxn')
    rxn = rdChemReactions.ReactionFromRxnFile(fileN) 
    self.failUnless(rxn)

    self.failUnless(rxn.GetNumReactantTemplates()==2)
    self.failUnless(rxn.GetNumProductTemplates()==1)
    
    reacts = (Chem.MolFromSmiles('C(=O)O'),Chem.MolFromSmiles('N'))
    ps = rxn.RunReactants(reacts)
    self.failUnless(len(ps)==1)
    self.failUnless(len(ps[0])==1)
    self.failUnless(ps[0][0].GetNumAtoms()==3)  
    
    rxnBlock = file(fileN,'r').read()
    rxn = rdChemReactions.ReactionFromRxnBlock(rxnBlock) 
    self.failUnless(rxn)

    self.failUnless(rxn.GetNumReactantTemplates()==2)
    self.failUnless(rxn.GetNumProductTemplates()==1)
    
    reacts = (Chem.MolFromSmiles('C(=O)O'),Chem.MolFromSmiles('N'))
    ps = rxn.RunReactants(reacts)
    self.failUnless(len(ps)==1)
    self.failUnless(len(ps[0])==1)
    self.failUnless(ps[0][0].GetNumAtoms()==3)  

  def test4ErrorHandling(self):
    self.failUnlessRaises(ValueError,lambda x='[C:1](=[O:2])Q.[N:3]>>[C:1](=[O:2])[N:3]':rdChemReactions.ReactionFromSmarts(x))
    self.failUnlessRaises(ValueError,lambda x='[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]Q':rdChemReactions.ReactionFromSmarts(x))
    self.failUnlessRaises(ValueError,lambda x='[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]>>CC':rdChemReactions.ReactionFromSmarts(x))

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
    self.failUnlessRaises(ValueError,lambda x=block:rdChemReactions.ReactionFromRxnBlock(x))
        
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
    #self.failUnlessRaises(ValueError,lambda x=block:rdChemReactions.ReactionFromRxnBlock(x))

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
    #self.failUnlessRaises(ValueError,lambda x=block:rdChemReactions.ReactionFromRxnBlock(x))

  def test5Validation(self):
    rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]')
    self.failUnless(rxn)
    self.failUnless(rxn.Validate()==(0,0))

    rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:1])O.[N:3]>>[C:1](=[O:2])[N:3]')
    self.failUnless(rxn)
    self.failUnless(rxn.Validate()==(0,2))

    rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])[O:4].[N:3]>>[C:1](=[O:2])[N:3]')
    self.failUnless(rxn)
    self.failUnless(rxn.Validate()==(1,0))

    rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3][C:5]')
    self.failUnless(rxn)
    self.failUnless(rxn.Validate()==(0,1))

  def test6Exceptions(self):
    rxn = rdChemReactions.ReactionFromSmarts('[C:1]Cl>>[C:1]')
    self.failUnless(rxn)
    self.failUnlessRaises(ValueError,lambda x=rxn:x.RunReactants(()))
    self.failUnlessRaises(ValueError,lambda x=rxn:x.RunReactants((Chem.MolFromSmiles('CC'),Chem.MolFromSmiles('C'))))
    ps=rxn.RunReactants((Chem.MolFromSmiles('CCCl'),))
    self.failUnless(len(ps)==1)
    self.failUnless(len(ps[0])==1)

  def _test7Leak(self):
    rxn = rdChemReactions.ReactionFromSmarts('[C:1]Cl>>[C:1]')
    self.failUnless(rxn)
    print 'running: '
    for i in range(1e5):
      ps=rxn.RunReactants((Chem.MolFromSmiles('CCCl'),))
      self.failUnless(len(ps)==1)
      self.failUnless(len(ps[0])==1)
      if not i%1000: print i

  def test8Properties(self):
    rxn = rdChemReactions.ReactionFromSmarts('[O:1]>>[O:1][3#0]')
    self.failUnless(rxn)
    ps=rxn.RunReactants((Chem.MolFromSmiles('CO'),))
    self.failUnless(len(ps)==1)
    self.failUnless(len(ps[0])==1)
    Chem.SanitizeMol(ps[0][0])
    self.failUnless(ps[0][0].GetAtomWithIdx(1).GetMass()==3);

        
if __name__ == '__main__':
  unittest.main()
