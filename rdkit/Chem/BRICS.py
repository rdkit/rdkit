# $Id$
#
#  Copyright (c) 2009, Novartis Institutes for BioMedical Research Inc.
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
# Created by Greg Landrum, Nov 2008

""" Implementation of the BRICS algorithm from Degen et al. ChemMedChem *3* 1503-7 (2008)

"""
from rdkit import Chem
from rdkit.Chem import rdChemReactions as Reactions
import sys,re,random

# These are the definitions that will be applied to fragment molecules:
environs = {
  'L1':'C([#0,#6,#7,#8])(=[O])',
  'L2':'[N;!R;!$(N=*)]-;!@[#0,#6]',
  # this one turned out to be too tricky to define above, so we set it off
  # in its own definition:
  'L2a':'[N;D3;R;$(N(@[C;!$(C=*)])@[C;!$(C=*)])]',
  'L3':'O-;!@[#0,#6]',
  'L4':'[C;!$(C=*)]-;!@[#6]',
  'L5':'[N;!$(N*!-*);!$(N=*);!$(N-[!C;!#0])]-[#0,C]',
  'L6':'C(=[O])-;!@[#0,#6,#7,#8]',
  'L7a':'C-;!@[#6]',
  'L7b':'C-;!@[#6]',
  'L8':'C-;!@C-;!@[#6]',
  'L9':'[n;$(n(:[c,n,o,s]):[c,n,o,s])]',
  'L10':'[N;R;$(N(@C(=O))@[C,N,O,S])]',
  'L11':'S(-;!@[#0,#6])',
  'L12':'S([#6,#0])(=O)(=O)',
  'L13':'[C;$(C(-;@[C,N,O,S])-;@[N,O,S])]',
  'L14':'[c;$(c(:[c,n,o,s]):[n,o,s])]',
  'L14b':'[c;$(c(:[c,n,o,s]):[n,o,s])]',
  'L15':'[C;$(C(-;@C)-;@C)]',
  'L16':'[c;$(c(:c):c)]',
  'L16b':'[c;$(c(:c):c)]',
  }
reactionDefs = (
  # L1
  [
    ('1','3','-'),
    ('1','2','-'),
    ('1','2a','-'),
    ('1','10','-'),
  ],

   # L2 
   [
    ('2','12','-'),
    ('2','14','-'),
    ('2','16','-'),
   ],
  
   # L2a 
   [
    ('2a','12','-'),
    ('2a','14','-'),
    ('2a','16','-'),
   ],
  
   # L3 
   [
    ('3','4','-'),
    ('3','13','-'),
    ('3','14','-'),
    ('3','15','-'),
    ('3','16','-'),
   ],
  
   # L4
   [
    ('4','5','-'),
    ('4','11','-'),
   ],

   # L5
   [
    ('5','13','-'),
    ('5','15','-'),
   ],
  
    # L6
    [
    ('6','13','-'),
    ('6','14','-'),
    ('6','15','-'),
    ('6','16','-'),
    ],
  
    # L7
    [
    ('7a','7b','='),
    ],

    # L8
    [
    ('8','9','-'),
    ('8','10','-'),
    ('8','13','-'),
    ('8','14','-'),
    ('8','15','-'),
    ('8','16','-'),
    ],
  
    # L9
    [
    ('9','15','-'),
    ('9','16','-'),
    ],
  
    # L10
    [
    ('10','13','-'),
    ('10','14','-'),
    ('10','15','-'),
    ('10','16','-'),
    ],
  
    # L11
    [
    ('11','13','-'),
    ('11','14','-'),
    ('11','15','-'),
    ('11','16','-'),
    ],
  
    # L12
    # none left

    # L13
    [
    ('13','14','-'),
    ('13','15','-'),
    ('13','16','-'),
    ],
  
    # L14
    [
    ('14','14','-'),# not in original paper
    ('14','15','-'),
    ('14','16','-'),
    ],
  
    # L15
    [
    ('15','16','-'),
    ],

    # L16
    [
    ('16','16','-'), # not in original paper
    ],
  )
import copy
smartsGps=copy.deepcopy(reactionDefs)
for gp in smartsGps:
  for j,defn in enumerate(gp):
    g1,g2,bnd = defn
    r1=environs['L'+g1]
    r2=environs['L'+g2]
    g1 = re.sub('[a-z,A-Z]','',g1)
    g2 = re.sub('[a-z,A-Z]','',g2)
    sma='[$(%s):1]%s;!@[$(%s):2]>>[%s*]-[*:1].[%s*]-[*:2]'%(r1,bnd,r2,g1,g2)
    gp[j] =sma

for gp in smartsGps:
  for defn in gp:
    try:
      t=Reactions.ReactionFromSmarts(defn)
      t.Initialize()
    except:
      print defn
      raise

reactions = tuple([[Reactions.ReactionFromSmarts(y) for y in x] for x in smartsGps])
reverseReactions = []
for i,rxnSet in enumerate(smartsGps):
  for j,sma in enumerate(rxnSet):
    rs,ps = sma.split('>>')
    sma = '%s>>%s'%(ps,rs)
    rxn = Reactions.ReactionFromSmarts(sma)
    labels = re.findall(r'\[([0-9]+?)\*\]',ps)
    rxn._matchers=[Chem.MolFromSmiles('[%s*]'%x) for x in labels]
    reverseReactions.append(rxn)

def FindBRICSBonds(mol,randomizeOrder=False,silent=True):
  """ returns the bonds in a molecule that BRICS would cleave

  >>> from rdkit import Chem
  >>> m = Chem.MolFromSmiles('CCCOCC')
  >>> res = list(FindBRICSBonds(m))
  >>> res
  [((3, 2), ('3', '4')), ((3, 4), ('3', '4'))]

  a more complicated case:
  >>> m = Chem.MolFromSmiles('CCCOCCC(=O)c1ccccc1')
  >>> res = list(FindBRICSBonds(m))
  >>> res
  [((3, 2), ('3', '4')), ((3, 4), ('3', '4')), ((6, 8), ('6', '16'))]

  we can also randomize the order of the results:
  >>> random.seed(23)
  >>> res = list(FindBRICSBonds(m,randomizeOrder=True))
  >>> res
  [((6, 8), ('6', '16')), ((3, 2), ('3', '4')), ((3, 4), ('3', '4'))]

  Note that this is a generator function :
  >>> res = FindBRICSBonds(m)
  >>> res
  <generator object at 0x...>
  >>> res.next()
  ((3, 2), ('3', '4'))

  >>> m = Chem.MolFromSmiles('CC=CC')
  >>> res = list(FindBRICSBonds(m))
  >>> res.sort()
  >>> res
  [((1, 2), ('7', '7'))]
  
  """
  letter = re.compile('[a-z,A-Z]')
  indices = range(len(reactionDefs))
  bondsDone=set()
  if randomizeOrder: random.shuffle(indices)
  for gpIdx in indices:
    if randomizeOrder:
      compats =reactionDefs[gpIdx][:]
      random.shuffle(compats)
    else:
      compats =reactionDefs[gpIdx]
    for i1,i2,bType in compats:
      e1 = environs['L%s'%i1]
      e2 = environs['L%s'%i2]
      patt = '[$(%s)]%s[$(%s)]'%(e1,bType,e2)
      patt = Chem.MolFromSmarts(patt)
      matches = mol.GetSubstructMatches(patt)
      i1 = letter.sub('',i1)
      i2 = letter.sub('',i2)
      for match in matches:
        if match not in bondsDone:
          bondsDone.add(match)
          yield(((match[0],match[1]),(i1,i2)))

def BreakBRICSBonds(mol,bonds=None,sanitize=True,silent=True):
  """ breaks the BRICS bonds in a molecule and returns the results

  >>> from rdkit import Chem
  >>> m = Chem.MolFromSmiles('CCCOCC')
  >>> m2=BreakBRICSBonds(m)
  >>> Chem.MolToSmiles(m2,True)
  '[4*]CC.[4*]CCC.[3*]O[3*]'

  a more complicated case:
  >>> m = Chem.MolFromSmiles('CCCOCCC(=O)c1ccccc1')
  >>> m2=BreakBRICSBonds(m)
  >>> Chem.MolToSmiles(m2,True)
  '[16*]c1ccccc1.[3*]O[3*].[4*]CCC.[6*]C(=O)CC[4*]'

  can also specify a limited set of bonds to work with:
  >>> m = Chem.MolFromSmiles('CCCOCC')
  >>> m2 = BreakBRICSBonds(m,[((3, 2), ('3', '4'))])
  >>> Chem.MolToSmiles(m2,True)
  '[3*]OCC.[4*]CCC'
  
  this can be used as an alternate approach for doing a BRICS decomposition by
  following BreakBRICSBonds with a call to Chem.GetMolFrags:
  >>> m = Chem.MolFromSmiles('CCCOCC')
  >>> m2=BreakBRICSBonds(m)
  >>> frags = Chem.GetMolFrags(m2,asMols=True)
  >>> [Chem.MolToSmiles(x,True) for x in frags]
  ['[4*]CCC', '[3*]O[3*]', '[4*]CC']

  """
  if not bonds:
    bonds = FindBRICSBonds(mol)
  if not bonds:
    return Chem.Mol(mol.ToBinary())
  eMol = Chem.EditableMol(mol)
  nAts = mol.GetNumAtoms()
  for indices,dummyTypes in bonds:
    ia,ib = indices
    obond = mol.GetBondBetweenAtoms(ia,ib)
    bondType=obond.GetBondType()
    eMol.RemoveBond(ia,ib)

    da,db = dummyTypes
    atoma = Chem.Atom(0)
    atoma.SetMass(int(da))
    atoma.SetNoImplicit(True)
    idxa = nAts
    nAts+=1
    eMol.AddAtom(atoma)
    eMol.AddBond(ia,idxa,bondType)
    
    atomb = Chem.Atom(0)
    atomb.SetMass(int(db))
    atomb.SetNoImplicit(True)
    idxb = nAts
    nAts+=1
    eMol.AddAtom(atomb)
    eMol.AddBond(ib,idxb,bondType)

  res = eMol.GetMol()
  if sanitize:
    Chem.SanitizeMol(res)
  return res

def BRICSDecompose(mol,allNodes=None,minFragmentSize=1,onlyUseReactions=None,
                   silent=True,keepNonLeafNodes=False,singlePass=False):
  """ returns the BRICS decomposition for a molecule

  >>> from rdkit import Chem
  >>> m = Chem.MolFromSmiles('CCCOCc1cc(c2ncccc2)ccc1')
  >>> res = list(BRICSDecompose(m))
  >>> res.sort()
  >>> res
  ['[14*]c1ncccc1', '[16*]c1cccc(C[4*])c1', '[3*]O[3*]', '[4*]CCC']

  nexavar, an example from the paper (corrected):
  >>> m = Chem.MolFromSmiles('CNC(=O)C1=NC=CC(OC2=CC=C(NC(=O)NC3=CC(=C(Cl)C=C3)C(F)(F)F)C=C2)=C1')
  >>> res = list(BRICSDecompose(m))
  >>> res.sort()
  >>> res
  ['[1*]C([1*])=O', '[1*]C([6*])=O', '[14*]c1nccc([16*])c1', '[16*]c1ccc(Cl)c(C(F)(F)F)c1', '[16*]c1ccc([16*])cc1', '[2*]NC', '[2*]N[2*]', '[3*]O[3*]']

  """
  global reactions
  mSmi = Chem.MolToSmiles(mol,1)
  
  if allNodes is None:
    allNodes=set()

  if mSmi in allNodes:
    return set()

  activePool={mSmi:mol}
  allNodes.add(mSmi)
  for gpIdx,reactionGp in enumerate(reactions):
    newPool = {}
    while activePool:
      matched=False
      nSmi = activePool.keys()[0]
      mol = activePool.pop(nSmi)
      for rxnIdx,reaction in enumerate(reactionGp):
        if onlyUseReactions and (gpIdx,rxnIdx) not in onlyUseReactions:
          continue
        if not silent:
          print '--------'
          print smartsGps[gpIdx][rxnIdx]
        ps = reaction.RunReactants((mol,))
        if ps:
          if not silent: print  nSmi,'->',len(ps),'products'
          for prodSeq in ps:
            seqOk=True
            # we want to disqualify small fragments, so sort the product sequence by size
            prodSeq = [(prod.GetNumAtoms(onlyHeavy=True),prod) for prod in prodSeq]
            prodSeq.sort()
            for nats,prod in prodSeq:
              pSmi = Chem.MolToSmiles(prod,1)
              if minFragmentSize>0:
                nDummies = pSmi.count('*')
                if nats-nDummies<minFragmentSize:
                  seqOk=False
                  break
              prod.pSmi = pSmi
            if seqOk:
              matched=True
              for nats,prod in prodSeq:
                pSmi = prod.pSmi
                #print '\t',nats,pSmi
                if pSmi not in allNodes:
                  if not singlePass:
                    activePool[pSmi] = prod
                  allNodes.add(pSmi)
      if singlePass or not matched:
        newPool[nSmi]=mol
    activePool = newPool
  if not singlePass:
    res = set(activePool.keys())
  else:
    res = allNodes
  return res


import random
dummyPattern=Chem.MolFromSmiles('[*]')
def BRICSBuild(fragments,onlyCompleteMols=True,seeds=None,uniquify=True,
               scrambleReagents=True,maxDepth=3):
  seen = set()
  if not seeds:
    seeds = list(fragments)
  if scrambleReagents:
    seeds = list(seeds)
    random.shuffle(seeds)
  if scrambleReagents:
    tempReactions = list(reverseReactions)
    random.shuffle(tempReactions)
  else:
    tempReactions=reverseReactions
  for seed in seeds:
    seedIsR1=False
    seedIsR2=False
    nextSteps=[]
    for rxn in tempReactions:
      if seed.HasSubstructMatch(rxn._matchers[0]):
        seedIsR1=True
      if seed.HasSubstructMatch(rxn._matchers[1]):
        seedIsR2=True
      for fragment in fragments:
        ps = None
        if fragment.HasSubstructMatch(rxn._matchers[0]):
          if seedIsR2:
            ps = rxn.RunReactants((fragment,seed))
        if fragment.HasSubstructMatch(rxn._matchers[1]):
          if seedIsR1:
            ps = rxn.RunReactants((seed,fragment))
        if ps:
          for p in ps:
            if uniquify:
              pSmi =Chem.MolToSmiles(p[0],True)
              if pSmi in seen:
                continue
              else:
                seen.add(pSmi)
            if p[0].HasSubstructMatch(dummyPattern):
              nextSteps.append(p[0])
              if not onlyCompleteMols:
                yield p[0]
            else:
              yield p[0]
    if nextSteps and maxDepth>0:
      for p in BRICSBuild(fragments,onlyCompleteMols=onlyCompleteMols,
                          seeds=nextSteps,uniquify=uniquify,
                          maxDepth=maxDepth-1):
        if uniquify:
          pSmi =Chem.MolToSmiles(p,True)
          if pSmi in seen:
            continue
          else:
            seen.add(pSmi)
        yield p
  


# ------- ------- ------- ------- ------- ------- ------- -------
# Begin testing code


#------------------------------------
#
#  doctest boilerplate
#
def _test():
  import doctest,sys
  return doctest.testmod(sys.modules["__main__"],
                         optionflags=doctest.ELLIPSIS+doctest.NORMALIZE_WHITESPACE)

if __name__=='__main__':
  import unittest
  class TestCase(unittest.TestCase):
    def test1(self):
      m = Chem.MolFromSmiles('CC(=O)OC')
      res = BRICSDecompose(m)
      self.failUnless(res)
      self.failUnless(len(res)==2)

      m = Chem.MolFromSmiles('CC(=O)N1CCC1=O')
      res = BRICSDecompose(m)
      self.failUnless(res)
      self.failUnless(len(res)==2,res)

      m = Chem.MolFromSmiles('c1ccccc1N(C)C')
      res = BRICSDecompose(m)
      self.failUnless(res)
      self.failUnless(len(res)==2,res)

      m = Chem.MolFromSmiles('c1cccnc1N(C)C')
      res = BRICSDecompose(m)
      self.failUnless(res)
      self.failUnless(len(res)==2,res)

      m = Chem.MolFromSmiles('o1ccnc1N(C)C')
      res = BRICSDecompose(m)
      self.failUnless(res)
      self.failUnless(len(res)==2)

      m = Chem.MolFromSmiles('c1ccccc1OC')
      res = BRICSDecompose(m)
      self.failUnless(res)
      self.failUnless(len(res)==2)

      m = Chem.MolFromSmiles('o1ccnc1OC')
      res = BRICSDecompose(m)
      self.failUnless(res)
      self.failUnless(len(res)==2)

      m = Chem.MolFromSmiles('O1CCNC1OC')
      res = BRICSDecompose(m)
      self.failUnless(res)
      self.failUnless(len(res)==2)

      m = Chem.MolFromSmiles('CCCSCC')
      res = BRICSDecompose(m)
      self.failUnless(res)
      self.failUnless(len(res)==3,res)
      self.failUnless('[11*]S[11*]' in res,res)

      m = Chem.MolFromSmiles('CCNC(=O)C1CC1')
      res = BRICSDecompose(m)
      self.failUnless(res)
      self.failUnless(len(res)==4,res)
      self.failUnless('[2*]N[5*]' in res,res)

    def test2(self):
      # example from the paper, nexavar: 
      m = Chem.MolFromSmiles('CNC(=O)C1=NC=CC(OC2=CC=C(NC(=O)NC3=CC(=C(Cl)C=C3)C(F)(F)F)C=C2)=C1')
      res = BRICSDecompose(m)
      self.failUnless(res)
      self.failUnless(len(res)==8,res)

    def test3(self):
      m = Chem.MolFromSmiles('FC(F)(F)C1=C(Cl)C=CC(NC(=O)NC2=CC=CC=C2)=C1')
      res = BRICSDecompose(m)
      self.failUnless(res)
      self.failUnless(len(res)==4,res)
      self.failUnless('[2*]N[2*]' in res,res)
      self.failUnless('[16*]c1ccccc1' in res,res)


    def test4(self):
      allNodes = set()
      m = Chem.MolFromSmiles('c1ccccc1OCCC')
      res = BRICSDecompose(m,allNodes=allNodes)
      self.failUnless(res)
      leaves=res
      self.failUnless(len(leaves)==3,leaves)
      self.failUnless(len(allNodes)==6,allNodes)
      res = BRICSDecompose(m,allNodes=allNodes)
      self.failIf(res)
      self.failUnless(len(allNodes)==6,allNodes)

      m = Chem.MolFromSmiles('c1ccccc1OCCCC')
      res = BRICSDecompose(m,allNodes=allNodes)
      self.failUnless(res)
      leaves.update(res)
      self.failUnless(len(allNodes)==9,allNodes)
      self.failUnless(len(leaves)==4,leaves)
      
      
      m = Chem.MolFromSmiles('c1cc(C(=O)NCC)ccc1OCCC')
      res = BRICSDecompose(m,allNodes=allNodes)
      self.failUnless(res)
      leaves.update(res)
      self.failUnless(len(leaves)==8,leaves)
      self.failUnless(len(allNodes)==18,allNodes)

    def test5(self):
      allNodes = set()
      frags = [
        '[14*]c1ncncn1',
        '[16*]c1ccccc1',
        '[14*]c1ncccc1',
        ]
      frags = [Chem.MolFromSmiles(x) for x in frags]
      res = BRICSBuild(frags)
      self.failUnless(res)
      res= list(res)
      self.failUnless(len(res)==6)
      smis = [Chem.MolToSmiles(x,True) for x in res]
      self.failUnless('c1ccc(-c2ccccc2)cc1' in smis)
      self.failUnless('c1ccc(-c2ncccc2)cc1' in smis)

    def test5a(self):
      allNodes = set()
      frags = [
        '[3*]O[3*]',
        '[16*]c1ccccc1',
        ]
      frags = [Chem.MolFromSmiles(x) for x in frags]
      res = BRICSBuild(frags)
      self.failUnless(res)
      res=list(res)
      smis = [Chem.MolToSmiles(x,True) for x in res]
      self.failUnless(len(smis)==2,smis)
      self.failUnless('c1ccc(Oc2ccccc2)cc1' in smis)
      self.failUnless('c1ccc(-c2ccccc2)cc1' in smis)

      
      
    def test6(self):
      allNodes = set()
      frags = [
        '[16*]c1ccccc1',
        '[3*]OC',
        '[9*]n1cccc1',
        ]
      frags = [Chem.MolFromSmiles(x) for x in frags]
      res = BRICSBuild(frags)
      self.failUnless(res)
      res= list(res)
      self.failUnless(len(res)==3)
      smis = [Chem.MolToSmiles(x,True) for x in res]
      self.failUnless('c1ccc(-c2ccccc2)cc1' in smis)
      self.failUnless('COc1ccccc1' in smis)
      self.failUnless('c1ccc(-n2cccc2)cc1' in smis)

    def test7(self):
      allNodes = set()
      frags = [
        '[16*]c1ccccc1',
        '[3*]OC',
        '[3*]OCC(=O)[6*]',
        ]
      frags = [Chem.MolFromSmiles(x) for x in frags]
      res = BRICSBuild(frags)
      self.failUnless(res)
      res= list(res)
      smis = [Chem.MolToSmiles(x,True) for x in res]
      self.failUnless(len(res)==3)
      self.failUnless('c1ccc(-c2ccccc2)cc1' in smis)
      self.failUnless('COc1ccccc1' in smis)
      self.failUnless('O=C(COc1ccccc1)c1ccccc1' in smis)

    def test8(self):
      random.seed(23)
      base = Chem.MolFromSmiles("n1cncnc1OCC(C1CC1)OC1CNC1")
      catalog = BRICSDecompose(base)
      self.failUnless(len(catalog)==4,catalog)
      catalog = [Chem.MolFromSmiles(x) for x in catalog]
      ms = [Chem.MolToSmiles(x) for x in BRICSBuild(catalog)]
      self.failUnless(len(ms)==9,ms)

      ts = ['n1cnc(C2CNC2)nc1','n1cnc(-c2ncncn2)nc1','C(OC1CNC1)C(C1CC1)OC1CNC1',
            'n1cnc(OC(COC2CNC2)C2CC2)nc1']
      ts = [Chem.MolToSmiles(Chem.MolFromSmiles(x),True) for x in ts]
      for t in ts:
        self.failUnless(t in ms,(t,ms))
        
    def test9(self):
      m = Chem.MolFromSmiles('CCOc1ccccc1c1ncc(c2nc(NCCCC)ncn2)cc1')
      res=BRICSDecompose(m)
      self.failUnlessEqual(len(res),7)
      self.failUnless('[3*]O[3*]' in res)
      self.failIf('[14*]c1ncnc(NCCCC)n1' in res)
      res = BRICSDecompose(m,singlePass=True)
      self.failUnlessEqual(len(res),11)
      self.failUnless('[3*]OCC' in res)
      self.failUnless('[14*]c1ncnc(NCCCC)n1' in res)

    def test10(self):
      m = Chem.MolFromSmiles('C1CCCCN1c1ccccc1')
      res=BRICSDecompose(m)
      self.failUnlessEqual(len(res),2,res)


  failed,tried = _test()
  if failed:
    sys.exit(failed)

  unittest.main()

