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
import sys,re

# These are the definitions that will be applied to fragment molecules:
environs = {
  'L1':('[#0:1,#6,#7,#8][C:2](=[O:3])','[*:1][C:2](=[O:3])-[1*]'),
  'L2':('[N:21;!$(N=*)]-;!@[#0:22,#6]','[2*]-[N:21][*:22]'),
  # this one turned out to be too tricky to define above, so we set it off
  # in its own definition:
  'L2a':('[N:21;D3;R;$(N(@[C;!$(C=*)])@[C;!$(C=*)])]','[2*]-[N:21]'),
  'L3':('[O:31]-;!@[#0:32,#6]','[3*]-[O:31][*:32]'),
  'L4':('[#6:41]-;!@[C:42;!$(C=*)]','[*:41]-[C:42]-[4*]'),
  'L5':('[N:51;!$(N*!-*);!$(N=*);!$(N-[!C;!#0])]-[#0:52,C]','[5*]-[N:51][*:52]'),
  'L6':('[C:61](=[O:62])-;!@[#0:63,#6,#7,#8]','[6*]-[C:61](=[O:62])[*:63]'), 
  'L7a':('[#6:71]-;!@[C:72]','[#6:71]-;!@[C:72]=[7*]'),
  'L7b':('[C:74]-;!@[#6:73]','[#6:73]-;!@[C:74]=[7*]'),
  'L8':('[C:81]-;!@[C:82]-;!@[#6:83]','[8*]-[C:81]-;!@[C:82]-;!@[*:83]'),
  'L9':('[n:91;$(n(:[c,n,o,s]):[c,n,o,s])]','[9*]-[*:91]'),
  'L10':('[N:101;R;$(N(@C(=O))@[C,N,O,S])]','[10*]-[N:101]'),
  'L11':('[S:112](-;!@[#0:111,#6])','[*:111][S:112]-[11*]'),
  'L12':('[#6:121,#0][S:122](=[O:123])(=[O:124])','[*:121][S:122](=[O:123])(=[O:124])-[12*]'),
  'L13':('[C:131;$(C(-;@[C,N,O,S])-;@[N,O,S])]','[*:131]-[13*]'),
  'L14':('[c:141;$(c(:[c,n,o,s]):[n,o,s])]','[*:141]-[14*]'),
  'L14b':('[c:142;$(c(:[c,n,o,s]):[n,o,s])]','[*:142]-[14*]'),
  'L15':('[C:151;$(C(-;@C)-;@C)]','[*:151]-[15*]'),
  'L16':('[c:161;$(c(:c):c)]','[*:161]-[16*]'),
  'L16b':('[c:162;$(c(:c):c)]','[*:162]-[16*]'),

  }
reactionDefs = (
  # L1
  [
  "{L1}-;!@{L3}>>{L1p}.{L3p}",
  "{L1}-;!@{L2}>>{L1p}.{L2p}",
  "{L1}-;!@{L2a}>>{L1p}.{L2ap}",
  "{L1}-;!@{L10}>>{L1p}.{L10p}",
  ],

  # L2 
  [
  "{L12}-;!@{L2}>>{L12p}.{L2p}",
  "{L14}-;!@{L2}>>{L14p}.{L2p}",
  "{L16}-;!@{L2}>>{L16p}.{L2p}",
  ],
  
  # L2a 
  [
  "{L12}-;!@{L2a}>>{L12p}.{L2ap}",
  "{L14}-;!@{L2a}>>{L14p}.{L2ap}",
  "{L16}-;!@{L2a}>>{L16p}.{L2ap}",
  ],
  
  # L3 
  [
  "{L4}-;!@{L3}>>{L4p}.{L3p}",
  "{L13}-;!@{L3}>>{L13p}.{L3p}",
  "{L14}-;!@{L3}>>{L14p}.{L3p}",
  "{L15}-;!@{L3}>>{L15p}.{L3p}",
  "{L16}-;!@{L3}>>{L16p}.{L3p}",
  ],
  
  # L4
  [
  "{L4}-;!@{L5}>>{L4p}.{L5p}",
  "{L4}-;!@{L3}>>{L4p}.{L3p}",
  "{L4}-;!@{L11}>>{L4p}.{L11p}",
  ],

  # L5
  [
  "{L13}-;!@{L5}>>{L13p}.{L5p}",
  "{L15}-;!@{L5}>>{L15p}.{L5p}",
  ],
  
  # L6
  [
  "{L13}-;!@{L6}>>{L13p}.{L6p}",
  "{L14}-;!@{L6}>>{L14p}.{L6p}",
  "{L15}-;!@{L6}>>{L15p}.{L6p}",
  "{L16}-;!@{L6}>>{L16p}.{L6p}",
  ],
  
  # L7
  [
  "{L7a}=;!@{L7b}>>{L7ap}.{L7bp}",
  ],

  # L8
  [
  "{L9}-;!@{L8}>>{L9p}.{L8p}",
  "{L10}-;!@{L8}>>{L10p}.{L8p}",
  "{L13}-;!@{L8}>>{L13p}.{L8p}",
  "{L14}-;!@{L8}>>{L14p}.{L8p}",
  "{L15}-;!@{L8}>>{L15p}.{L8p}",
  "{L16}-;!@{L8}>>{L16p}.{L8p}",
  ],
  
  # L9
  [
  "{L9}-;!@{L15}>>{L9p}.{L15p}",
  "{L9}-;!@{L16}>>{L9p}.{L16p}",
  ],
  
  # L10
  [
  "{L10}-;!@{L13}>>{L10p}.{L13p}",
  "{L10}-;!@{L14}>>{L10p}.{L14p}",
  "{L10}-;!@{L15}>>{L10p}.{L15p}",
  "{L10}-;!@{L16}>>{L10p}.{L16p}",
  ],
  
  # L11
  [
  "{L11}-;!@{L13}>>{L11p}.{L13p}",
  "{L11}-;!@{L14}>>{L11p}.{L14p}",
  "{L11}-;!@{L15}>>{L11p}.{L15p}",
  "{L11}-;!@{L16}>>{L11p}.{L16p}",
  ],
  
  # L12
  # none left

  # L13
  [
  "{L13}-;!@{L14}>>{L13p}.{L14p}",
  "{L13}-;!@{L15}>>{L13p}.{L15p}",
  "{L13}-;!@{L16}>>{L13p}.{L16p}",
  ],
  
  # L14
  [
  "{L14}-;!@{L14b}>>{L14p}.{L14bp}", # not in original paper
  "{L14}-;!@{L15}>>{L14p}.{L15p}",
  "{L14}-;!@{L16}>>{L14p}.{L16p}",
  ],
  
  # L15
  [
  "{L15}-;!@{L16}>>{L15p}.{L16p}",
  ],

  # L16
  [
  "{L16}-;!@{L16b}>>{L16p}.{L16bp}",  # not in original paper
  ],
  )

smartsGps=reactionDefs
for k,v in environs.iteritems():
  for gp in smartsGps:
    for j,sma in enumerate(gp):
      sma = sma.replace('{%sp}'%k,v[1]).replace('{%s}'%k,v[0])
      gp[j] =sma
for gp in smartsGps:
  for defn in gp:
    try:
      Reactions.ReactionFromSmarts(defn)
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
          print reactionDefs[gpIdx][rxnIdx]
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
  return doctest.testmod(sys.modules["__main__"])

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

