#  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
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
#       products derived from this software without specific prior written
#       permission.
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
import itertools, time

class TestCase(unittest.TestCase) :
  def setUp(self):
    self.dataDir = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','ChemReactions','testData')

  def testCartesianProduct(self):
    rxn = rdChemReactions.ChemicalReaction();
    rgroups = [[Chem.MolFromSmiles("C")]*10,
               [Chem.MolFromSmiles("N")]*5,
               [Chem.MolFromSmiles("O")]*6]

    cartProd = rdChemReactions.CartesianProductStrategy()
    cartProd.Initialize(rxn, rgroups)
    self.assertEquals(cartProd.GetNumPermutations(), 10*5*6)
    groups = []
    while cartProd:
      groups.append(tuple(cartProd.Next()))
    self.assertEquals(len(groups), 10*5*6)
    # see if we are equal to the Python implementation
    g = list(itertools.product( list(range(10)), list(range(5)), list(range(6)) ))
    self.assertEquals(set(g), set(groups))

  def testRandomSample(self):
    rgroups = [[Chem.MolFromSmiles("C")]*10,
               [Chem.MolFromSmiles("N")]*5,
               [Chem.MolFromSmiles("O")]*6]
    rxn = rdChemReactions.ChemicalReaction();

    randProd = rdChemReactions.RandomSampleStrategy()
    randProd.Initialize(rxn, rgroups)
    self.assertEquals(randProd.GetNumPermutations(), 10*5*6)
    groups = []
    for i in range(10*5*6):
      groups.append(tuple(randProd.Next()))
    print( len(set(groups)), "out of", 10*5*6 )

    randProd = rdChemReactions.RandomSampleStrategy()
    randProd.Initialize(rxn, rgroups)
    self.assertEquals(randProd.GetNumPermutations(), 10*5*6)
    groups = []
    for i in range(10):
      groups.append(tuple(randProd.Next()))

    for i in range(3):
      print( i, len(set([g[i] for g in groups])), "out of", [10,5,6][i] )

  def testRandomSampleAllBBs(self):
    rxn = rdChemReactions.ChemicalReaction();
    rgroups = [[Chem.MolFromSmiles("C")]*10,
               [Chem.MolFromSmiles("N")]*5,
               [Chem.MolFromSmiles("O")]*6]

    randProd = rdChemReactions.RandomSampleAllBBsStrategy()
    randProd.Initialize(rxn, rgroups)
    self.assertEquals(randProd.GetNumPermutations(), 10*5*6)
    groups = []
    for i in range(10*5*6):
      groups.append(tuple(randProd.Next()))

    print( len(set(groups)), "out of", 10*5*6 )

    randProd = rdChemReactions.RandomSampleAllBBsStrategy()
    randProd.Initialize(rxn, rgroups)
    self.assertEquals(randProd.GetNumPermutations(), 10*5*6)
    groups = []
    for i in range(10):
      groups.append(tuple(randProd.Next()))

    for i in range(3):
      print( i, len(set([g[i] for g in groups])), "out of", [10,5,6][i] )
      self.assertEquals(len(set([g[i] for g in groups])), [10,5,6][i])
    
  def testTimings(self):
    rxn = rdChemReactions.ChemicalReaction();

    rgroups = [[Chem.MolFromSmiles("C")]*17000,
               [Chem.MolFromSmiles("N")]*50000,
               [Chem.MolFromSmiles("O")]*4000]
    cartProd = rdChemReactions.CartesianProductStrategy()
    randProd = rdChemReactions.RandomSampleStrategy()
    randAllBBs = rdChemReactions.RandomSampleAllBBsStrategy()
    for r in [cartProd, randProd, randAllBBs]:
      r.Initialize(rxn, rgroups)
      num = 10000000
      t1 = time.time()
      r.Skip(num)
      t2 = time.time()
      print("%s Skipped %s in %s seconds"%(r, num, t2-t1))
            
  def testEnumerateLibrary(self):
    smirks_thiourea = "[N;$(N-[#6]):3]=[C;$(C=S):1].[N;$(N[#6]);!$(N=*);!$([N-]);!$(N#*);!$([ND3]);!$([ND4]);!$(N[O,N]);!$(N[C,S]=[S,O,N]):2]>>[N:3]-[C:1]-[N+0:2]"
    rxn = rdChemReactions.ReactionFromSmarts(smirks_thiourea)
    reagents = [
      [Chem.MolFromSmiles('C=CCN=C=S'), Chem.MolFromSmiles('CC=CCN=C=S')],
      [Chem.MolFromSmiles('NCc1ncc(Cl)cc1Br'),
       Chem.MolFromSmiles('NCCc1ncc(Cl)cc1Br'),
       Chem.MolFromSmiles('NCCCc1ncc(Cl)cc1Br'),
     ]
    ]

    enumerator = rdChemReactions.EnumerateLibrary(rxn, reagents)
    self.assertTrue(enumerator)
    smiresults = ['C=CCNC(=S)NCc1ncc(Cl)cc1Br',
                  'CC=CCNC(=S)NCc1ncc(Cl)cc1Br',
                  'C=CCNC(=S)NCCc1ncc(Cl)cc1Br',
                  'CC=CCNC(=S)NCCc1ncc(Cl)cc1Br',
                  'C=CCNC(=S)NCCCc1ncc(Cl)cc1Br',
                  'CC=CCNC(=S)NCCCc1ncc(Cl)cc1Br']
    results = [Chem.MolToSmiles(Chem.MolFromSmiles(smi)) for smi in smiresults]
    
    pickle = enumerator.Serialize()
    enumerator2 = rdChemReactions.EnumerateLibrary()
    enumerator2.InitFromString(pickle)
    
    print("==", enumerator.GetEnumerator().Type(), enumerator2.GetEnumerator().Type())
    self.assertEquals(enumerator.GetEnumerator().Type(), enumerator2.GetEnumerator().Type())

    out = []
    for en in [enumerator, enumerator2]:
      i = 0
      for i, prods in enumerate(en):
        for mols in prods:
          self.assertEquals(len(mols), 1)
          smi = Chem.MolToSmiles(mols[0])
          if en is enumerator:
            out.append(smi)
          self.assertEquals(smi, results[i])

        if en is enumerator and i == 1:
          # save the state not at the start
          pickle_at_2 = enumerator.Serialize()
      self.assertEquals(i, 5)

    # see if we can restore the enumeration from the middle
    out3 = []
    enumerator3 = rdChemReactions.EnumerateLibrary()
    enumerator3.InitFromString(pickle_at_2)
    for prods in enumerator3:
      for mols in prods:
        self.assertEquals(len(mols), 1)
        smi = Chem.MolToSmiles(mols[0])
        out3.append(smi)

    self.assertEquals(out[2:], out3)
    # test smiles interface
    enumerator = rdChemReactions.EnumerateLibrary(rxn, reagents)
    i = 0
    while enumerator:
      for mols in enumerator.nextSmiles():
        self.assertEquals(len(mols), 1)
        self.assertEquals(mols[0], results[i])
      i += 1
    self.assertEquals(i, 6)
      
  def testRandomEnumerateLibrary(self):
    smirks_thiourea = "[N;$(N-[#6]):3]=[C;$(C=S):1].[N;$(N[#6]);!$(N=*);!$([N-]);!$(N#*);!$([ND3]);!$([ND4]);!$(N[O,N]);!$(N[C,S]=[S,O,N]):2]>>[N:3]-[C:1]-[N+0:2]"
    rxn = rdChemReactions.ReactionFromSmarts(smirks_thiourea)
    reagents = [
      [Chem.MolFromSmiles('C=CCN=C=S'), Chem.MolFromSmiles('CC=CCN=C=S')],
      [Chem.MolFromSmiles('NCc1ncc(Cl)cc1Br'),
       Chem.MolFromSmiles('NCCc1ncc(Cl)cc1Br'),
       Chem.MolFromSmiles('NCCCc1ncc(Cl)cc1Br'),
     ]
    ]

    enumerator = rdChemReactions.EnumerateLibrary(rxn, reagents, rdChemReactions.RandomSampleStrategy())
    self.assertTrue(enumerator)
    smiresults = ['C=CCNC(=S)NCc1ncc(Cl)cc1Br',
                  'CC=CCNC(=S)NCc1ncc(Cl)cc1Br',
                  'C=CCNC(=S)NCCc1ncc(Cl)cc1Br',
                  'CC=CCNC(=S)NCCc1ncc(Cl)cc1Br',
                  'C=CCNC(=S)NCCCc1ncc(Cl)cc1Br',
                  'CC=CCNC(=S)NCCCc1ncc(Cl)cc1Br']
    results = [Chem.MolToSmiles(Chem.MolFromSmiles(smi)) for smi in smiresults]
    
    pickle = enumerator.Serialize()
    enumerator2 = rdChemReactions.EnumerateLibrary()
    enumerator2.InitFromString(pickle)
    
    self.assertEquals(enumerator.GetEnumerator().Type(), enumerator2.GetEnumerator().Type())

    iteren = iter(enumerator)
    iteren2 = iter(enumerator2)

    outsmiles = []
    for i in range(10):
      prods1 = iteren.next()
      prods2 = iteren2.next()
      self.assertEquals(len(prods1), len(prods2))
      for mols1, mols2 in zip(prods1, prods2):
          self.assertEquals(len(mols1), 1)
          smi1 = Chem.MolToSmiles(mols1[0])
          self.assertEquals(smi1, Chem.MolToSmiles(mols2[0]))
          outsmiles.append(smi1)

      if i == 1:
        pickle_at_2 = enumerator.Serialize()

    # make sure we can pickle the state as well
    enumerator3 = rdChemReactions.EnumerateLibrary()
    enumerator3.InitFromString(pickle_at_2)
    iteren3 = iter(enumerator3)
    outsmiles2 = []
    for i in range(8):
      prods3 = iteren3.next()
      for mols3 in prods3:
          self.assertEquals(len(mols3), 1)
          smi1 = Chem.MolToSmiles(mols3[0])
          self.assertEquals(smi1, Chem.MolToSmiles(mols3[0]))
          outsmiles2.append(smi1)

    self.assertEquals(outsmiles2, outsmiles[2:])
    
  def testRandomEnumerateAllBBsLibrary(self):
    smirks_thiourea = "[N;$(N-[#6]):3]=[C;$(C=S):1].[N;$(N[#6]);!$(N=*);!$([N-]);!$(N#*);!$([ND3]);!$([ND4]);!$(N[O,N]);!$(N[C,S]=[S,O,N]):2]>>[N:3]-[C:1]-[N+0:2]"
    rxn = rdChemReactions.ReactionFromSmarts(smirks_thiourea)
    reagents = [
      [Chem.MolFromSmiles('C=CCN=C=S'), Chem.MolFromSmiles('CC=CCN=C=S')],
      [Chem.MolFromSmiles('NCc1ncc(Cl)cc1Br'),
       Chem.MolFromSmiles('NCCc1ncc(Cl)cc1Br'),
       Chem.MolFromSmiles('NCCCc1ncc(Cl)cc1Br'),
     ]
    ]
    enumerator = rdChemReactions.EnumerateLibrary(rxn, reagents, rdChemReactions.RandomSampleAllBBsStrategy())
    self.assertTrue(enumerator)
    smiresults = ['C=CCNC(=S)NCc1ncc(Cl)cc1Br',
                  'CC=CCNC(=S)NCc1ncc(Cl)cc1Br',
                  'C=CCNC(=S)NCCc1ncc(Cl)cc1Br',
                  'CC=CCNC(=S)NCCc1ncc(Cl)cc1Br',
                  'C=CCNC(=S)NCCCc1ncc(Cl)cc1Br',
                  'CC=CCNC(=S)NCCCc1ncc(Cl)cc1Br']
    results = [Chem.MolToSmiles(Chem.MolFromSmiles(smi)) for smi in smiresults]
    
    pickle = enumerator.Serialize()
    enumerator2 = rdChemReactions.EnumerateLibrary()
    enumerator2.InitFromString(pickle)
    
    self.assertEquals(enumerator.GetEnumerator().Type(), enumerator2.GetEnumerator().Type())
    iteren = iter(enumerator)
    iteren2 = iter(enumerator2)

    outsmiles = []
    for i in range(10):
      prods1 = iteren.next()
      prods2 = iteren2.next()
      self.assertEquals(len(prods1), len(prods2))
      for mols1, mols2 in zip(prods1, prods2):
          self.assertEquals(len(mols1), 1)
          smi1 = Chem.MolToSmiles(mols1[0])
          self.assertEquals(smi1, Chem.MolToSmiles(mols2[0]))
          outsmiles.append(smi1)

      if i == 1:
        pickle_at_2 = enumerator.Serialize()

    # make sure we can pickle the state as well
    enumerator3 = rdChemReactions.EnumerateLibrary()
    enumerator3.InitFromString(pickle_at_2)
    self.assertEquals(enumerator.GetEnumerator().Type(), enumerator3.GetEnumerator().Type())

    iteren3 = iter(enumerator3)
    outsmiles2 = []
    for i in range(8):
      prods3 = iteren3.next()
      for mols3 in prods3:
          self.assertEquals(len(mols3), 1)
          smi1 = Chem.MolToSmiles(mols3[0])
          self.assertEquals(smi1, Chem.MolToSmiles(mols3[0]))
          outsmiles2.append(smi1)

    self.assertEquals(outsmiles2, outsmiles[2:])

            
  def testRGroupState(self):
    smirks_thiourea = "[N;$(N-[#6]):3]=[C;$(C=S):1].[N;$(N[#6]);!$(N=*);!$([N-]);!$(N#*);!$([ND3]);!$([ND4]);!$(N[O,N]);!$(N[C,S]=[S,O,N]):2]>>[N:3]-[C:1]-[N+0:2]"
    rxn = rdChemReactions.ReactionFromSmarts(smirks_thiourea)
    reagents = [
      [Chem.MolFromSmiles('C=CCN=C=S'), Chem.MolFromSmiles('CC=CCN=C=S')],
      [Chem.MolFromSmiles('NCc1ncc(Cl)cc1Br'),
       Chem.MolFromSmiles('NCCc1ncc(Cl)cc1Br'),
       Chem.MolFromSmiles('NCCCc1ncc(Cl)cc1Br'),
     ]
    ]

    def tostr(l):
      return [[str(x) for x in v] for v in l]
    enumerator = rdChemReactions.EnumerateLibrary(rxn, reagents)
    state = enumerator.GetState()

    enumerator = rdChemReactions.EnumerateLibrary(rxn, reagents,
                                                  rdChemReactions.RandomSampleStrategy())
    state = enumerator.GetState()
    p = enumerator.nextSmiles()
    p2 = enumerator.nextSmiles()
    enumerator.SetState(state)
    self.assertEquals(tostr(enumerator.nextSmiles()), tostr(p))
    self.assertEquals(tostr(enumerator.nextSmiles()), tostr(p2))
    state = enumerator.GetState()
    
    enumerator = rdChemReactions.EnumerateLibrary(rxn, reagents,
                                                  rdChemReactions.RandomSampleAllBBsStrategy())
    state = enumerator.GetState()
    p = enumerator.nextSmiles()
    p2 = enumerator.nextSmiles()
    enumerator.SetState(state)
    self.assertEquals(tostr(enumerator.nextSmiles()), tostr(p))
    self.assertEquals(tostr(enumerator.nextSmiles()), tostr(p2))
    state = enumerator.GetState()

if __name__ == '__main__':
  unittest.main()
