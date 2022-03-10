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


import unittest
import os,sys, copy

import pickle

from rdkit import rdBase
from rdkit import Chem
from rdkit.Chem import AllChem,rdChemReactions
from rdkit import Geometry
from rdkit import RDConfig
import itertools, time
import numpy as np

def log(s):
  rdBase.LogErrorMsg("== " + s)

class TestCase(unittest.TestCase) :
  def setUp(self):
    self.dataDir = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','ChemReactions','testData')

  def testCartesianProduct(self):
    log("testCartesianProduct")
    rxn = rdChemReactions.ChemicalReaction()
    rgroups = [[Chem.MolFromSmiles("C")]*10,
               [Chem.MolFromSmiles("N")]*5,
               [Chem.MolFromSmiles("O")]*6]

    cartProd = rdChemReactions.CartesianProductStrategy()
    cartProd.Initialize(rxn, rgroups)
    self.assertEquals(cartProd.GetNumPermutations(), 10*5*6)
    groups = []
    count = 0
    print (cartProd.__bool__())
    while cartProd:
      groups.append(tuple(cartProd.next()))
#      count += 1
#      assert count <= cartProd.GetNumPermutations()
    self.assertEquals(len(groups), 10*5*6)
    # see if we are equal to the Python implementation
    g = list(itertools.product( list(range(10)), list(range(5)), list(range(6)) ))
    self.assertEquals(set(g), set(groups))
    copy.copy(cartProd)
    
  def testRandomSample(self):
    log("testRandomSample")
    rgroups = [[Chem.MolFromSmiles("C")]*10,
               [Chem.MolFromSmiles("N")]*5,
               [Chem.MolFromSmiles("O")]*6]
    rxn = rdChemReactions.ChemicalReaction()

    randProd = rdChemReactions.RandomSampleStrategy()
    randProd.Initialize(rxn, rgroups)
    self.assertEquals(randProd.GetNumPermutations(), 10*5*6)
    groups = []
    for i in range(10*5*6):
      groups.append(tuple(randProd.next()))
    print( len(set(groups)), "out of", 10*5*6 )

    randProd = rdChemReactions.RandomSampleStrategy()
    randProd.Initialize(rxn, rgroups)
    self.assertEquals(randProd.GetNumPermutations(), 10*5*6)
    groups = []
    for i in range(10):
      groups.append(tuple(randProd.next()))

    for i in range(3):
      print( i, len(set([g[i] for g in groups])), "out of", [10,5,6][i] )
    copy.copy(randProd)
    
  def testRandomSampleAllBBs(self):
    log("testRandomSampleAllBBs")
    rxn = rdChemReactions.ChemicalReaction()
    rgroups = [[Chem.MolFromSmiles("C")]*10,
               [Chem.MolFromSmiles("N")]*5,
               [Chem.MolFromSmiles("O")]*6]

    randProd = rdChemReactions.RandomSampleAllBBsStrategy()
    randProd.Initialize(rxn, rgroups)
    self.assertEquals(randProd.GetNumPermutations(), 10*5*6)
    groups = []
    for i in range(10*5*6):
      groups.append(tuple(randProd.next()))

    print( len(set(groups)), "out of", 10*5*6 )

    randProd = rdChemReactions.RandomSampleAllBBsStrategy()
    randProd.Initialize(rxn, rgroups)
    self.assertEquals(randProd.GetNumPermutations(), 10*5*6)
    groups = []
    for i in range(10):
      groups.append(tuple(randProd.next()))

    for i in range(3):
      print( i, len(set([g[i] for g in groups])), "out of", [10,5,6][i] )
      self.assertEquals(len(set([g[i] for g in groups])), [10,5,6][i])
    copy.copy(randProd)
    
  def testTimings(self):
    log("testTimings")
    rxn = rdChemReactions.ChemicalReaction()

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

  def testEvenPairsSampling(self):
    rxn = rdChemReactions.ChemicalReaction()
    
    rgroups = [[Chem.MolFromSmiles("C")]*10,
               [Chem.MolFromSmiles("N")]*10,
               [Chem.MolFromSmiles("O")]*10]

    rxn = rdChemReactions.ChemicalReaction()
    count = 0
    pairs01 = {}
    pairs12 = {}
    pairs02 = {}
    
    strategy = rdChemReactions.EvenSamplePairsStrategy()
    strategy.Initialize(rxn, rgroups)
    # try 100 samples
    while count < 100:
      v = strategy.next()
      p01 = (v[0], v[1])
      p12 = (v[1], v[2])
      p02 = (v[0], v[2])
      pairs01[p01] = pairs01.get(p01, 0) + 1
      pairs12[p01] = pairs12.get(p12, 0) + 1
      pairs02[p01] = pairs02.get(p02, 0) + 1
      count += 1

    # each pair should be used roughly once
    self.assertEquals(np.median(list(pairs01.values())), 1.0)
    self.assertEquals(np.median(list(pairs02.values())), 1.0)
    self.assertEquals(np.median(list(pairs12.values())), 1.0)

    # now try 1000
    pairs01 = {}
    pairs12 = {}
    pairs02 = {}
    strategy = rdChemReactions.EvenSamplePairsStrategy()
    strategy.Initialize(rxn, rgroups)
    count = 0
    while count < 1000:
      v = strategy.next()
      p01 = (v[0], v[1])
      p12 = (v[1], v[2])
      p02 = (v[0], v[2])
      pairs01[p01] = pairs01.get(p01, 0) + 1
      pairs12[p01] = pairs12.get(p12, 0) + 1
      pairs02[p01] = pairs02.get(p02, 0) + 1
      count += 1

    # each pair should be used roughly 10 times
    self.assertTrue( 9 <= np.median(list(pairs01.values())) <= 11)
    self.assertTrue( 9 <= np.median(list(pairs02.values())) <= 11)
    self.assertTrue( 9 <= np.median(list(pairs12.values())) <= 11)

    # now try 500
    pairs01 = {}
    pairs12 = {}
    pairs02 = {}
    strategy = rdChemReactions.EvenSamplePairsStrategy()
    strategy.Initialize(rxn, rgroups)
    count = 0
    while count < 500:
      v = strategy.next()
      p01 = (v[0], v[1])
      p12 = (v[1], v[2])
      p02 = (v[0], v[2])
      pairs01[p01] = pairs01.get(p01, 0) + 1
      pairs12[p01] = pairs12.get(p12, 0) + 1
      pairs02[p01] = pairs02.get(p02, 0) + 1
      count += 1

    # each pair should be used roughly 5 times      
    self.assertTrue( 4 <= np.median(list(pairs01.values())) <= 6)
    self.assertTrue( 4 <= np.median(list(pairs02.values())) <= 6)
    self.assertTrue( 4 <= np.median(list(pairs12.values())) <= 6)
    
    
    self.assertTrue("PAIRSTAT" in strategy.Stats())

  def testEnumerateLibrary(self):
    log("testEnumerateLibrary")
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

    # need to initialize the reaction before getting the binary serialization
    rxn.Initialize()
    self.assertEquals(rxn.ToBinary(), enumerator.GetReaction().ToBinary())

    bbs = enumerator.GetReagents()
    for i in range(len(bbs)):
      for j in range(len(bbs[i])):
        self.assertTrue(Chem.MolToSmiles(reagents[i][j]) == Chem.MolToSmiles(bbs[i][j]))
    
    smiresults = ['C=CCNC(=S)NCc1ncc(Cl)cc1Br',
                  'CC=CCNC(=S)NCc1ncc(Cl)cc1Br',
                  'C=CCNC(=S)NCCc1ncc(Cl)cc1Br',
                  'CC=CCNC(=S)NCCc1ncc(Cl)cc1Br',
                  'C=CCNC(=S)NCCCc1ncc(Cl)cc1Br',
                  'CC=CCNC(=S)NCCCc1ncc(Cl)cc1Br']
    results = [Chem.MolToSmiles(Chem.MolFromSmiles(smi)) for smi in smiresults]

    enumerators = [enumerator]

    # add serialized enumerators as well for testing if possible
    if rdChemReactions.EnumerateLibraryCanSerialize():
      pickle = enumerator.Serialize()
      enumerator2 = rdChemReactions.EnumerateLibrary()
      enumerator2.InitFromString(pickle)
      
      # make sure old pickles work
      enumerator3 = rdChemReactions.EnumerateLibrary()
      enumerator3.InitFromString(open(os.path.join(self.dataDir, "enumeration.pickle"), 'rb').read())
    
      print("==", enumerator.GetEnumerator().Type(), enumerator2.GetEnumerator().Type())
      self.assertEquals(enumerator.GetEnumerator().Type(), enumerator2.GetEnumerator().Type())
      enumerators.append(enumerator2)
      enumerators.append(enumerator3)

    # check for fully sampled and deterministic ordering in final index values
    expected_positions = [[0, 0],[1, 0],[0, 1],[1, 1],[0, 2],[1, 2]]
    
    out = []
    for en in enumerators:
      i = 0
      positions = []
      for i, prods in enumerate(en):
        positions.append( list(en.GetPosition()) )
        for mols in prods:
          self.assertEquals(len(mols), 1)
          smi = Chem.MolToSmiles(mols[0])
          if en is enumerator:
            out.append(smi)
          self.assertEquals(smi, results[i])

        if en is enumerator and i == 1 and rdChemReactions.EnumerateLibraryCanSerialize():
          # save the state not at the start
          pickle_at_2 = enumerator.Serialize()
      self.assertEquals(i, 5)
      self.assertEquals(positions, expected_positions)
      
    if rdChemReactions.EnumerateLibraryCanSerialize():      
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
    log("testRandomEnumerateLibrary")
    smirks_thiourea = "[N;$(N-[#6]):3]=[C;$(C=S):1].[N;$(N[#6]);!$(N=*);!$([N-]);!$(N#*);!$([ND3]);!$([ND4]);!$(N[O,N]);!$(N[C,S]=[S,O,N]):2]>>[N:3]-[C:1]-[N+0:2]"
    rxn = rdChemReactions.ReactionFromSmarts(smirks_thiourea)
    reagents = [
      [Chem.MolFromSmiles('C=CCN=C=S'), Chem.MolFromSmiles('CC=CCN=C=S')],
      [Chem.MolFromSmiles('NCc1ncc(Cl)cc1Br'),
       Chem.MolFromSmiles('NCCc1ncc(Cl)cc1Br'),
       Chem.MolFromSmiles('NCCCc1ncc(Cl)cc1Br'),
     ]
    ]

    enumerator = rdChemReactions.EnumerateLibrary(rxn, reagents,
                                                  rdChemReactions.RandomSampleStrategy())
    self.assertTrue(enumerator)
    smiresults = ['C=CCNC(=S)NCc1ncc(Cl)cc1Br',
                  'CC=CCNC(=S)NCc1ncc(Cl)cc1Br',
                  'C=CCNC(=S)NCCc1ncc(Cl)cc1Br',
                  'CC=CCNC(=S)NCCc1ncc(Cl)cc1Br',
                  'C=CCNC(=S)NCCCc1ncc(Cl)cc1Br',
                  'CC=CCNC(=S)NCCCc1ncc(Cl)cc1Br']
    results = [Chem.MolToSmiles(Chem.MolFromSmiles(smi)) for smi in smiresults]

    enumerator = rdChemReactions.EnumerateLibrary(rxn, reagents,
                                                  rdChemReactions.RandomSampleStrategy())
    iteren = iter(enumerator)
    res = set()
    count = 0
    while res != set(results):
      count += 1
      if count > 100000:
        print("Unable to find enumerate set with 100,000 random samples!", file=sys.stderr)
        self.assertEquals(res,set(results))

      prod = iteren.next()
      for mols in prod:
        smi1 = Chem.MolToSmiles(mols[0])
        res.add(smi1)
        
    if rdChemReactions.EnumerateLibraryCanSerialize():
      enumerator = rdChemReactions.EnumerateLibrary(rxn, reagents,
                                                    rdChemReactions.RandomSampleStrategy())
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
    log("testRandomEnumerateAllBBsLibrary")
    smirks_thiourea = "[N;$(N-[#6]):3]=[C;$(C=S):1].[N;$(N[#6]);!$(N=*);!$([N-]);!$(N#*);!$([ND3]);!$([ND4]);!$(N[O,N]);!$(N[C,S]=[S,O,N]):2]>>[N:3]-[C:1]-[N+0:2]"
    rxn = rdChemReactions.ReactionFromSmarts(smirks_thiourea)
    reagents = [
      [Chem.MolFromSmiles('C=CCN=C=S'), Chem.MolFromSmiles('CC=CCN=C=S')],
      [Chem.MolFromSmiles('NCc1ncc(Cl)cc1Br'),
       Chem.MolFromSmiles('NCCc1ncc(Cl)cc1Br'),
       Chem.MolFromSmiles('NCCCc1ncc(Cl)cc1Br'),
     ]
    ]
    enumerator = rdChemReactions.EnumerateLibrary(rxn, reagents,
                                                  rdChemReactions.RandomSampleAllBBsStrategy())
    self.assertTrue(enumerator)

    # test the BB sampling here
    strategy = iter(enumerator)
    r1 = set()
    r2 = set()
    strategy.next()    
    groups = strategy.GetPosition()
    print("**", list(groups), file=sys.stderr)
    r1.add(groups[0])
    r2.add(groups[1])
    strategy.next()
    groups = strategy.GetPosition()
    print("**", list(groups),file=sys.stderr)
    r1.add(groups[0])
    r2.add(groups[1])
    self.assertEquals(r1, set([0,1])) # two bbs at reagent one all sampled at one iteration
    strategy.next()
    groups = strategy.GetPosition()
    print("**", list(groups),file=sys.stderr)
    r1.add(groups[0])
    r2.add(groups[1])
    self.assertEquals(r2, set([0,1,2])) # three bbs at reagent one all sampled in three iterations
    
    smiresults = ['C=CCNC(=S)NCc1ncc(Cl)cc1Br',
                  'CC=CCNC(=S)NCc1ncc(Cl)cc1Br',
                  'C=CCNC(=S)NCCc1ncc(Cl)cc1Br',
                  'CC=CCNC(=S)NCCc1ncc(Cl)cc1Br',
                  'C=CCNC(=S)NCCCc1ncc(Cl)cc1Br',
                  'CC=CCNC(=S)NCCCc1ncc(Cl)cc1Br']
    results = [Chem.MolToSmiles(Chem.MolFromSmiles(smi)) for smi in smiresults]


    if rdChemReactions.EnumerateLibraryCanSerialize():
      enumerator = rdChemReactions.EnumerateLibrary(rxn, reagents,
                                                    rdChemReactions.RandomSampleAllBBsStrategy())
      self.assertTrue(enumerator)
      
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
    if not rdChemReactions.EnumerateLibraryCanSerialize():
      print("-- Skipping testRGroupState, serialization of EnumerateLibrary not enabled", file=sys.stderr)
      return
    
    log("testRGroupState")
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
    p = enumerator.nextSmiles()
    p2 = enumerator.nextSmiles()
    enumerator.SetState(state)
    self.assertEquals(tostr(enumerator.nextSmiles()), tostr(p))
    self.assertEquals(tostr(enumerator.nextSmiles()), tostr(p2))

    enumerator = rdChemReactions.EnumerateLibrary(rxn, reagents,
                                                  rdChemReactions.RandomSampleStrategy())
    
    state = enumerator.GetState()
    p = enumerator.nextSmiles()
    p2 = enumerator.nextSmiles()
    enumerator.SetState(state)
    self.assertEquals(tostr(enumerator.nextSmiles()), tostr(p))
    self.assertEquals(tostr(enumerator.nextSmiles()), tostr(p2))

    enumerator = rdChemReactions.EnumerateLibrary(rxn, reagents,
                                                  rdChemReactions.RandomSampleAllBBsStrategy())
    state = enumerator.GetState()
    p = enumerator.nextSmiles()
    p2 = enumerator.nextSmiles()
    enumerator.SetState(state)
    self.assertEquals(tostr(enumerator.nextSmiles()), tostr(p))
    self.assertEquals(tostr(enumerator.nextSmiles()), tostr(p2))


    enumerator = rdChemReactions.EnumerateLibrary(rxn, reagents)
    smiresults = ['C=CCNC(=S)NCc1ncc(Cl)cc1Br',
                  'CC=CCNC(=S)NCc1ncc(Cl)cc1Br',
                  'C=CCNC(=S)NCCc1ncc(Cl)cc1Br',
                  'CC=CCNC(=S)NCCc1ncc(Cl)cc1Br',
                  'C=CCNC(=S)NCCCc1ncc(Cl)cc1Br',
                  'CC=CCNC(=S)NCCCc1ncc(Cl)cc1Br']
    smiresults = [Chem.MolToSmiles(Chem.MolFromSmiles(smi)) for smi in smiresults]
    enumerator.GetEnumerator().Skip(10)
    enumerator.ResetState()

    results  = []
    for result in enumerator:
      for prodSet in result:
        for mol in prodSet:
          results.append( Chem.MolToSmiles(mol) )

    self.assertEquals(results, smiresults)

  def testRemovingBadMatches(self):
    log("testRemoveBadMatches")
    smirks_thiourea = "[N;$(N-[#6]):3]=[C;$(C=S):1].[N;$(N[#6]);!$(N=*);!$([N-]);!$(N#*);!$([ND3]);!$([ND4]);!$(N[O,N]);!$(N[C,S]=[S,O,N]):2]>>[N:3]-[C:1]-[N+0:2]"
    
    rxn = rdChemReactions.ReactionFromSmarts(smirks_thiourea)
    # invert matches so nothing matches
    reagents = [
      [Chem.MolFromSmiles('NCc1ncc(Cl)cc1Br'),
       Chem.MolFromSmiles('NCCc1ncc(Cl)cc1Br'),
       Chem.MolFromSmiles('NCCCc1ncc(Cl)cc1Br'),
     ],

      [Chem.MolFromSmiles('C=CCN=C=S'),
       Chem.MolFromSmiles('CC=CCN=C=S'),
       Chem.MolFromSmiles('CCC'),
       Chem.MolFromSmiles('CCCCC'),
     ],
    ]

    enumerator = rdChemReactions.EnumerateLibrary(rxn, reagents)
    self.assertEquals([], list(enumerator))

  def testRemoveInsaneReagents(self):
    rxndata = "$RXN\nUntitled Document-1\n  ChemDraw10291618492D\n\n  3  1\n$MOL\n\n\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.4125    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0\n   -0.4125    0.0000    0.0000 R2  0  0  0  0  0  0  0  0  0  2  0  0\n  1  2  1  0        0\nM  END\n$MOL\n\n\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n   -0.4125    0.0000    0.0000 R1  0  0  0  0  0  0  0  0  0  1  0  0\n    0.4125    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0        0\nM  END\n$MOL\n\n\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.4125    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  5  0  0\n   -0.4125    0.0000    0.0000 R4  0  0  0  0  0  0  0  0  0  4  0  0\n  1  2  1  0        0\nM  END\n$MOL\n\n\n\n 14 15  0  0  0  0  0  0  0  0999 V2000\n    0.5072   -0.5166    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.5072    0.3084    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2949   -0.7616    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7817   -0.0880    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2967    0.5794    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5558   -1.5443    0.0000 R1  0  0  0  0  0  0  0  0  0  1  0  0\n   -0.2073    0.7208    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9218    0.3083    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9217   -0.5167    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.2073   -0.9292    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.6362    0.7208    0.0000 N   0  0  0  0  0  0  0  0  0  3  0  0\n    1.5452    1.3661    0.0000 N   0  0  0  0  0  0  0  0  0  5  0  0\n    2.3507    1.5443    0.0000 R4  0  0  0  0  0  0  0  0  0  4  0  0\n   -2.3507    0.3083    0.0000 R2  0  0  0  0  0  0  0  0  0  2  0  0\n  1  2  2  0        0\n  1  3  1  0        0\n  3  4  1  0        0\n  4  5  1  0        0\n  5  2  1  0        0\n  3  6  1  0        0\n  2  7  1  0        0\n  7  8  2  0        0\n  8  9  1  0        0\n  9 10  2  0        0\n 10  1  1  0        0\n  8 11  1  0        0\n 12 13  1  0        0\n 11 14  1  0        0\n 12  5  1  0        0\nM  END\n"

    rxn = AllChem.ReactionFromRxnBlock(rxndata)
    bbs = []
    r1 = [ Chem.MolFromSmiles("CCNCC"),
           Chem.MolFromSmiles("NCC"),
           ]
    r2 = [ Chem.MolFromSmiles("ClC1CCCC1"),
           Chem.MolFromSmiles("ClC1CCCC1Cl"),
           ]
    r3 = [ Chem.MolFromSmiles("CCNCC"),
           Chem.MolFromSmiles("NCC"),
           ]
    bbs = [r1, r2, r3]

    # nothing matches!
    for i,reagent in enumerate(rxn.GetReactants()):
      for bb in bbs[i]:
        self.assertFalse(bb.HasSubstructMatch(reagent))

    # everything matches - yay sanitization!
    rdChemReactions.SanitizeRxn(rxn)
    for i,reagent in enumerate(rxn.GetReactants()):
      for bb in bbs[i]:
        self.assertTrue(bb.HasSubstructMatch(reagent))

    en = rdChemReactions.EnumerateLibrary(rxn, bbs)
    self.assertTrue(len(en.GetReagents()[0]) == 2)
    self.assertTrue(len(en.GetReagents()[1]) == 2)
    self.assertTrue(len(en.GetReagents()[2]) == 2)

    #####################################################################################
    # Match only at rgroups (ChemDraw style)
    rxn = AllChem.ReactionFromRxnBlock(rxndata)
    expected_matches = [[False,True], [True,True],[False, True] ]
    rdChemReactions.SanitizeRxn(rxn, params=rdChemReactions.GetChemDrawRxnAdjustParams())
    for i,(reagent, expected) in enumerate(zip(rxn.GetReactants(), expected_matches)):
      match = [bb.HasSubstructMatch(reagent) for reagent in bbs[i]]
      self.assertTrue(match, expected)

    # Now try EnumerateLibrary
    en = rdChemReactions.EnumerateLibrary(rxn, bbs)
    self.assertTrue(len(en.GetReagents()[0]) == 1)
    self.assertTrue(len(en.GetReagents()[1]) == 2)
    self.assertTrue(len(en.GetReagents()[2]) == 1)


    #####################################################################################
    # now set the removal options ot only make one product per reagent set
    rxn = AllChem.ReactionFromRxnBlock(rxndata)      
    rdChemReactions.SanitizeRxn(rxn)

    opts = rdChemReactions.EnumerationParams()
    opts.reagentMaxMatchCount = 1
    en = rdChemReactions.EnumerateLibrary(rxn, bbs, params=opts)
    self.assertTrue(len(en.GetReagents()[0]) == 1)
    self.assertTrue(len(en.GetReagents()[1]) == 1)
    self.assertTrue(len(en.GetReagents()[2]) == 1)

    #####################################################################################
    # now set the removal options ot only make one product per reagent set
    #  but wt
    rxn = AllChem.ReactionFromRxnBlock(rxndata)
    rdChemReactions.SanitizeRxn(rxn,
                                params=rdChemReactions.GetChemDrawRxnAdjustParams())      


    opts = rdChemReactions.EnumerationParams()
    opts.reagentMaxMatchCount = 1
    en = rdChemReactions.EnumerateLibrary(rxn, bbs, params=opts)
    self.assertTrue(len(en.GetReagents()[0]) == 1)
    self.assertTrue(len(en.GetReagents()[1]) == 1)
    self.assertTrue(len(en.GetReagents()[2]) == 1)
    

if __name__ == '__main__':
  unittest.main()
