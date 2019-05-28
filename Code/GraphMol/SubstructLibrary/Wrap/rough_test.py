#  Copyright (C) 2017-2019  Novartis Institute of BioMedical Research
#         All Rights Reserved
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
""" This is a rough coverage test of the python wrapper for the SubstructLibrary

it is intended to be shallow but broad.
"""


import doctest, unittest, os, sys

from rdkit import RDConfig
from rdkit.RDLogger import logger
logger = logger()
from rdkit import Chem
from rdkit.Chem import rdSubstructLibrary
import time
import pickle

def load_tests(loader, tests, ignore):
  tests.addTests(doctest.DocTestSuite(rdSubstructLibrary))
  return tests

def makeStereoExamples():
  el = "NO"
  mols = []
  for e in el:
      for e2 in el:
          if e != e2:
              smi = "C1CCO[C@@](%s)(%s)1"%(e,e2)
              m = Chem.MolFromSmiles(smi)
              if m:
                  mols.append(m)
              smi = "C1CCO[C@](%s)(%s)1"%(e,e2)
              m = Chem.MolFromSmiles(smi)
              if m:
                  mols.append(m)
                  
  return mols
                
class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def test0SubstructLibrary(self):
    for fpholderCls in [None, rdSubstructLibrary.PatternHolder]:
      for holder in [rdSubstructLibrary.MolHolder(), rdSubstructLibrary.CachedMolHolder(),
                     rdSubstructLibrary.CachedSmilesMolHolder()]:
        if fpholderCls: fpholder = fpholderCls()
        else: fpholder = None
        slib_ = rdSubstructLibrary.SubstructLibrary(holder, fpholder)
        for i in range(100):
            m = Chem.MolFromSmiles("c1ccccc1")
            self.assertEqual(slib_.AddMol(m), i)

        libs = [slib_]
        if rdSubstructLibrary.SubstructLibraryCanSerialize():
          serialized1 = pickle.loads(pickle.dumps(slib_))
          serialized2 = rdSubstructLibrary.SubstructLibrary(slib_.Serialize())
          libs.append(serialized1)
          libs.append(serialized2)
          
        for slib in libs:
          res = slib.GetMatches(m)
          t2 = time.time()
          self.assertTrue(len(res) == 100)

          res = slib.GetMatches(m)

          self.assertEqual(len(res), 100)
          self.assertTrue(set(res) == set(list(range(100))))

          res = slib.GetMatches(m, maxResults=100);
          self.assertEqual(len(res), 100)
          self.assertEqual(len(slib.GetMatches(m, startIdx=0, endIdx=100)), 100)

          self.assertTrue(slib.HasMatch(m))
          self.assertEqual(slib.CountMatches(m), 100)

  def test1SubstructLibrary(self):
    for fpholderCls in [None, rdSubstructLibrary.PatternHolder]:
      for holder in [rdSubstructLibrary.MolHolder(), rdSubstructLibrary.CachedMolHolder(),
                     rdSubstructLibrary.CachedSmilesMolHolder()]:
        if fpholderCls: fpholder = fpholderCls()
        else: fpholder = None
        
        slib_ = rdSubstructLibrary.SubstructLibrary(holder, fpholder)
        mols = []
        for i in range(100):
            m = Chem.MolFromSmiles("c1ccccc1")
            self.assertEqual(slib_.AddMol(m), i*2)
            mols.append(m)
            m2 = Chem.MolFromSmiles("CCCC")
            self.assertEqual(slib_.AddMol(m2), i*2+1)
            mols.append(m2)

        libs = [slib_]
        if rdSubstructLibrary.SubstructLibraryCanSerialize():
          serialized1 = pickle.loads(pickle.dumps(slib_))
          serialized2 = rdSubstructLibrary.SubstructLibrary(slib_.Serialize())
          libs.append(serialized1)
          libs.append(serialized2)
          
        for slib in libs:
          res = slib.GetMatches(m)
          self.assertEqual(len(res), 100)
          self.assertEqual(set(res), set(list(range(0,200,2))))

          res = slib.GetMatches(m2)
          self.assertEqual(len(res), 100)
          self.assertTrue(set(res) == set(list(range(1,200,2))))

          res = slib.GetMatches(m)
          self.assertEqual(len(res), 100)

          res = slib.GetMatches(m, maxResults=100);
          self.assertEqual(len(res), 100)

          self.assertEqual(len(slib.GetMatches(m, startIdx=0, endIdx=50*2)), 50)
          self.assertEqual(len(slib.GetMatches(m2, startIdx=1, endIdx=50*2+1)), 50)

          self.assertTrue(slib.HasMatch(m))
          self.assertTrue(slib.HasMatch(m2))
          self.assertEqual(slib.CountMatches(m), 100)
          self.assertEqual(slib.CountMatches(m2), 100)        

  def testOptions(self):
    mols = makeStereoExamples() * 10

    for holderCls in [
        rdSubstructLibrary.MolHolder,
        rdSubstructLibrary.CachedMolHolder,
        rdSubstructLibrary.CachedSmilesMolHolder,
        rdSubstructLibrary.CachedTrustedSmilesMolHolder,
    ]:
      holder = holderCls()
      slib_ = rdSubstructLibrary.SubstructLibrary(holder, None)

      for mol in mols:
        slib_.AddMol(mol)

      libs = [slib_]
      if rdSubstructLibrary.SubstructLibraryCanSerialize():
        serialized1 = pickle.loads(pickle.dumps(slib_))
        serialized2 = rdSubstructLibrary.SubstructLibrary(slib_.Serialize())
        libs.append(serialized1)
        libs.append(serialized2)
          
      for slib in libs:
        core = Chem.MolFromSmarts("C-1-C-C-O-C(-*)(-*)1")          
        res = slib.GetMatches(core)
        self.assertEqual(len(res),
                         len([x for x in mols if x.HasSubstructMatch(core, useChirality=True)]))

        core = Chem.MolFromSmarts("C-1-C-C-O-C(-[O])(-[N])1")
        core.SetProp("core", "core")
        res = slib.GetMatches(core, useChirality=False)
        self.assertEqual(len(res),
                         len([x for x in mols if x.HasSubstructMatch(core, useChirality=False)]))

        core = Chem.MolFromSmarts("C-1-C-C-O-[C@@](-[O])(-[N])1")          
        res = slib.GetMatches(core, useChirality=False)
        self.assertEqual(len(res),
                         len([x for x in mols if x.HasSubstructMatch(core, useChirality=False)]))

        core = Chem.MolFromSmarts("C-1-C-C-O-[C@@](-[O])(-[N])1")          
        res = slib.GetMatches(core)
        self.assertEqual(len(res),
                         len([x for x in mols if x.HasSubstructMatch(core, useChirality=True)]))

  def testSmilesCache(self):
    mols = makeStereoExamples() * 10
    holder = rdSubstructLibrary.CachedSmilesMolHolder()

    slib_ = rdSubstructLibrary.SubstructLibrary(holder, None)

    for mol in mols:
      holder.AddSmiles(Chem.MolToSmiles(mol, isomericSmiles=True))

    libs = [slib_]
    if rdSubstructLibrary.SubstructLibraryCanSerialize():
      serialized1 = pickle.loads(pickle.dumps(slib_))
      serialized2 = rdSubstructLibrary.SubstructLibrary(slib_.Serialize())
      libs.append(serialized1)
      libs.append(serialized2)
          
    for slib in libs:
      core = Chem.MolFromSmarts("C-1-C-C-O-C(-*)(-*)1")          
      res = slib.GetMatches(core)
      self.assertEqual(len(res),
                       len([x for x in mols if x.HasSubstructMatch(core, useChirality=True)]))

      core = Chem.MolFromSmarts("C-1-C-C-O-C(-[O])(-[N])1")
      core.SetProp("core", "core")
      res = slib.GetMatches(core, useChirality=False)
      self.assertEqual(len(res),
                       len([x for x in mols if x.HasSubstructMatch(core, useChirality=False)]))

      core = Chem.MolFromSmarts("C-1-C-C-O-[C@@](-[O])(-[N])1")          
      res = slib.GetMatches(core, useChirality=False)
      self.assertEqual(len(res),
                       len([x for x in mols if x.HasSubstructMatch(core, useChirality=False)]))

      core = Chem.MolFromSmarts("C-1-C-C-O-[C@@](-[O])(-[N])1")          
      res = slib.GetMatches(core)
      self.assertEqual(len(res),
                       len([x for x in mols if x.HasSubstructMatch(core, useChirality=True)]))


  def testTrustedSmilesCache(self):
    mols = makeStereoExamples() * 10
    holder = rdSubstructLibrary.CachedTrustedSmilesMolHolder()

    slib_ = rdSubstructLibrary.SubstructLibrary(holder, None)

    for mol in mols:
      holder.AddSmiles(Chem.MolToSmiles(mol, isomericSmiles=True))

    libs = [slib_]
    if rdSubstructLibrary.SubstructLibraryCanSerialize():
      serialized1 = pickle.loads(pickle.dumps(slib_))
      serialized2 = rdSubstructLibrary.SubstructLibrary(slib_.Serialize())
      libs.append(serialized1)
      libs.append(serialized2)
      
    for slib in libs:
      core = Chem.MolFromSmarts("C-1-C-C-O-C(-*)(-*)1")          
      res = slib.GetMatches(core)
      self.assertEqual(len(res),
                       len([x for x in mols if x.HasSubstructMatch(core, useChirality=True)]))

      core = Chem.MolFromSmarts("C-1-C-C-O-C(-[O])(-[N])1")
      core.SetProp("core", "core")
      res = slib.GetMatches(core, useChirality=False)
      self.assertEqual(len(res),
                       len([x for x in mols if x.HasSubstructMatch(core, useChirality=False)]))

      core = Chem.MolFromSmarts("C-1-C-C-O-[C@@](-[O])(-[N])1")          
      res = slib.GetMatches(core, useChirality=False)
      self.assertEqual(len(res),
                       len([x for x in mols if x.HasSubstructMatch(core, useChirality=False)]))

      core = Chem.MolFromSmarts("C-1-C-C-O-[C@@](-[O])(-[N])1")          
      res = slib.GetMatches(core)
      self.assertEqual(len(res),
                       len([x for x in mols if x.HasSubstructMatch(core, useChirality=True)]))
    
  def testBinaryCache(self):
    mols = makeStereoExamples() * 10
    holder = rdSubstructLibrary.CachedMolHolder()

    slib_ = rdSubstructLibrary.SubstructLibrary(holder, None)

    for mol in mols:
      holder.AddBinary(mol.ToBinary())

    libs = [slib_]
    if rdSubstructLibrary.SubstructLibraryCanSerialize():
      serialized1 = pickle.loads(pickle.dumps(slib_))
      serialized2 = rdSubstructLibrary.SubstructLibrary(slib_.Serialize())
      libs.append(serialized1)
      libs.append(serialized2)
      
    for slib in libs:
      core = Chem.MolFromSmarts("C-1-C-C-O-C(-*)(-*)1")          
      res = slib.GetMatches(core)
      self.assertEqual(len(res),
                       len([x for x in mols if x.HasSubstructMatch(core, useChirality=True)]))

      core = Chem.MolFromSmarts("C-1-C-C-O-C(-[O])(-[N])1")
      core.SetProp("core", "core")
      res = slib.GetMatches(core, useChirality=False)
      self.assertEqual(len(res),
                       len([x for x in mols if x.HasSubstructMatch(core, useChirality=False)]))

      core = Chem.MolFromSmarts("C-1-C-C-O-[C@@](-[O])(-[N])1")          
      res = slib.GetMatches(core, useChirality=False)
      self.assertEqual(len(res),
                       len([x for x in mols if x.HasSubstructMatch(core, useChirality=False)]))

      core = Chem.MolFromSmarts("C-1-C-C-O-[C@@](-[O])(-[N])1")          
      res = slib.GetMatches(core)
      self.assertEqual(len(res),
                       len([x for x in mols if x.HasSubstructMatch(core, useChirality=True)]))
        
if __name__ == '__main__':
  unittest.main()
