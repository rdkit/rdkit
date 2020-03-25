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
import tempfile

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

  def testRingSmartsWithTrustedSmiles(self):
    pat = Chem.MolFromSmarts("[C&R1]")
    pat2 = Chem.MolFromSmarts("C@C") # ring bond
    holder = rdSubstructLibrary.CachedTrustedSmilesMolHolder()
    lib = rdSubstructLibrary.SubstructLibrary(holder)
    lib.AddMol(Chem.MolFromSmiles("C1CC1"))

    # make sure we can get an unsanitized molecule that fails (no ring info)
    print("Testing atom rings")
    with self.assertRaises(RuntimeError):
      holder.GetMol(0).HasSubstructMatch(pat)
    print("testing bond rings")
    with self.assertRaises(RuntimeError):
      holder.GetMol(0).HasSubstructMatch(pat2)

    # shouldn't throw
    print("searching atom rings")
    self.assertEqual(len(lib.GetMatches(pat)), 1)
    self.assertEqual(lib.CountMatches(pat), 1)
    print("searching bond rings")
    self.assertEqual(len(lib.GetMatches(pat2)), 1)
    self.assertEqual(lib.CountMatches(pat2), 1)
    print("done")

  def test_init_from_and_to_stream(self):
    mols = makeStereoExamples() * 10
    holder = rdSubstructLibrary.CachedSmilesMolHolder()

    # one day I'll fix this, but we need to write text but read binary
    #  grrr....  something about the python_streambuf handler.
    slib = rdSubstructLibrary.SubstructLibrary(holder, None)

    for mol in mols:
      holder.AddSmiles(Chem.MolToSmiles(mol, isomericSmiles=True))

    if rdSubstructLibrary.SubstructLibraryCanSerialize():
      fd, path = tempfile.mkstemp()
      with open(path, 'w') as file:
        slib.ToStream(file)

      with open(path, 'rb') as file:
        slib2 = rdSubstructLibrary.SubstructLibrary()
        slib2.InitFromStream(file)
        self.assertEqual(len(slib), len(slib2))

    from io import StringIO, BytesIO
    s = StringIO()
    slib.ToStream(s)

    sb = BytesIO(s.getvalue().encode("ascii"))
    self.assertTrue(len(sb.getvalue()) > 0)
    slib3 = rdSubstructLibrary.SubstructLibrary()
    slib3.InitFromStream(sb)
    self.assertEqual(len(slib), len(slib2))


  def test_addpatterns(self):
    pdb_ligands = [
      "CCS(=O)(=O)c1ccc(OC)c(Nc2ncc(-c3cccc(-c4ccccn4)c3)o2)c1",
      "COc1ccc(S(=O)(=O)NCC2CC2)cc1Nc1ncc(-c2cccc(-c3cccnc3)c2)o1",
      "COc1ccc(-c2oc3ncnc(N)c3c2-c2ccc(NC(=O)Nc3cc(C(F)(F)F)ccc3F)cc2)cc1",
      "COC(=O)Nc1nc2ccc(Oc3ccc(NC(=O)Nc4cc(C(F)(F)F)ccc4F)cc3)cc2[nH]1",
      "COc1cc(Nc2ncnc(-c3cccnc3Nc3ccccc3)n2)cc(OC)c1OC",
      "O=C(Nc1ccc(Oc2ccccc2)cc1)c1cccnc1NCc1ccncc1",
      "O=C(Nc1ccc(Oc2ccccc2)cc1)c1cccnc1NCc1ccncc1",
      "CNC(=O)c1cc(Oc2ccc3[nH]c(Nc4ccc(Cl)c(C(F)(F)F)c4)nc3c2)ccn1",
      "CNC(=O)c1cc(Oc2ccc3oc(Nc4ccc(Cl)c(OCC5CCC[NH+]5C)c4)nc3c2)ccn1",
      "CNC(=O)c1cc(Oc2ccc3oc(Nc4ccc(Cl)c(OCC5CCC[NH+]5C)c4)nc3c2)ccn1",
      "COc1cc2nccc(Oc3ccc4c(c3)OCCN4C(=O)Nc3ccc(Cl)cc3)c2cc1OC",
      "CNC(=O)c1c(C)oc2cc(Oc3cc[nH+]c4cc(OCCN5CCOCC5)ccc34)ccc12",
      "COc1cc2[nH+]ccc(Oc3ccc4c(C(=O)Nc5ccc(Cl)cc5)cccc4c3)c2cc1OC",
      "COc1cc2[nH+]ccc(Oc3ccc4c(C(=O)Nc5ccc(Cl)cc5)cccc4c3)c2cc1OC",
      "COc1cc2[nH+]ccc(Oc3ccc4c(C(=O)NC5CC5)cccc4c3)c2cc1OC",
      "COc1cc2[nH+]ccc(Oc3ccc4c(C(=O)NC5CC5)cccc4c3)c2cc1OC",
      "Cc1ccc(C(=O)Nc2cc(CCC[NH+](C)C)cc(C(F)(F)F)c2)cc1Nc1ncccc1-c1ccncn1",
      "COc1cc(Nc2nccc(Nc3ccc4c(C)n[nH]c4c3)n2)cc(OC)c1OC",
      "COc1cc(Nc2nccc(N(C)c3ccc4c(C)n[nH]c4c3)n2)cc(OC)c1OC",
      "Cc1ccn(-c2ccc3c(c2)NCC3(C)C)c(=O)c1-c1ccc2nc(N)ncc2c1",
      "Cc1ccn(-c2ccc3c(c2)NCC3(C)C)c(=O)c1-c1ccc2nc(N)ncc2c1",
      "Cc1ccc(C(=O)NCCC2CCCC2)cc1C(=O)Nc1ccc(N)nc1",
      "Cc1ccc(C(=O)NCCC2CCCC2)cc1C(=O)Nc1ccc(N)nc1",
      "Cc1ccn(-c2cccc(C(F)(F)F)c2)c(=O)c1-c1ccc2nc(N)ncc2c1",
      "Cc1ccn(-c2cccc(C(F)(F)F)c2)c(=O)c1-c1ccc2nc(N)ncc2c1",
      "O=C(Nc1cncnc1)c1c(Cl)ccc2c(Nc3cccc(C(F)(F)F)c3)noc12",
      "O=C(Nc1cncnc1)c1c(Cl)ccc2c(Nc3cccc(C(F)(F)F)c3)noc12",
      "CC1(C)CNc2cc(NC(=O)c3cccnc3NCc3ccncc3)ccc21",
      "CC1(C)CNc2cc(NC(=O)c3cccnc3NCc3ccncc3)ccc21"
    ]
    
    mols = [Chem.MolFromSmiles(smi) for smi in pdb_ligands]
    holder = rdSubstructLibrary.CachedMolHolder()
    patterns = rdSubstructLibrary.PatternHolder()
    
    slib_with_patterns = rdSubstructLibrary.SubstructLibrary(holder, patterns)
    

    for mol in mols:
      slib_with_patterns.AddMol(mol)

    for nthreads in [1, 2, 0]:
      slib_without_patterns = rdSubstructLibrary.SubstructLibrary(holder, None)
      rdSubstructLibrary.AddPatterns(slib_without_patterns, nthreads)
      # check for seg fault
      #  were the fingerprints really created
      slib_without_patterns.GetFpHolder().GetFingerprint(0)
      for mol in mols:
        l1 = slib_with_patterns.CountMatches(mol)
        l2 = slib_without_patterns.CountMatches(mol)
        self.assertTrue(l1)
        self.assertEqual(l1,l2)
                         

  def test_basic_addpatterns(self):
    # add mols
    pdb_ligands = [
      "CCS(=O)(=O)c1ccc(OC)c(Nc2ncc(-c3cccc(-c4ccccn4)c3)o2)c1",
      "COc1ccc(S(=O)(=O)NCC2CC2)cc1Nc1ncc(-c2cccc(-c3cccnc3)c2)o1",
      "COc1ccc(-c2oc3ncnc(N)c3c2-c2ccc(NC(=O)Nc3cc(C(F)(F)F)ccc3F)cc2)cc1",
      "COC(=O)Nc1nc2ccc(Oc3ccc(NC(=O)Nc4cc(C(F)(F)F)ccc4F)cc3)cc2[nH]1",
      "COc1cc(Nc2ncnc(-c3cccnc3Nc3ccccc3)n2)cc(OC)c1OC",
      "O=C(Nc1ccc(Oc2ccccc2)cc1)c1cccnc1NCc1ccncc1",
      "O=C(Nc1ccc(Oc2ccccc2)cc1)c1cccnc1NCc1ccncc1",
      "CNC(=O)c1cc(Oc2ccc3[nH]c(Nc4ccc(Cl)c(C(F)(F)F)c4)nc3c2)ccn1",
      "CNC(=O)c1cc(Oc2ccc3oc(Nc4ccc(Cl)c(OCC5CCC[NH+]5C)c4)nc3c2)ccn1",
      "CNC(=O)c1cc(Oc2ccc3oc(Nc4ccc(Cl)c(OCC5CCC[NH+]5C)c4)nc3c2)ccn1",
      "COc1cc2nccc(Oc3ccc4c(c3)OCCN4C(=O)Nc3ccc(Cl)cc3)c2cc1OC",
      "CNC(=O)c1c(C)oc2cc(Oc3cc[nH+]c4cc(OCCN5CCOCC5)ccc34)ccc12",
      "COc1cc2[nH+]ccc(Oc3ccc4c(C(=O)Nc5ccc(Cl)cc5)cccc4c3)c2cc1OC",
      "COc1cc2[nH+]ccc(Oc3ccc4c(C(=O)Nc5ccc(Cl)cc5)cccc4c3)c2cc1OC",
      "COc1cc2[nH+]ccc(Oc3ccc4c(C(=O)NC5CC5)cccc4c3)c2cc1OC",
      "COc1cc2[nH+]ccc(Oc3ccc4c(C(=O)NC5CC5)cccc4c3)c2cc1OC",
      "Cc1ccc(C(=O)Nc2cc(CCC[NH+](C)C)cc(C(F)(F)F)c2)cc1Nc1ncccc1-c1ccncn1",
      "COc1cc(Nc2nccc(Nc3ccc4c(C)n[nH]c4c3)n2)cc(OC)c1OC",
      "COc1cc(Nc2nccc(N(C)c3ccc4c(C)n[nH]c4c3)n2)cc(OC)c1OC",
      "Cc1ccn(-c2ccc3c(c2)NCC3(C)C)c(=O)c1-c1ccc2nc(N)ncc2c1",
      "Cc1ccn(-c2ccc3c(c2)NCC3(C)C)c(=O)c1-c1ccc2nc(N)ncc2c1",
      "Cc1ccc(C(=O)NCCC2CCCC2)cc1C(=O)Nc1ccc(N)nc1",
      "Cc1ccc(C(=O)NCCC2CCCC2)cc1C(=O)Nc1ccc(N)nc1",
      "Cc1ccn(-c2cccc(C(F)(F)F)c2)c(=O)c1-c1ccc2nc(N)ncc2c1",
      "Cc1ccn(-c2cccc(C(F)(F)F)c2)c(=O)c1-c1ccc2nc(N)ncc2c1",
      "O=C(Nc1cncnc1)c1c(Cl)ccc2c(Nc3cccc(C(F)(F)F)c3)noc12",
      "O=C(Nc1cncnc1)c1c(Cl)ccc2c(Nc3cccc(C(F)(F)F)c3)noc12",
      "CC1(C)CNc2cc(NC(=O)c3cccnc3NCc3ccncc3)ccc21",
      "CC1(C)CNc2cc(NC(=O)c3cccnc3NCc3ccncc3)ccc21"
    ]

    for holder in [rdSubstructLibrary.CachedSmilesMolHolder(),
                   rdSubstructLibrary.CachedTrustedSmilesMolHolder()]:
      for smi in pdb_ligands:
        holder.AddSmiles(smi)

      lib = rdSubstructLibrary.SubstructLibrary(holder)
      rdSubstructLibrary.AddPatterns(lib, numThreads=-1)
      self.assertEqual(len(lib.GetMolHolder()), len(lib.GetFpHolder()))
      for smi in pdb_ligands:
        self.assertTrue( lib.CountMatches(Chem.MolFromSmiles(smi)) )

if __name__ == '__main__':
  unittest.main()
