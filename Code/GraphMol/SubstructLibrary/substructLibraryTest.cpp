//  Copyright (c) 2017-2019, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

// std bits
#include <RDGeneral/test.h>
#include <iostream>

// RD bits
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SubstructLibrary/SubstructLibrary.h>
#include <GraphMol/SubstructLibrary/PatternFactory.h>

#include <GraphMol/Substruct/SubstructMatch.h>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>

using namespace RDKit;

namespace {
boost::dynamic_bitset<> runTest(SubstructLibrary &ssslib, const ROMol &pattern,
                                int nThreads) {
  std::vector<unsigned int> libMatches = ssslib.getMatches(pattern, nThreads);
  boost::dynamic_bitset<> hasMatch(ssslib.size());
  for (auto idx : libMatches) {
    hasMatch[idx] = 1;
  }

  for (unsigned int i = 0; i < ssslib.size(); ++i) {
    MatchVectType match;
    bool matched = SubstructMatch(*ssslib.getMol(i), pattern, match);
    // std::cerr << MolToSmiles(*ssslib.getMol(i), true) << " " << hasMatch[i]
    //           << " " << matched << std::endl;
    TEST_ASSERT(hasMatch[i] == matched);
  }
  return hasMatch;
};

void runTest(SubstructLibrary &ssslib, const ROMol &pattern, int nThreads,
             const boost::dynamic_bitset<> &hasMatch) {
  std::vector<unsigned int> libMatches = ssslib.getMatches(pattern, nThreads);
  boost::dynamic_bitset<> hasMatch2(ssslib.size());
  for (auto idx : libMatches) {
    hasMatch2[idx] = 1;
  }
  TEST_ASSERT(hasMatch == hasMatch2);

  for (unsigned int i = 0; i < ssslib.size(); ++i) {
    MatchVectType match;
    bool matched = SubstructMatch(*ssslib.getMol(i), pattern, match);
    // std::cerr << MolToSmiles(*ssslib.getMol(i), true) << " " << hasMatch[i]
    //           << " " << matched << std::endl;
    TEST_ASSERT(hasMatch[i] == matched);
  }
};

}  // namespace

void test1() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test1" << std::endl;

  std::string fName = getenv("RDBASE");
  fName += "/Data/NCI/first_200.props.sdf";
  SDMolSupplier suppl(fName);
  SubstructLibrary ssslib;
  while (!suppl.atEnd()) {
    ROMol *mol = nullptr;
    try {
      mol = suppl.next();
    } catch (...) {
      continue;
    }
    if (!mol) {
      continue;
    }
    ssslib.addMol(*mol);
    delete mol;
  }

  std::vector<SubstructLibrary *> libs;
  libs.push_back(&ssslib);

#ifdef RDK_USE_BOOST_SERIALIZATION
  std::string pickle = ssslib.Serialize();
  SubstructLibrary serialized;
  serialized.initFromString(pickle);
  TEST_ASSERT(serialized.size() == ssslib.size());
  libs.push_back(&serialized);
#endif

  boost::dynamic_bitset<> hasMatch;

  int i = 0;
  for (auto lib : libs) {
    ROMol *query = SmartsToMol("[#6;$([#6]([#6])[!#6])]");
    if (i == 0) {
      hasMatch = runTest(*lib, *query, 1);
    } else {
      runTest(*lib, *query, 1, hasMatch);
    }

#ifdef RDK_TEST_MULTITHREADED
    runTest(*lib, *query, -1, hasMatch);
#endif
    delete query;
    ++i;
  }

  i = 0;
  for (auto lib : libs) {
    ROMol *query = SmartsToMol("[$([O,S]-[!$(*=O)])]");
    if (i == 0) {
      hasMatch = runTest(*lib, *query, 1);
    } else {
      runTest(*lib, *query, 1, hasMatch);
    }

#ifdef RDK_TEST_MULTITHREADED
    runTest(*lib, *query, -1, hasMatch);
#endif
    delete query;
    ++i;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void test2() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test2" << std::endl;

  std::string fName = getenv("RDBASE");
  fName += "/Data/NCI/first_200.props.sdf";
  SDMolSupplier suppl(fName);
  auto *mols = new MolHolder();
  auto *fps = new PatternHolder();
  boost::shared_ptr<MolHolder> mols_ptr(mols);
  boost::shared_ptr<PatternHolder> fps_ptr(fps);

  SubstructLibrary ssslib(mols_ptr, fps_ptr);
  while (!suppl.atEnd()) {
    ROMol *mol = nullptr;
    try {
      mol = suppl.next();
    } catch (...) {
      continue;
    }
    if (!mol) {
      continue;
    }
    ssslib.addMol(*mol);
    delete mol;
  }

  std::vector<SubstructLibrary *> libs;
  libs.push_back(&ssslib);

#ifdef RDK_USE_BOOST_SERIALIZATION
  std::string pickle = ssslib.Serialize();
  SubstructLibrary serialized;
  serialized.initFromString(pickle);
  TEST_ASSERT(serialized.size() == ssslib.size());

  // check to see if we are still the right base type
  MolHolderBase *_holder = serialized.getMolHolder().get();
  TEST_ASSERT(_holder != nullptr);
  TEST_ASSERT(dynamic_cast<MolHolder *>(_holder) != nullptr);
  try {
    serialized.getFingerprints();
  } catch (...) {
    TEST_ASSERT(0);
  }

  libs.push_back(&serialized);
#endif

  for (auto lib : libs) {
    ROMol *query = SmartsToMol("[#6]([#6])[!#6]");
    runTest(*lib, *query, 1);
#ifdef RDK_TEST_MULTITHREADED
    runTest(*lib, *query, -1);
#endif
    delete query;
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void test3() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test3 (stereo options)" << std::endl;

  SubstructLibrary ssslib(boost::make_shared<MolHolder>());
  for (int i = 0; i < 10; ++i) {
    ROMol *m1 = SmilesToMol("C1CCO[C@@](N)(O)1");
    ROMol *m2 = SmilesToMol("C1CCO[C@](N)(O)1");
    ROMol *m3 = SmilesToMol("C1CCO[C@@](O)(N)1");
    ROMol *m4 = SmilesToMol("C1CCO[C@](O)(N)1");
    ssslib.addMol(*m1);
    ssslib.addMol(*m2);
    ssslib.addMol(*m3);
    ssslib.addMol(*m4);
    delete m1;
    delete m2;
    delete m3;
    delete m4;
  }

  std::vector<SubstructLibrary *> libs;
  libs.push_back(&ssslib);

#ifdef RDK_USE_BOOST_SERIALIZATION
  std::string pickle = ssslib.Serialize();
  SubstructLibrary serialized;
  serialized.initFromString(pickle);
  TEST_ASSERT(serialized.size() == ssslib.size());
  libs.push_back(&serialized);
  // check to see if we are still the right base type
  MolHolderBase *_holder = serialized.getMolHolder().get();
  TEST_ASSERT(_holder != nullptr);
  TEST_ASSERT(dynamic_cast<MolHolder *>(_holder) != nullptr);
#endif

  for (auto lib : libs) {
    ROMol *query = SmartsToMol("C-1-C-C-O-C(-[O])(-[N])1");
    std::vector<unsigned int> res = lib->getMatches(*query, true, false);
    TEST_ASSERT(res.size() == 40);

    delete query;
    query = SmartsToMol("C-1-C-C-O-[C@@](-[O])(-[N])1");

    res = lib->getMatches(*query, true, true);
    TEST_ASSERT(res.size() == 20);

    res = lib->getMatches(*query, true, false);
    TEST_ASSERT(res.size() == 40);

    delete query;
  }
  BOOST_LOG(rdErrorLog) << "    Done (stereo options)" << std::endl;
}

void test4() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test4 (trusted smiles)" << std::endl;

  boost::shared_ptr<CachedSmilesMolHolder> holder =
      boost::make_shared<CachedSmilesMolHolder>();
  SubstructLibrary ssslib(holder);

  for (int i = 0; i < 10; ++i) {
    holder->addSmiles("C1CCO[C@@](N)(O)1");
    holder->addSmiles("C1CCO[C@](N)(O)1");
    holder->addSmiles("C1CCO[C@@](O)(N)1");
    holder->addSmiles("C1CCO[C@](O)(N)1");
  }

  std::vector<SubstructLibrary *> libs;
  libs.push_back(&ssslib);

#ifdef RDK_USE_BOOST_SERIALIZATION
  std::string pickle = ssslib.Serialize();
  SubstructLibrary serialized;
  serialized.initFromString(pickle);
  TEST_ASSERT(serialized.size() == ssslib.size());
  libs.push_back(&serialized);
  // check to see if we are still the right base type
  MolHolderBase *_holder = serialized.getMolHolder().get();
  TEST_ASSERT(_holder != nullptr);
  TEST_ASSERT(dynamic_cast<CachedSmilesMolHolder *>(_holder) != nullptr);
#endif

  for (auto lib : libs) {
    ROMol *query = SmartsToMol("C-1-C-C-O-C(-[O])(-[N])1");

    std::vector<unsigned int> res = lib->getMatches(*query, true, false);
    TEST_ASSERT(res.size() == 40);

    delete query;
    query = SmartsToMol("C-1-C-C-O-[C@@](-[O])(-[N])1");

    res = lib->getMatches(*query, true, true);
    TEST_ASSERT(res.size() == 20);

    res = lib->getMatches(*query, true, false);
    TEST_ASSERT(res.size() == 40);
    delete query;
  }

  BOOST_LOG(rdErrorLog) << "    Done (trusted smiles)" << std::endl;
}

/// Tests the code in the docs
//   to make sure it compiles.
void docTest() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing C++ docs" << std::endl;

  ROMol *q = SmartsToMol("C-1-C-C-O-C(-[O])(-[N])1");
  ROMol *m = SmilesToMol("C1CCO[C@@](N)(O)1");
  ROMol &query = *q;
  ROMol &mol = *m;

  {
    SubstructLibrary lib;
    lib.addMol(mol);
    std::vector<unsigned int> results = lib.getMatches(query);
    for (std::vector<unsigned int>::const_iterator matchIndex = results.begin();
         matchIndex != results.end(); ++matchIndex) {
      boost::shared_ptr<ROMol> match = lib.getMol(*matchIndex);
    }
  }

  {
    boost::shared_ptr<CachedTrustedSmilesMolHolder> molHolder =
        boost::make_shared<CachedTrustedSmilesMolHolder>();
    boost::shared_ptr<PatternHolder> patternHolder =
        boost::make_shared<PatternHolder>();

    SubstructLibrary lib(molHolder, patternHolder);
    lib.addMol(mol);
  }

  {
    boost::shared_ptr<CachedTrustedSmilesMolHolder> molHolder =
        boost::make_shared<CachedTrustedSmilesMolHolder>();
    boost::shared_ptr<PatternHolder> patternHolder =
        boost::make_shared<PatternHolder>();

    // the PatternHolder instance is able to make fingerprints.
    //  These, of course, can be read from a file.  For demonstration
    //   purposes we construct them here.
    const std::string trustedSmiles = "c1ccccc1";
    ROMol *m = SmilesToMol(trustedSmiles);
    const ExplicitBitVect *bitVector = patternHolder->makeFingerprint(*m);

    // The trusted smiles and bitVector can be read from any source.
    //  This is the fastest way to load a substruct library.
    molHolder->addSmiles(trustedSmiles);
    patternHolder->addFingerprint(*bitVector);
    SubstructLibrary lib(molHolder, patternHolder);
    delete m;
    delete bitVector;
  }

  delete q;
  delete m;
  BOOST_LOG(rdErrorLog) << "    Done (C++ doc tests)" << std::endl;
}

template <class Holder>
void ringTest(const std::string &name) {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing C++ ring query: " << name << std::endl;

  std::unique_ptr<ROMol> q(SmartsToMol("[C&R1]"));
  std::unique_ptr<ROMol> q2(SmartsToMol("C@C"));

  std::unique_ptr<ROMol> m(SmilesToMol("C1CCO[C@@](N)(O)1"));

  boost::shared_ptr<CachedTrustedSmilesMolHolder> molHolder =
      boost::make_shared<CachedTrustedSmilesMolHolder>();
  boost::shared_ptr<Holder> patternHolder = boost::make_shared<Holder>();

  SubstructLibrary lib(molHolder, patternHolder);
  lib.addMol(*m.get());
  std::vector<unsigned int> results = lib.getMatches(*q.get());
  TEST_ASSERT(results.size() == 1);
  results = lib.getMatches(*q2.get());
  TEST_ASSERT(results.size() == 1);

  BOOST_LOG(rdErrorLog) << "    Done (C++ ring query tests)" << std::endl;
}

void testAddPatterns() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "   Add Patterns " << std::endl;
  std::vector<std::string> pdb_ligands = {
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
      "CC1(C)CNc2cc(NC(=O)c3cccnc3NCc3ccncc3)ccc21"};

  boost::shared_ptr<CachedSmilesMolHolder> holder =
      boost::make_shared<CachedSmilesMolHolder>();

  for (auto s : pdb_ligands) {
    holder->addSmiles(s);
  }

  SubstructLibrary ssslib(holder);
  std::vector<int> num_threads = {1, 0};
  for (auto nthreads : num_threads) {
    SubstructLibrary ssslib_with_patterns(holder);
    SubstructLibrary ssslib_with_taut_patterns(holder);
    addPatterns(ssslib_with_patterns, nthreads);
    boost::shared_ptr<TautomerPatternHolder> patterns(
        new TautomerPatternHolder);
    addPatterns(ssslib_with_taut_patterns, patterns, nthreads);
    for (unsigned int i = 0; i < ssslib.size(); ++i) {
      TEST_ASSERT(ssslib.countMatches(*ssslib.getMol(i).get()) ==
                  ssslib_with_patterns.countMatches(*ssslib.getMol(i).get()));
      TEST_ASSERT(
          ssslib.countMatches(*ssslib.getMol(i).get()) ==
          ssslib_with_taut_patterns.countMatches(*ssslib.getMol(i).get()));
    }
  }
}

void testMaxResultsNumThreads() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "   Results do not depend on numThreads "
                        << std::endl;

  std::string fName = getenv("RDBASE");
  fName += "/Data/NCI/first_5K.smi";
  SmilesMolSupplier suppl(fName, "\t", 0, 1, false);
  auto *mols = new MolHolder();
  auto *fps = new PatternHolder();
  boost::shared_ptr<MolHolder> mols_ptr(mols);
  boost::shared_ptr<PatternHolder> fps_ptr(fps);

  SubstructLibrary ssslib(mols_ptr, fps_ptr);
  boost::logging::disable_logs("rdApp.error");
  while (!suppl.atEnd()) {
    ROMol *mol = nullptr;
    try {
      mol = suppl.next();
    } catch (...) {
      continue;
    }
    if (!mol) {
      continue;
    }
    ssslib.addMol(*mol);
    delete mol;
  }
  boost::logging::enable_logs("rdApp.error");
  std::vector<std::vector<unsigned int>> resVect;
  ROMOL_SPTR query(SmartsToMol("N"));
  TEST_ASSERT(query);
  for (auto numThreads : {1, 2, 4, 8}) {
    resVect.emplace_back(
        ssslib.getMatches(*query, true, false, false, numThreads));
  }
  for (auto it = resVect.begin() + 1; it != resVect.end(); ++it) {
    TEST_ASSERT(resVect.front().size() == it->size());
    for (size_t i = 0; i < resVect.front().size(); ++i) {
      TEST_ASSERT(resVect.front().at(i) == it->at(i));
    }
  }
  size_t results60 = resVect.front().size() * 0.6;
  size_t results99 = resVect.front().size() * 0.99;
  for (auto maxRes : {results60, results99}) {
    std::vector<std::vector<unsigned int>> resVectPartial;
    for (auto numThreads : {1, 2, 4, 8}) {
      resVectPartial.emplace_back(
          ssslib.getMatches(*query, true, false, false, numThreads, maxRes));
    }
    for (auto it = resVectPartial.begin(); it != resVectPartial.end(); ++it) {
      TEST_ASSERT(it->size() == maxRes);
      for (size_t i = 0; i < maxRes; ++i) {
        TEST_ASSERT(resVect.front().at(i) == it->at(i));
      }
    }
  }
}

void testMaxResultsAllSameNumThreads() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "   Results do not depend on numThreads (all same) "
                        << std::endl;

  auto *mols = new MolHolder();
  auto *fps = new PatternHolder();
  boost::shared_ptr<MolHolder> mols_ptr(mols);
  boost::shared_ptr<PatternHolder> fps_ptr(fps);

  SubstructLibrary ssslib(mols_ptr, fps_ptr);
  boost::logging::disable_logs("rdApp.error");
  auto mol = "N"_smiles;
  for (int i = 0; i < 999; ++i) {
    ssslib.addMol(*mol);
  }

  boost::logging::enable_logs("rdApp.error");
  std::vector<std::vector<unsigned int>> resVect;
  ROMOL_SPTR query(SmartsToMol("N"));
  TEST_ASSERT(query);
  for (auto numThreads : {1, 2, 4, 8}) {
    resVect.emplace_back(
        ssslib.getMatches(*query, true, false, false, numThreads));
    TEST_ASSERT(resVect.back().size() == 999);
  }
  for (auto it = resVect.begin() + 1; it != resVect.end(); ++it) {
    TEST_ASSERT(resVect.front().size() == it->size());
    for (size_t i = 0; i < resVect.front().size(); ++i) {
      TEST_ASSERT(resVect.front().at(i) == it->at(i));
    }
  }
  size_t results60 = resVect.front().size() * 0.6;
  size_t results99 = resVect.front().size() * 0.99;
  for (auto maxRes : {results60, results99}) {
    std::vector<std::vector<unsigned int>> resVectPartial;
    for (auto numThreads : {1, 2, 4, 8}) {
      resVectPartial.emplace_back(
          ssslib.getMatches(*query, true, false, false, numThreads, maxRes));
    }
    for (auto it = resVectPartial.begin(); it != resVectPartial.end(); ++it) {
      TEST_ASSERT(it->size() == maxRes);
      for (size_t i = 0; i < maxRes; ++i) {
        TEST_ASSERT(resVect.front().at(i) == it->at(i));
      }
    }
  }
}

template <class Holder>
void testPatternHolder(const std::string &name) {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "   testing " << name << std::endl;

  std::string fName = getenv("RDBASE");
  fName += "/Data/NCI/first_5K.smi";
  SmilesMolSupplier suppl(fName, "\t", 0, 1, false);
  boost::shared_ptr<CachedTrustedSmilesMolHolder> mols1(
      new CachedTrustedSmilesMolHolder());
  boost::shared_ptr<Holder> fps1(new Holder());
  SubstructLibrary ssslib1(mols1, fps1);
  boost::shared_ptr<CachedTrustedSmilesMolHolder> mols2(
      new CachedTrustedSmilesMolHolder());
  boost::shared_ptr<Holder> fps2(new Holder());
  SubstructLibrary ssslib2(mols2, fps2);

  boost::logging::disable_logs("rdApp.error");
  for (unsigned int i = 0; i < 1000; i += 10) {
    ROMol *mol = nullptr;
    try {
      mol = suppl[i];
    } catch (...) {
      continue;
    }
    if (!mol) {
      continue;
    }
    mols1->addSmiles(MolToSmiles(*mol));
    fps1->addFingerprint(fps1->makeFingerprint(*mol));
    ssslib2.addMol(*mol);
    delete mol;
  }
  boost::logging::enable_logs("rdApp.error");
  ROMOL_SPTR query(SmartsToMol("N"));
  TEST_ASSERT(query);
  {
    auto matches1 = ssslib1.getMatches(*query);
    std::sort(matches1.begin(), matches1.end());
    auto matches2 = ssslib2.getMatches(*query);
    std::sort(matches2.begin(), matches2.end());
    TEST_ASSERT(matches1.size() == matches2.size());
    for (size_t i = 0; i < matches1.size(); ++i) {
      TEST_ASSERT(matches1.at(i) == matches2.at(i));
    }
  }
#ifdef RDK_USE_BOOST_SERIALIZATION
  std::string pickle = ssslib1.Serialize();
  SubstructLibrary serialized;
  serialized.initFromString(pickle);
  TEST_ASSERT(serialized.size() == ssslib1.size());
  SubstructLibrary serializedLegacy;
  std::string pklName = getenv("RDBASE");
  TEST_ASSERT(!pklName.empty());
  pklName += "/Code/GraphMol/test_data/substructLibV1.pkl";
  std::ifstream pickle_istream(pklName.c_str(), std::ios_base::binary);
  serializedLegacy.initFromStream(pickle_istream);
  pickle_istream.close();
  TEST_ASSERT(serializedLegacy.size() == serialized.size());
  {
    auto matches1 = serializedLegacy.getMatches(*query);
    std::sort(matches1.begin(), matches1.end());
    auto matches2 = serialized.getMatches(*query);
    std::sort(matches2.begin(), matches2.end());
    TEST_ASSERT(matches1.size() == matches2.size());
    for (size_t i = 0; i < matches1.size(); ++i) {
      TEST_ASSERT(matches1.at(i) == matches2.at(i));
    }
  }
  for (size_t i = 0; i < 2; ++i) {
    auto serialized_pattern_holder =
        dynamic_cast<Holder *>(serialized.getFpHolder().get());
    TEST_ASSERT(serialized_pattern_holder);
    auto orig_pattern_holder =
        dynamic_cast<Holder *>(ssslib1.getFpHolder().get());
    TEST_ASSERT(orig_pattern_holder);
    TEST_ASSERT(serialized_pattern_holder->getNumBits() ==
                orig_pattern_holder->getNumBits());
    if (i) {
      break;
    }
    orig_pattern_holder->getNumBits() = 1024;
    pickle = ssslib1.Serialize();
    serialized.initFromString(pickle);
  }
#endif
}

void testSegFaultInHolder() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "   testSegFaultInHolder" << std::endl;

  boost::shared_ptr<CachedTrustedSmilesMolHolder> mols1(
      new CachedTrustedSmilesMolHolder());
  boost::shared_ptr<CachedSmilesMolHolder> mols2(new CachedSmilesMolHolder());
  for (int i = 0; i < 100; ++i) {
    if (i % 2 == 0) {
      mols1->addSmiles("dsafsdf");
      mols2->addSmiles("dsafsdf");
    } else {
      mols1->addSmiles("c1ccccc1");
      mols2->addSmiles("c1ccccc1");
    }
  }
  SubstructLibrary sss(mols1);
  SubstructLibrary sss2(mols2);
  ROMOL_SPTR query(SmartsToMol("c1ccccc1"));
  auto matches1 = sss.getMatches(*query);
  TEST_ASSERT(matches1.size() == 50);
  matches1 = sss2.getMatches(*query);
  TEST_ASSERT(matches1.size() == 50);

  // Check that we don't segfault when adding patterns
  addPatterns(sss, 2);
  addPatterns(sss2, 2);
}

void testTautomerQueries() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "   testTautomerQueries" << std::endl;

  boost::shared_ptr<CachedTrustedSmilesMolHolder> mols1(
      new CachedTrustedSmilesMolHolder());
  mols1->addSmiles("CN1C2=C(C(=O)Nc3ccccc3)C(=O)CCN2c2ccccc21");
  SubstructLibrary sss(mols1);
  auto query = "Cc1nc2ccccc2[nH]1"_smiles;
  // auto matches1 = sss.getMatches(*query);
  // TEST_ASSERT(matches1.size() == 0);
  std::unique_ptr<TautomerQuery> tq(TautomerQuery::fromMol(*query));
  auto matches2 = sss.getMatches(*tq);
  TEST_ASSERT(matches2.size() == 1);

  SubstructLibrary sss2(sss);
  addPatterns(sss, boost::make_shared<TautomerPatternHolder>());
  matches2 = sss.getMatches(*tq);
  TEST_ASSERT(matches2.size() == 1);

  // should work but throw logging errors
  addPatterns(sss2);
  matches2 = sss2.getMatches(*tq);
  TEST_ASSERT(matches2.size() == 1);
}

void github3881() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "  github3881 recursive smarts with rings "
                        << std::endl;
  boost::shared_ptr<CachedTrustedSmilesMolHolder> mols(
      new CachedTrustedSmilesMolHolder());
  mols->addSmiles("c1ccccc1S(=O)(=O)Cl");
  SubstructLibrary sss(mols);
  auto pat = "[$(S-!@[#6]):2](=O)(=O)(Cl)"_smarts;
  TEST_ASSERT(sss.getMatches(*pat).size() == 1);
}

int main() {
  RDLog::InitLogs();
#if 1
  test1();
  test2();
  test3();
  test4();
  docTest();
  ringTest<PatternHolder>("PatternHolder");
  ringTest<TautomerPatternHolder>("TautomerPatternHolder");
  testAddPatterns();
  testPatternHolder<PatternHolder>("PatternHolder");
  testPatternHolder<TautomerPatternHolder>("TautomerPatternHolder");
  testSegFaultInHolder();
#ifdef RDK_TEST_MULTITHREADED
  testMaxResultsNumThreads();
  testMaxResultsAllSameNumThreads();
  testTautomerQueries();
#endif
  github3881();
#endif
  return 0;
}
