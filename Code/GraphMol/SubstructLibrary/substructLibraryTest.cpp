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

#include <GraphMol/Substruct/SubstructMatch.h>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>

using namespace RDKit;

namespace {
boost::dynamic_bitset<> runTest(SubstructLibrary &ssslib, const ROMol &pattern,
                                 int nThreads) {
  std::vector<unsigned int> libMatches = ssslib.getMatches(pattern, true, true, false, nThreads);
  boost::dynamic_bitset<> hasMatch(ssslib.size());
  BOOST_FOREACH (unsigned int idx, libMatches) { hasMatch[idx] = 1; }

  for (unsigned int i = 0; i < ssslib.size(); ++i) {
    MatchVectType match;
    bool matched = SubstructMatch(*ssslib.getMol(i), pattern, match);
    // std::cerr << MolToSmiles(*ssslib.getMol(i), true) << " " << hasMatch[i]
    //           << " " << matched << std::endl;
    TEST_ASSERT(hasMatch[i] == matched);
  }
  return hasMatch;
};

void runTest(SubstructLibrary &ssslib,
             const ROMol &pattern,
             int nThreads,
             const boost::dynamic_bitset<> &hasMatch
             ) {
  std::vector<unsigned int> libMatches = ssslib.getMatches(pattern, true, true, false, nThreads);
  boost::dynamic_bitset<> hasMatch2(ssslib.size());
  BOOST_FOREACH (unsigned int idx, libMatches) { hasMatch2[idx] = 1; }
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
    if (!mol) continue;
    ssslib.addMol(*mol);
    delete mol;
  }

  std::vector<SubstructLibrary*> libs;
  libs.push_back(&ssslib);

#ifdef RDK_USE_BOOST_SERIALIZATION  
  std::string pickle = ssslib.Serialize();
  SubstructLibrary serialized;
  serialized.initFromString(pickle);
  TEST_ASSERT(serialized.size() == ssslib.size());
  libs.push_back(&serialized);
#endif

  boost::dynamic_bitset<> hasMatch;
  
  int i=0;
  for(auto lib: libs) {
    ROMol *query = SmartsToMol("[#6;$([#6]([#6])[!#6])]");
    if (i == 0)
      hasMatch = runTest(*lib, *query, 1);
    else
      runTest(*lib, *query, 1, hasMatch);
    
#ifdef RDK_TEST_MULTITHREADED
    runTest(*lib, *query, -1, hasMatch);
#endif
    delete query;
    ++i;
  }

  i = 0;
  for(auto lib: libs) {
    ROMol *query = SmartsToMol("[$([O,S]-[!$(*=O)])]");
    if (i == 0)
      hasMatch = runTest(*lib, *query, 1);
    else
      runTest(*lib, *query, 1, hasMatch);

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
    if (!mol) continue;
    ssslib.addMol(*mol);
    delete mol;
  }

  std::vector<SubstructLibrary*> libs;
  libs.push_back(&ssslib);

#ifdef RDK_USE_BOOST_SERIALIZATION  
  std::string pickle = ssslib.Serialize();
  SubstructLibrary serialized;
  serialized.initFromString(pickle);
  TEST_ASSERT(serialized.size() == ssslib.size());

  // check to see if we are still the right base type
  MolHolderBase *_holder = serialized.getMolHolder().get();
  TEST_ASSERT(_holder != nullptr);
  TEST_ASSERT(dynamic_cast<MolHolder*>(_holder) != nullptr);
  try { serialized.getFingerprints(); }
  catch(...) { TEST_ASSERT(0); }
  
  libs.push_back(&serialized);
#endif

  for(auto lib: libs) {  
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

  std::vector<SubstructLibrary*> libs;
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
  TEST_ASSERT(dynamic_cast<MolHolder*>(_holder) != nullptr);  
#endif


  for(auto lib: libs) {  
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

  std::vector<SubstructLibrary*> libs;
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
  TEST_ASSERT(dynamic_cast<CachedSmilesMolHolder*>(_holder) != nullptr);
#endif

  for(auto lib: libs) {
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

int main() {
  RDLog::InitLogs();
#if 1
  test1();
  test2();
  test3();
  test4();
  docTest();
#endif
  return 0;
}
