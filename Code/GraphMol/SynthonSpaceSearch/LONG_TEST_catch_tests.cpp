//
//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>

#include <GraphMol/SubstructLibrary/SubstructLibrary.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <GraphMol/Fingerprints/RDKitFPGenerator.h>
#include <GraphMol/RascalMCES/RascalMCES.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>
#include <GraphMol/SynthonSpaceSearch/SearchResults.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <boost/multiprecision/number.hpp>
#include <boost/parameter/aux_/pp_impl/match.hpp>

#include <catch2/catch_all.hpp>

using namespace RDKit;
using namespace RDKit::SynthonSpaceSearch;
using namespace RDKit::SynthonSpaceSearch::details;
using namespace RDKit::RascalMCES;

const char *rdbase = getenv("RDBASE");

int num5567MReads = 0;
int num5567RReads = 0;
int num5567NReads = 0;

// The Synthons_5567.csv is in the test data.  It is used multiple times
// so convert to the binary if that isn't already there.  This makes the
// tests faster and also tests the conversion process.
void convert5567FileR() {
  num5567RReads++;
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.csv";
  std::string binName1 =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567_rdkit.spc";
  if (!std::filesystem::exists(binName1)) {
    bool cancelled = false;
    std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGen(
        RDKitFP::getRDKitFPGenerator<std::uint64_t>());
    convertTextToDBFile(libName, binName1, cancelled, fpGen.get());
  }
}

void convert5567FileM() {
  num5567MReads++;
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.csv";
  std::string binName2 =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567_morgan.spc";
  if (!std::filesystem::exists(binName2)) {
    bool cancelled = false;
    std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGen(
        MorganFingerprint::getMorganGenerator<std::uint64_t>(2));
    convertTextToDBFile(libName, binName2, cancelled, fpGen.get());
  }
}

void convert5567FileN() {
  num5567NReads++;
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.csv";
  std::string binName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.spc";
  if (!std::filesystem::exists(binName)) {
    bool cancelled = false;
    convertTextToDBFile(libName, binName, cancelled);
  }
}

void tidy5567Binary() {
  std::string fName(rdbase);
  std::string binName1 =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567_rdkit.spc";
  std::string binName2 =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567_morgan.spc";
  if (num5567RReads == 1 && std::filesystem::exists(binName1)) {
    std::remove(binName1.c_str());
  }
  if (num5567MReads == 3 && std::filesystem::exists(binName2)) {
    std::remove(binName2.c_str());
  }
  std::string binName3 =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.spc";
  if (num5567NReads == 4 && std::filesystem::exists(binName3)) {
    std::cout << "removing " << binName3 << std::endl;
    std::remove(binName3.c_str());
  }
}

TEST_CASE("S Biggy") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  convert5567FileN();
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.spc";
  SynthonSpace synthonspace;
  synthonspace.readDBFile(libName);
  const std::vector<std::string> smis{"c1ccccc1C(=O)N1CCCC1",
                                      "c1ccccc1NC(=O)C1CCN1",
                                      "c12ccccc1c(N)nc(N)n2",
                                      "c12ccc(C)cc1[nH]nc2C(=O)NCc1cncs1",
                                      "c1nncn1",
                                      "C(=O)NC(CC)C(=O)N(CC)C"};
  const std::vector<size_t> numRes{6785, 4544, 48892, 1, 29147, 5651};
  const std::vector<size_t> maxRes{6785, 4544, 48893, 1, 29312, 5869};
  SubstructMatchParameters matchParams;
  SynthonSpaceSearchParams params;
  params.maxHits = -1;
  for (auto numThreads : std::vector<int>{1, 2, -1}) {
    params.numThreads = numThreads;
    for (size_t i = 0; i < smis.size(); ++i) {
      auto queryMol = v2::SmilesParse::MolFromSmarts(smis[i]);
      auto results =
          synthonspace.substructureSearch(*queryMol, matchParams, params);
      CHECK(results.getHitMolecules().size() == numRes[i]);
      CHECK(results.getMaxNumResults() == maxRes[i]);
    }
  }
  tidy5567Binary();
}

TEST_CASE("S Random Hits") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.spc";
  SynthonSpace synthonspace;
  convert5567FileN();
  synthonspace.readDBFile(libName);

  auto queryMol = "c1ccccc1C(=O)N1CCCC1"_smiles;
  SubstructMatchParameters matchParams;
  SynthonSpaceSearchParams params;
  params.maxHits = 100;
  params.randomSample = true;
  params.randomSeed = 1;
  params.numThreads = -1;
  auto results =
      synthonspace.substructureSearch(*queryMol, matchParams, params);
  std::map<std::string, int> libCounts;
  for (const auto &m : results.getHitMolecules()) {
    auto molName = m->getProp<std::string>(common_properties::_Name);
    std::string lib(molName.substr(molName.length() - 2));
    if (const auto &c = libCounts.find(lib); c == libCounts.end()) {
      libCounts.insert(std::make_pair(lib, 1));
    } else {
      c->second++;
    }
  }
  CHECK(results.getHitMolecules().size() == 100);
  // std::shuffle gives different results on macOS, Linux and Windows, so
  // just check we have some hits from each of a1, a6 and a7.
  CHECK(libCounts.find("a1") != libCounts.end());
  CHECK(libCounts.find("a6") != libCounts.end());
  CHECK(libCounts.find("a7") != libCounts.end());
  tidy5567Binary();
}

TEST_CASE("S Later hits") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  // Test use of params.hitStart
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.spc";
  SynthonSpace synthonspace;
  convert5567FileN();
  synthonspace.readDBFile(libName);

  auto queryMol = "c1ccccc1C(=O)N1CCCC1"_smiles;
  SubstructMatchParameters matchParams;
  SynthonSpaceSearchParams params;
  params.maxHits = 200;
  params.numThreads = -1;
  auto results =
      synthonspace.substructureSearch(*queryMol, matchParams, params);
  std::vector<std::string> hitNames1;
  for (const auto &m : results.getHitMolecules()) {
    hitNames1.push_back(m->getProp<std::string>(common_properties::_Name));
  }

  params.maxHits = 100;
  params.hitStart = 100;
  results = synthonspace.substructureSearch(*queryMol, matchParams, params);
  std::vector<std::string> hitNames2;
  for (const auto &m : results.getHitMolecules()) {
    hitNames2.push_back(m->getProp<std::string>(common_properties::_Name));
  }
  CHECK(hitNames1.size() == 200);
  CHECK(hitNames2.size() == 100);
  for (int i = 0; i < 100; ++i) {
    CHECK(hitNames1[100 + i] == hitNames2[i]);
  }

  params.hitStart = 6780;
  results = synthonspace.substructureSearch(*queryMol, matchParams, params);
  CHECK(results.getHitMolecules().size() == 5);

  params.hitStart = 7000;
  results = synthonspace.substructureSearch(*queryMol, matchParams, params);
  CHECK(results.getHitMolecules().empty());
  tidy5567Binary();
}

TEST_CASE("S Complex query") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  // Just to demonstrate that a complex query works.
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.spc";
  SynthonSpace synthonspace;
  convert5567FileN();
  synthonspace.readDBFile(libName);

  auto queryMol = v2::SmilesParse::MolFromSmarts(
      "[$(c1ccccc1),$(c1ccncc1),$(c1cnccc1)]C(=O)N1[C&!$(CC(=O))]CCC1");
  REQUIRE(queryMol);
  SubstructMatchParameters matchParams;
  SynthonSpaceSearchParams params;
  params.maxHits = -1;
  params.numThreads = -1;
  auto results =
      synthonspace.substructureSearch(*queryMol, matchParams, params);
  CHECK(results.getHitMolecules().size() == 7649);
  // The screenout is poor for a complex query, so a lot of things
  // will be identified as possible that aren't.
  CHECK(results.getMaxNumResults() == 32814);
  tidy5567Binary();
}

TEST_CASE("FP Biggy") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  convert5567FileM();
  std::string binName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567_morgan.spc";
  SynthonSpace synthonspace;
  synthonspace.readDBFile(binName);

  std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGen(
      MorganFingerprint::getMorganGenerator<std::uint64_t>(2));

  std::map<std::string, std::unique_ptr<RWMol>> mols;

  const std::vector<std::string> smis{
      "c1ccccc1C(=O)N1CCCC1", "c1ccccc1NC(=O)C1CCN1",
      "c12ccccc1c(N)nc(N)n2", "c12ccc(C)cc1[nH]nc2C(=O)NCc1cncs1",
      "c1n[nH]cn1",           "C(=O)NC(CC)C(=O)N(CC)C"};
  const std::vector<size_t> numRes{46, 2, 0, 123, 0, 0};
  const std::vector<size_t> maxRes{2408, 197, 0, 833, 0, 4};
  SynthonSpaceSearchParams params;
  params.approxSimilarityAdjuster = 0.2;
  params.maxHits = -1;
  params.numThreads = -1;
  for (size_t i = 0; i < smis.size(); ++i) {
    auto queryMol = v2::SmilesParse::MolFromSmiles(smis[i]);
    auto results = synthonspace.fingerprintSearch(*queryMol, *fpGen, params);
    CHECK(results.getHitMolecules().size() == numRes[i]);
    CHECK(results.getMaxNumResults() == maxRes[i]);
  }
  tidy5567Binary();
}

TEST_CASE("FP Random Hits") {
  REQUIRE(rdbase);
  convert5567FileM();
  std::string fName(rdbase);
  std::string binName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567_morgan.spc";
  SynthonSpace synthonspace;
  synthonspace.readDBFile(binName);

  auto queryMol = "c12ccc(C)cc1[nH]nc2C(=O)NCc1cncs1"_smiles;
  SynthonSpaceSearchParams params;
  params.maxHits = 100;
  params.randomSample = true;
  params.randomSeed = 1;
  params.numThreads = -1;
  params.approxSimilarityAdjuster = 0.1;
  std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGen(
      MorganFingerprint::getMorganGenerator<std::uint64_t>(2));
  auto results = synthonspace.fingerprintSearch(*queryMol, *fpGen, params);
  std::map<std::string, int> libCounts;
  for (const auto &m : results.getHitMolecules()) {
    std::string molName = m->getProp<std::string>(common_properties::_Name);
    std::string lib(molName.substr(molName.length() - 2));
    if (const auto &c = libCounts.find(lib); c == libCounts.end()) {
      libCounts.insert(std::make_pair(lib, 1));
    } else {
      c->second++;
    }
  }
  CHECK(results.getHitMolecules().size() == 100);
  std::map<std::string, int> expCounts{{"a1", 100}};
  CHECK(expCounts == libCounts);
  tidy5567Binary();
}

TEST_CASE("Timeout") {
  REQUIRE(rdbase);
  convert5567FileM();
  std::string fName(rdbase);
  std::string binName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567_morgan.spc";
  SynthonSpace synthonspace;
  synthonspace.readDBFile(binName);
  SynthonSpaceSearchParams params;
  params.maxHits = -1;
  params.similarityCutoff = 0.3;
  params.fragSimilarityAdjuster = 0.3;
  params.timeOut = 2;
  params.numThreads = -1;
  std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGen(
      MorganFingerprint::getMorganGenerator<std::uint64_t>(2));

  auto queryMol = "c12ccc(C)cc1[nH]nc2C(=O)NCc1cncs1"_smiles;
  auto results = synthonspace.fingerprintSearch(*queryMol, *fpGen, params);
  CHECK(results.getTimedOut());

  // Make sure no timeout also works, but only on a short search.
  params.maxHits = 100;
  params.similarityCutoff = 0.3;
  params.fragSimilarityAdjuster = 0.2;
  params.timeOut = 0;
  auto results1 = synthonspace.fingerprintSearch(*queryMol, *fpGen, params);
  CHECK(!results1.getTimedOut());
  tidy5567Binary();
}

TEST_CASE("FP Approx Similarity") {
  REQUIRE(rdbase);
  convert5567FileR();
  std::string fName(rdbase);
  std::string binName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567_rdkit.spc";
  SynthonSpace synthonspace;
  synthonspace.readDBFile(binName);
  SynthonSpaceSearchParams params;
  params.randomSeed = 1;
  params.similarityCutoff = 0.5;
  params.timeOut = 0;
  params.maxHits = 1000;
  params.numThreads = -1;

  std::unique_ptr<FingerprintGenerator<std::uint64_t>> fpGen(
      RDKitFP::getRDKitFPGenerator<std::uint64_t>());
  auto queryMol = "c12ccc(C)cc1[nH]nc2C(=O)NCc1cncs1"_smiles;

  // With RDKit fingerprints, 0.05 gives a reasonable compromise
  // between speed and hits missed.
  params.approxSimilarityAdjuster = 0.05;
  auto results = synthonspace.fingerprintSearch(*queryMol, *fpGen, params);
  CHECK(results.getHitMolecules().size() == 546);
  CHECK(results.getMaxNumResults() == 1478);

  // A tighter adjuster misses more hits.
  params.approxSimilarityAdjuster = 0.01;
  results = synthonspace.fingerprintSearch(*queryMol, *fpGen, params);
  CHECK(results.getHitMolecules().size() == 187);

  // This is the actual number of hits achievable.
  params.approxSimilarityAdjuster = 0.25;
  params.maxHits = -1;
  results = synthonspace.fingerprintSearch(*queryMol, *fpGen, params);
  CHECK(results.getHitMolecules().size() == 981);
  tidy5567Binary();
}

TEST_CASE("Rascal Biggy") {
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName + "/Code/GraphMol/SynthonSpaceSearch/data/Syntons_5567.csv";
  SynthonSpace synthonspace;
  bool cancelled = false;
  synthonspace.readTextFile(libName, cancelled);

  const std::vector<std::string> smis{
      "c1ccccc1C(=O)N1CCCC1", "c1ccccc1NC(=O)C1CCN1",
      "c12ccccc1c(N)nc(N)n2", "c12ccc(C)cc1[nH]nc2C(=O)NCc1cncs1",
      "c1n[nH]cn1",           "C(=O)NC(CC)C(=O)N(CC)C"};
  const std::vector<size_t> numRes{254, 89, 2, 34, 0, 14};
  const std::vector<size_t> maxRes{376110, 278747, 79833, 34817, 190, 45932};
  SynthonSpaceSearchParams params;
  params.maxHits = -1;
  params.numThreads = -1;
  RascalOptions rascalOptions;

  for (size_t i = 0; i < smis.size(); ++i) {
    auto queryMol = v2::SmilesParse::MolFromSmiles(smis[i]);
    auto results = synthonspace.rascalSearch(*queryMol, rascalOptions, params);
    CHECK(results.getHitMolecules().size() == numRes[i]);
    CHECK(results.getMaxNumResults() == maxRes[i]);
  }
}
