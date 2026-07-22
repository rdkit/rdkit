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
#include <GraphMol/GaussianShape/GaussianShape.h>
#include <GraphMol/RascalMCES/RascalMCES.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>
#include <GraphMol/SynthonSpaceSearch/SearchResults.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearch_details.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
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
  const std::vector<size_t> maxRes{376417, 279307, 79833, 35604, 190, 46017};
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

unsigned int countFileLines(const std::string &filename) {
  std::ifstream ifs(filename.c_str());
  ifs.unsetf(std::ios_base::skipws);
  unsigned int numLines = std::count(std::istream_iterator<char>(ifs),
                                     std::istream_iterator<char>(), '\n');
  ifs.close();
  return numLines;
}

TEST_CASE("Shape More Small tests") {
  // A couple more full shape searches of small spaces.  Originally
  // They were in shape_search_catch_tests.cpp but made that too
  // slow for routine testing.
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string fullRoot(fName + "/Code/GraphMol/SynthonSpaceSearch/data/");
#if 0
  // These are the source files for the shape databases.  Useful to keep
  // around in case the databases ever need updating.
  std::vector<std::string> libNames{
      fullRoot + "triazole_space.txt",
      fullRoot + "urea_space.txt",
  };
  std::vector<std::string> enumLibNames{
      fullRoot + "triazole_space_enum.smi",
      fullRoot + "urea_space_enum.smi",
  };
#endif
  std::vector<std::string> dbNames{
      fullRoot + "triazole_space_shapes.spc",
      fullRoot + "urea_space_shapes.spc",
  };

  // The search of the enumerated libraries give 8 and 4 hits
  // respectively.
  std::vector<std::string> querySmis{
      "C[C@H]1CCN(c2nnc(CO)n2C2CCCC2)C1 |(3.88187,-1.9608,1.02401;3.40947,-0.556473,0.685633;3.49763,-0.278493,-0.787477;2.18424,0.313597,-1.21035;1.39217,0.480256,0.0110759;0.36454,1.42337,0.278193;0.656593,2.66561,0.766052;-0.477075,3.33073,0.923413;-1.48168,2.52297,0.540741;-2.93641,2.8664,0.551747;-3.33354,3.43671,-0.657732;-0.965924,1.32594,0.134935;-1.71735,0.224133,-0.333621;-1.11549,-0.643316,-1.37025;-2.26431,-1.65079,-1.54776;-2.58926,-1.96617,-0.0975393;-2.00229,-0.836476,0.740437;1.90383,-0.551243,0.906815),wD:1.0|",
      "C[C@@H]1C[C@H](NC(=O)NC2COC2)CN(C(=O)c2nccnc2F)C1 |(-0.346914,-0.986206,-4.28744;-0.686863,-0.0357247,-3.13265;0.429505,-0.1946,-2.14134;0.21099,0.659676,-0.907145;1.06526,0.0812473,0.104663;2.29297,0.75201,0.373712;2.5837,1.80373,-0.246936;3.23325,0.27478,1.33777;4.47647,0.99197,1.57602;4.94347,1.01294,2.99117;5.59284,-0.21541,2.82613;5.71049,0.107583,1.47157;-1.19052,0.721766,-0.417623;-2.14086,-0.0964663,-1.12904;-3.25312,-0.745252,-0.540367;-3.95877,-1.43507,-1.38825;-3.7227,-0.763533,0.803759;-4.9481,-0.204581,1.05395;-5.52492,-0.18107,2.24654;-4.86554,-0.748644,3.30451;-3.63585,-1.32731,3.13734;-3.08234,-1.32608,1.89052;-1.84839,-1.9059,1.69441;-2.02702,-0.329978,-2.57998),wD:1.0,wU:3.3|"};
  // The synthon search gives 1 hit for the urea space, where the
  // brute-force search gives 4 because the fragment similarities fall
  // below the threshold.  For example, comparing [2*]c1nccnc1F from
  // the query with synthon N#CCc(cncc1)c1[2*] (689988332-107515102)
  // when the dummy atoms are aligned, which they should be for a
  // good synthon match, the feature score is low because the nitrogen
  // acceptors don't align.  In the full molecule overlay, that is
  // compensated for by other things.
  std::vector<size_t> expNumHits{8, 1};
  std::vector<std::vector<double>> expScores{
      {0.9, 0.9, 0.9, 0.9, 0.8, 0.8, 0.8, 0.8}, {0.6}};
  std::vector<std::vector<std::string>> expNames{
      {"1-1;2-1;3-1;triazole-1", "1-1;2-2;3-1;triazole-1",
       "1-1;2-1;3-2;triazole-1", "1-1;2-2;3-2;triazole-1",
       "1-2;2-1;3-1;triazole-1", "1-2;2-1;3-2;triazole-1",
       "1-2;2-2;3-1;triazole-1", "1-2;2-2;3-2;triazole-1"},
      {"277310376-742385846;182115391-684092275;487354835-896308859;urea-3"}};
  ShapeBuildParams shapeBuildOptions;
  shapeBuildOptions.numConfs = 100;
  shapeBuildOptions.rmsThreshold = 0.5;
  shapeBuildOptions.numThreads = 1;
  shapeBuildOptions.shapeSimThreshold = 0.95;
  shapeBuildOptions.randomSeed = 0xdac;

  for (size_t i = 0; i < dbNames.size(); i++) {
    SynthonSpace synthonspace;
#if 0
    // In case the databases ever need updating.
    bool cancelled = false;
    synthonspace.readTextFile(libNames[i], cancelled);
    synthonspace.buildSynthonShapes(cancelled, shapeBuildOptions);
    synthonspace.writeDBFile(dbNames[i]);
#endif
    synthonspace.readDBFile(dbNames[i]);
    SynthonSpaceSearchParams params;
    params.similarityCutoff = 0.8;
    params.fragSimilarityAdjuster = 0.2;
    params.approxSimilarityAdjuster = 0.2;
    params.numConformers = shapeBuildOptions.numConfs;
    params.numThreads = shapeBuildOptions.numThreads;
    params.confRMSThreshold = shapeBuildOptions.rmsThreshold;
    params.timeOut = 0;
    params.randomSeed = shapeBuildOptions.randomSeed;
    params.bestHit = true;
    auto queryMol = v2::SmilesParse::MolFromSmiles(querySmis[i]);
    auto results = synthonspace.shapeSearch(*queryMol, params);
    unsigned int j = 0;
    CHECK(expNumHits[i] == results.getHitMolecules().size());
    for (const auto &mol : results.getHitMolecules()) {
      // Different machine/compiler combinations give slightly different results
      // for these tests which seems to be due to the conformer generator.
      CHECK(mol->getProp<double>("Similarity") > expScores[i][j]);
      CHECK(std::ranges::find(expNames[i], mol->getProp<std::string>(
                                               common_properties::_Name)) !=
            expNames[i].end());
      auto scores = GaussianShape::ScoreMolecule(*queryMol, *mol);
      CHECK_THAT(mol->getProp<double>("Similarity"),
                 Catch::Matchers::WithinAbs(scores[0], 0.001));
      ++j;
    }
#if 0
    // Leave this in for now, in case we need to check brute force search
    // in future.
    auto mols = loadLibrary(enumLibNames[i]);
    prepareMolecule(queryMol.get());
    // RDKit::SDWriter sdw2(enumOutputNames[i]);
    std::vector<float> matrix(12, 0.0);
    unsigned int numHits = 0;
    for (auto &[smiles, mol] : mols) {
      bool foundHit = false;
      for (unsigned int i = 0; i < queryMol->getNumConformers(); ++i) {
        for (unsigned int j = 0; j < mol->getNumConformers(); ++j) {
          auto [st, ct] = AlignMolecule(*queryMol, *mol, matrix, i, j);
          if (st + ct > params.similarityCutoff) {
            std::cout << mol->getProp<std::string>(common_properties::_Name)
                      << " hit at " << st << ", " << ct << " : " << st + ct
                      << " for " << i << ", " << j << std::endl;
            ++numHits;
            foundHit = true;
            // sdw2.write(*mol);
            break;
          }
        }
        if (foundHit) {
          break;
        }
      }
    }
#endif
  }
}

TEST_CASE("Shape Two piece query") {
  // Another, longer two piece query, originally in shape_search_catch_tests.cpp
  auto queryMol =
      "C1CCCC1.C[C@H]1CCNC1 |(-1.11549,-0.643316,-1.37025;-2.26431,-1.65079,-1.54776;-2.58926,-1.96617,-0.0975393;-2.00229,-0.836476,0.740437;-1.71735,0.224133,-0.333621;3.88187,-1.9608,1.02401;3.40947,-0.556473,0.685633;3.49763,-0.278493,-0.787477;2.18424,0.313597,-1.21035;1.39217,0.480256,0.0110759;1.90383,-0.551243,0.906815),wD:6.5|"_smiles;
  REQUIRE(queryMol);
  std::string fName(rdbase);
  std::string fullRoot(fName + "/Code/GraphMol/SynthonSpaceSearch/data/");
  std::string dbName = fullRoot + "triazole_space_shapes.spc";
  SynthonSpace synthonspace;
  synthonspace.readDBFile(dbName);
  SynthonSpaceSearchParams params;
  params.similarityCutoff = 0.7;
  params.fragSimilarityAdjuster = 0.2;
  params.numThreads = 1;
  params.confRMSThreshold = 0.5;
  params.timeOut = 0;
  params.randomSeed = 1;
  params.bestHit = true;
  params.shapeOverlayOptions.simAlpha = 0.95;
  params.shapeOverlayOptions.simBeta = 0.05;
  auto results = synthonspace.shapeSearch(*queryMol, params);
  CHECK(results.getHitMolecules().size() == 7);
  std::vector<double> expScores{0.801, 0.796, 0.778, 0.766,
                                0.744, 0.717, 0.711};
  for (unsigned int i = 0; i < results.getHitMolecules().size(); ++i) {
    auto &mol = results.getHitMolecules()[i];
    CHECK_THAT(mol->getProp<double>("Similarity"),
               Catch::Matchers::WithinAbs(expScores[i], 0.01));
  }
}

// There is only 1 conformer generator available, so just do it with
// different parameters.
static bool USED_generateConfs = false;  // To prove we went there.
std::unique_ptr<RWMol> generateConfs(const std::string &smiles,
                                     unsigned int numConformers) {
  USED_generateConfs = true;
  auto retMol = v2::SmilesParse::MolFromSmiles(smiles);
  MolOps::addHs(*retMol);
  auto dgParams = DGeomHelpers::ETKDGv3;
  dgParams.numThreads = 1;
  dgParams.pruneRmsThresh = -1;
  dgParams.randomSeed = 0xdac;
  auto cids =
      DGeomHelpers::EmbedMultipleConfs(*retMol, numConformers, dgParams);
  if (cids.empty() || !retMol->getNumConformers()) {
    retMol.reset();
  } else {
    MolOps::removeHs(*retMol);
  }
  return retMol;
}

TEST_CASE("Shape User-defined conformer generator") {
  // This is the urea space used elsewhere.
  std::string spaceText(
      R"(SMILES	synton_id	synton#	reaction_id
O=C(NC1COC1)[1*]	277310376-742385846	0	urea-3
#Oc(cc(cc1)NC([1*])=O)c1Cl	287123986-010598048	0	urea-3
O=C(Nc1c(CN[1*])cc[s]1)[2*]	584456271-623025187	1	urea-3
C[C@H](CC(C1)N[1*])CN1C([2*])=O	182115391-684092275	1	urea-3
FC(Oc(cc1)cnc1[2*])F	441848376-976230122	2	urea-3
Fc1nccnc1[2*]	487354835-896308859	2	urea-3
N#CCc(cncc1)c1[2*]	689988332-107515102	2	urea-3)");
  std::istringstream iss(spaceText);
  ShapeBuildParams shapeBuildParamsWith;
  shapeBuildParamsWith.numThreads = 1;
  shapeBuildParamsWith.numConfs = 10;
  shapeBuildParamsWith.userConformerGenerator = generateConfs;
  shapeBuildParamsWith.shapeSimThreshold = -1.0;
  bool cancelled = false;
  SynthonSpace space;
  space.readStream(iss, cancelled);
  space.buildSynthonShapes(cancelled, shapeBuildParamsWith);
  CHECK(USED_generateConfs);
  // There should be 10 conformers of each synthon except 182115391-684092275
  // which has an unspecified stereocentre so gets 20.
  auto react = space.getReaction(("urea-3"));
  for (auto sst : react->getSynthons()) {
    for (auto &[sn, s] : sst) {
      if (sn == "182115391-684092275") {
        CHECK(s->getShapes()->getShapes().getNumShapes() == 20);
      } else {
        CHECK(s->getShapes()->getShapes().getNumShapes() == 10);
      }
    }
  }
  // This is the same query as in Shape Small tests.
  auto queryMol =
      "C[C@@H]1C[C@H](NC(=O)NC2COC2)CN(C(=O)c2nccnc2F)C1 |(-0.346914,-0.986206,-4.28744;-0.686863,-0.0357247,-3.13265;0.429505,-0.1946,-2.14134;0.21099,0.659676,-0.907145;1.06526,0.0812473,0.104663;2.29297,0.75201,0.373712;2.5837,1.80373,-0.246936;3.23325,0.27478,1.33777;4.47647,0.99197,1.57602;4.94347,1.01294,2.99117;5.59284,-0.21541,2.82613;5.71049,0.107583,1.47157;-1.19052,0.721766,-0.417623;-2.14086,-0.0964663,-1.12904;-3.25312,-0.745252,-0.540367;-3.95877,-1.43507,-1.38825;-3.7227,-0.763533,0.803759;-4.9481,-0.204581,1.05395;-5.52492,-0.18107,2.24654;-4.86554,-0.748644,3.30451;-3.63585,-1.32731,3.13734;-3.08234,-1.32608,1.89052;-1.84839,-1.9059,1.69441;-2.02702,-0.329978,-2.57998),wD:1.0,wU:3.3|"_smiles;
  REQUIRE(queryMol);
  SynthonSpaceSearchParams spaceSearchParams;
  spaceSearchParams.userConformerGenerator = generateConfs;
  spaceSearchParams.useProgressBar = 0;
  spaceSearchParams.timeOut = 0;
  spaceSearchParams.similarityCutoff = 0.7;
  spaceSearchParams.bestHit = true;
  USED_generateConfs = false;
  auto results = space.shapeSearch(*queryMol, spaceSearchParams);
  CHECK(USED_generateConfs);
  CHECK(results.getHitMolecules().size() == 1);
  // This set of conformer parameters doesn't give a decent hit, but we
  // get something because of bestHit=true.
  for (const auto &mol : results.getHitMolecules()) {
    CHECK_THAT(mol->getProp<double>("Similarity"),
               Catch::Matchers::WithinAbs(0.728, 0.005));
    auto scores = GaussianShape::ScoreMolecule(*queryMol, *mol);
    CHECK_THAT(mol->getProp<double>("Similarity"),
               Catch::Matchers::WithinAbs(scores[0], 0.001));
  }
}

TEST_CASE("Shape Excluded volume") {
  // This is a piece of 4AJL, created by taking the 2 ligands from 4AJL and 4AJI
  // and dropping all atoms in 4AJL that were further than 5A from an atom
  // in either of the 2 ligands.  The RCSB PDB files of the 2 proteins are
  // aligned so the active sites are in the correct place.
  auto excVolMol =
      "C.C=O.CC.CC.CC(N)=O.CCC.CCC.CCC.CCC(=O)N[C@@H](C)C(=O)NCC=O.CCC(C)C.CCC(C)C.CCCC(=O)N[C@H](CN)CO.CCCCNC=O.CCCN[C@H](C=O)[C@@H](C)O.CCNC(=N)N.CNCC(N)=O.CNCCN.N.N.N=C(N)N.NC[C@H](CC(=O)O)NC=O.O.O.c1c[nH]cn1.CCC(C)C |(23.735,-5.551,4.831;10.558,-20.222,6.049;10.782,-19.593,7.082;25.756,-7.961,1.508;24.4,-8.595,1.992;7.488,-19.618,10.529;8.946,-19.661,10.098;16.988,-18.147,6.134;16.576,-18.196,7.564;15.432,-17.604,7.854;17.24,-18.815,8.408;15.229,-10.345,10.204;14.148,-9.909,9.204;13.16,-11.048,8.967;22.155,-3.456,-2.463;20.831,-2.784,-2.822;19.665,-3.765,-2.706;13.808,-17.366,17.137;12.337,-17.174,16.743;12.24,-16.466,15.411;18.52,-8.706,8.462;19.494,-9.564,7.605;19.097,-9.56,6.108;17.957,-9.849,5.775;20.037,-9.291,5.223;19.811,-9.168,3.789;20.952,-8.374,3.183;19.66,-10.496,3.038;20.161,-11.523,3.493;19.033,-10.416,1.856;18.853,-11.525,0.927;17.429,-11.999,0.744;16.5,-11.418,1.302;20.813,-11.741,-2.325;22.348,-11.427,-2.65;23.303,-11.616,-1.434;24.76,-11.817,-1.924;23.177,-10.46,-0.471;15.598,-14.535,14.266;14.867,-13.332,13.467;13.366,-13.102,13.804;13.147,-12.641,15.286;12.71,-12.136,12.79;18.012,-12.649,11.567;19.478,-12.385,11.172;20.187,-13.687,10.696;19.563,-14.252,9.404;18.599,-15.004,9.477;20.098,-13.86,8.229;19.648,-14.29,6.885;19.589,-15.803,6.809;18.634,-16.353,6.048;20.625,-13.781,5.83;20.759,-12.365,5.889;10.588,-11.836,3.288;11.37,-12.826,2.444;12.677,-13.19,3.171;13.48,-14.365,2.585;14.129,-13.976,1.342;15.412,-14.212,1.087;16.186,-14.727,1.907;7.421,-11.667,6.162;6.675,-12.99,6.581;6.654,-13.078,8.098;7.823,-13.191,8.731;7.885,-13.096,10.19;8.84,-11.926,10.495;9.77,-11.686,9.72;8.317,-14.401,10.877;7.872,-15.676,10.12;9.729,-14.357,11.106;13.973,-20.975,6.768;14.488,-21.058,8.204;13.543,-20.563,9.217;12.621,-21.314,9.805;11.859,-20.81,10.766;12.435,-22.565,9.423;14.659,-3.563,5.5;14.605,-4.881,5.293;13.35,-5.606,5.104;12.661,-5.673,6.457;11.536,-6.4,6.552;13.138,-5.069,7.429;21.154,-5.393,6.387;20.372,-4.647,5.612;19.099,-5.165,5.119;17.966,-4.53,5.886;17.138,-3.732,5.2;12.496,-8.081,8.665;12.367,-16.371,3.412;9.009,-17.463,14.046;9.374,-18.5,14.785;8.763,-18.749,15.936;10.373,-19.272,14.373;20.368,-1.798,-0.556;20.957,-1.177,0.463;20.595,-1.63,1.881;19.074,-1.555,2.174;18.319,-2.781,1.667;18.177,-2.926,0.433;17.852,-3.574,2.5;21.396,-0.863,2.83;21.647,-1.266,4.08;21.129,-2.256,4.588;26.022,-4.63,2.497;26.497,-10.537,-0.853;16.554,-20.899,12.518;16.264,-20.015,11.532;15.125,-19.344,11.93;14.762,-19.853,13.111;15.582,-20.814,13.498;27.715,-4.512,-2.12;27.605,-5.904,-2.016;28.301,-6.74,-2.887;29.098,-6.158,-3.87;28.228,-8.243,-2.741),wD:25.17,35.26,65.51,68.55,96.75,wU:40.30,49.38|"_smiles;
  // Query - ligands from 4AJI and 4AJ1, with proteins overlaid by RCSB.
  REQUIRE(excVolMol);
  GaussianShape::ShapeInputOptions excVolOptions;
  excVolOptions.useColors = false;
  auto excVolShape = std::make_unique<GaussianShape::ShapeInput>(*excVolMol, -1,
                                                                 excVolOptions);
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string fullRoot(fName + "/Code/GraphMol/SynthonSpaceSearch/data/");
  std::string dbName = fullRoot + "4ala_shapes.spc";
  SynthonSpace synthonSpace;
#if 0
  // This space is the ligand from 4al4, chopped up.  2-2 is the same as 2-1 but
  // with an extra phenyl group to ensure a big clash. 1-1 has an O instead of
  // an S which gave a bit less variability in the conformations.
  // Keeping it here in case it needs rebuilding.
  std::string spaceText(
      R"(SMILES	synton_id	synton#	reaction_id	release
Cc1nc2ccc(NC(=O)[1*])cc2o1	1-1	0	4al4
[1*]CCNC(=O)NCC[2*]	2-1	1	4al4
[1*]CCNC(=O)NC(c1ccccc1)C[2*]	2-2	1	4al4
[2*]Oc1ccc(cc1)CC(C(=O)O)C(=O)O	3-1	2	4al4
)");
  std::istringstream iss(spaceText);
  ShapeBuildParams shapeBuildParams;
  shapeBuildParams.numThreads = 1;
  shapeBuildParams.numConfs = 10;
  shapeBuildParams.shapeSimThreshold = -1.0;
  shapeBuildParams.randomSeed = 0xdac;
  bool cancelled = false;
  synthonSpace.readStream(iss, cancelled);
  synthonSpace.buildSynthonShapes(cancelled, shapeBuildParams);
  synthonSpace.writeDBFile(dbName);
#endif
  synthonSpace.readDBFile(dbName);

  auto comb_4aji_4aj1 =
      "CNc1nc2ccc(NC(C)=O)cc2s1.COc1ccc(CC(C(=O)O)C(=O)O)cc1OC |(23.956,-7.951,-4.055;24.081,-7.139,-2.841;23.068,-6.904,-1.967;23.237,-6.265,-0.838;22.03,-6.164,-0.135;21.811,-5.548,1.089;20.55,-5.531,1.66;19.465,-6.131,1.009;18.155,-6.121,1.583;17.318,-7.197,1.736;16.262,-7.032,2.809;17.43,-8.221,1.058;19.663,-6.743,-0.207;20.934,-6.758,-0.767;21.422,-7.467,-2.285;16.331,-12.235,6.625;15.378,-13.287,6.54;14.89,-13.803,7.711;15.669,-14.364,8.712;15.068,-14.88,9.87;13.688,-14.822,10.03;13.008,-15.402,11.226;12.32,-16.745,10.932;13.341,-17.807,10.61;14.399,-17.852,11.169;12.961,-18.68,9.696;11.481,-17.202,12.095;11.812,-18.168,12.768;10.397,-16.444,12.327;12.917,-14.275,9.014;13.494,-13.757,7.872;12.793,-13.19,6.852;11.466,-12.683,7.147)|"_smiles;

  SynthonSpaceSearchParams params;
  params.fragSimilarityAdjuster = 0.2;
  params.approxSimilarityAdjuster = 0.2;
  params.numConformers = 100;
  params.numThreads = 1;
  params.confRMSThreshold = 0.5;
  params.timeOut = 0;
  params.randomSeed = 0xdac;
  params.bestHit = true;
  params.maxHits = 10;
  params.similarityCutoff = 0.5;
  params.shapeOverlayOptions.simAlpha = 0.95;
  params.shapeOverlayOptions.simBeta = 0.05;
  params.excludedVolume = excVolShape.get();

  // This is one of those rare occasions where Mac and Linux give different
  // conformations even with the same parameters.  The results are similar
  // but ordered differently.
  std::vector<std::vector<double>> expVols{{71.2, 90.2}, {71.2, 95.1}};
  std::vector<std::vector<double>> expMeanVols{{3.0, 2.4}, {6.4, 2.3}};
  {
    params.possibleHitsFile = "poss_hits_1.txt";
    params.maxExcludedVolume = -1.0;
    auto results = synthonSpace.shapeSearch(*comb_4aji_4aj1, params);
    CHECK(results.getHitMolecules().size() == 2);
    unsigned int i = 0;
    for (const auto &mol : results.getHitMolecules()) {
      CHECK_THAT(mol->getProp<double>("ExcludedVolume"),
                 Catch::Matchers::WithinAbs(expVols[0][i], 0.1) ||
                     Catch::Matchers::WithinAbs(expVols[1][i], 0.1));
      CHECK_THAT(mol->getProp<double>("MeanExcludedVolume"),
                 Catch::Matchers::WithinAbs(expMeanVols[0][i], 0.1) ||
                     Catch::Matchers::WithinAbs(expMeanVols[1][i], 0.1));
      ++i;
    }
    CHECK(countFileLines("poss_hits_1.txt") == 2);
    // std::remove("poss_hits_1.txt");
  }

  {
    params.possibleHitsFile = "poss_hits_2.txt";
    params.maxExcludedVolume = 80.0;
    auto results = synthonSpace.shapeSearch(*comb_4aji_4aj1, params);
    CHECK(results.getHitMolecules().size() == 1);
    CHECK(results.getHitMolecules()[0]->getProp<std::string>(
              common_properties::_Name) == "1-1;2-1;3-1;4al4");
    CHECK_THAT(results.getHitMolecules()[0]->getProp<double>("ExcludedVolume"),
               Catch::Matchers::WithinAbs(expVols[0][0], 0.1) ||
                   Catch::Matchers::WithinAbs(expVols[1][1], 0.1));
    CHECK_THAT(
        results.getHitMolecules()[0]->getProp<double>("MeanExcludedVolume"),
        Catch::Matchers::WithinAbs(expMeanVols[0][0], 0.1) ||
            Catch::Matchers::WithinAbs(expMeanVols[1][1], 0.1));
    CHECK(countFileLines("poss_hits_2.txt") == 2);
    // std::remove("poss_hits_2.txt");
  }
}
