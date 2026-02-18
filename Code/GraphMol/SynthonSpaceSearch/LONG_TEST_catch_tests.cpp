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
  std::cout << results.getHitMolecules().front()->getProp<double>("Similarity")
            << " to "
            << results.getHitMolecules().back()->getProp<double>("Similarity")
            << std::endl;
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
  CHECK(results.getHitMolecules().size() == 1000);
  CHECK(results.getMaxNumResults() == 3319);

  // A tighter adjuster misses more hits.
  params.approxSimilarityAdjuster = 0.01;
  results = synthonspace.fingerprintSearch(*queryMol, *fpGen, params);
  CHECK(results.getHitMolecules().size() == 855);

  // This is the actual number of hits achievable.
  params.approxSimilarityAdjuster = 0.25;
  params.maxHits = -1;
  results = synthonspace.fingerprintSearch(*queryMol, *fpGen, params);
  CHECK(results.getHitMolecules().size() == 1007);
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

TEST_CASE("Shape Long Unspecified Stereo") {
  auto m1 = "C[C@H](Cl)CCOOC(Cl)F"_smiles;
  REQUIRE(m1);
  CHECK(details::hasUnspecifiedStereo(*m1) == true);
  CHECK(details::countChiralAtoms(*m1) == 2);

  auto m2 = "C[C@H](Cl)CCOO[C@@H](Cl)F"_smiles;
  REQUIRE(m2);
  CHECK(details::hasUnspecifiedStereo(*m2) == false);
  CHECK(details::countChiralAtoms(*m2) == 2);

  auto m3 = "C[C@H](Cl)CCOO[C@@](Cl)(F)CC=CC"_smiles;
  REQUIRE(m3);
  CHECK(details::hasUnspecifiedStereo(*m3) == true);

  auto m4 = R"(C[C@H](Cl)CCOO[C@@](Cl)(F)C\C=C/C)"_smiles;
  REQUIRE(m4);
  CHECK(details::hasUnspecifiedStereo(*m4) == false);

  SynthonSpace space;
  std::istringstream iss(R"(SMILES	synton_id	synton#	reaction_id
O=C(c1coc2cc(Cl)ccc12)[1*:1]	A	2	r1	3
OC(CC1CCCN1[1*:1])c1ccco1	B	1	r1	3
Cc1cnc([1*:1])nc1C	C	2	r2	3
Cn1cc(N2CCCC(N[1*:1])C2)cn1	D	1	r2	3
Cc1cc(N[1*:1])cc(C(F)(F)F)c1	E	1	r3	3
Cc1nc(-c2ccccn2)cc([1*:1])n1	F	2	r3	3
CC(C)(C)C1(C)CN([1*:1])CCO1	G	1	r4	3
CCCC1CN([2*:2])CCO1	H	3	r4	3
COC1C2CCCC2C1N[2*:2]	I	3	r4	3
COCc1nc([1*:1])cc([2*:2])n1	J	2	r4	3
FC1(F)CCC(N[1*:1])CC1	K	1	r4	3)");
  bool cancelled = false;
  space.readStream(iss, cancelled);

  ShapeBuildParams shapeOptions;
  shapeOptions.randomSeed = 1;
  space.buildSynthonShapes(cancelled, shapeOptions);

  SynthonSpaceSearchParams params;
  params.similarityCutoff = 1.6;
  params.enumerateUnspecifiedStereo = false;

  // This should bale with no results because there's unspecified
  // stereochem.
  auto results = space.shapeSearch(*m1, params);
  CHECK(results.getHitMolecules().empty());

  // This is one of the molecules in the library, so should always
  // be a hit.
  auto m5 = R"(Cc1cnc(NC2CCCN(c3cnn(C)c3)C2)nc1C)"_smiles;
  REQUIRE(m5);
  CHECK(details::hasUnspecifiedStereo(*m5) == true);

  params.enumerateUnspecifiedStereo = true;
  params.randomSeed = 1;
  results = space.shapeSearch(*m5, params);
  REQUIRE(results.getHitMolecules().size() == 1);
  auto &hitMol1 = results.getHitMolecules().front();
  double firstSim = hitMol1->getProp<double>("Similarity");
  params.bestHit = true;
  results = space.shapeSearch(*m5, params);
  REQUIRE(results.getHitMolecules().size() == 1);
  auto &hitMol2 = results.getHitMolecules().front();
  double bestSim = hitMol2->getProp<double>("Similarity");
  CHECK(bestSim > firstSim);
}

TEST_CASE("Shape Hits back onto query") {
  // Make sure the hits are properly translated to the reference
  // frame of the query.
  SynthonSpace synthonspace;
  REQUIRE(rdbase);
  std::string fName(rdbase);
  std::string libName =
      fName +
      "/Code/GraphMol/SynthonSpaceSearch/data/hits_back_onto_query_test_shapes.spc";
  synthonspace.readDBFile(libName);
  auto tagrisso_pdb_core =
      "c1cc(Nc2nccc(c3cn(C)c4ccccc34)n2)ccc1 |(-30.966,18.467,-10.003;-29.741,18.8,-10.881;-29.776,18.58,-12.402;-28.626,18.878,-13.264;-27.858,20.11,-13.139;-26.809,20.446,-14.135;-26.039,21.676,-14.006;-26.301,22.606,-12.864;-27.356,22.266,-11.866;-27.643,23.19,-10.674;-26.776,24.159,-10.172;-27.396,24.761,-9.099;-26.842,25.83,-8.286;-28.633,24.178,-8.929;-29.782,24.445,-7.884;-31.052,23.635,-7.939;-31.218,22.587,-8.984;-30.11,22.344,-9.979;-28.784,23.198,-9.912;-28.114,21.037,-12.005;-31.044,18.019,-13.045;-32.253,17.694,-12.176;-32.227,17.912,-10.676)|"_smiles;
  SynthonSpaceSearchParams params;
  params.maxHits = -1;
  params.numThreads = -1;
  params.similarityCutoff = 1.0;
  params.numConformers = 100;
  params.confRMSThreshold = 1.0;
  params.timeOut = 0;
  params.randomSeed = 0xdac;

  RDGeom::Point3D tag_centre;
  for (const auto atom : tagrisso_pdb_core->atoms()) {
    tag_centre += tagrisso_pdb_core->getConformer().getAtomPos(atom->getIdx());
  }
  tag_centre /= tagrisso_pdb_core->getNumAtoms();
  auto results = synthonspace.shapeSearch(*tagrisso_pdb_core, params);
  REQUIRE(!results.getHitMolecules().empty());
  for (const auto &m : results.getHitMolecules()) {
    RDGeom::Point3D hit_centre;
    for (const auto atom : m->atoms()) {
      hit_centre += m->getConformer().getAtomPos(atom->getIdx());
    }
    hit_centre /= m->getNumAtoms();
    CHECK((hit_centre - tag_centre).length() < 2.0);
  }
}
